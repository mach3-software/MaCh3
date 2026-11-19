// CudaReweightDevice: CpuReweight's cached pass and gradient on an NVIDIA GPU, see CudaReweightDevice.h
#include "Samples/CudaReweightDevice.h"

#include <cuda_runtime.h>
#include <cub/cub.cuh>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void Check(const cudaError_t Status, const char* What) {
  if (Status != cudaSuccess) throw std::runtime_error(std::string("CUDA ") + What + ": " + cudaGetErrorString(Status));
}

/// @brief A device array of T that frees itself
template <typename T>
struct DeviceArray {
  T* Ptr = nullptr;
  std::size_t N = 0;
  DeviceArray() = default;
  DeviceArray(const DeviceArray&) = delete;
  DeviceArray& operator=(const DeviceArray&) = delete;
  ~DeviceArray() { if (Ptr) cudaFree(Ptr); }
  void Allocate(const std::size_t Count) {
    if (Ptr) cudaFree(Ptr);
    Ptr = nullptr;
    N = Count;
    Check(cudaMalloc(&Ptr, std::max<std::size_t>(1, Count) * sizeof(T)), "malloc");
  }
  void Upload(const std::vector<T>& Values) {
    Allocate(Values.size());
    if (!Values.empty()) Check(cudaMemcpy(Ptr, Values.data(), Values.size() * sizeof(T), cudaMemcpyHostToDevice), "upload");
  }
  std::size_t Bytes() const { return N * sizeof(T); }
};

/// @brief Pinned host memory of T that frees itself
template <typename T>
struct PinnedArray {
  T* Ptr = nullptr;
  ~PinnedArray() { if (Ptr) cudaFreeHost(Ptr); }
  void Allocate(const std::size_t Count) {
    if (Ptr) cudaFreeHost(Ptr);
    Check(cudaMallocHost(&Ptr, std::max<std::size_t>(1, Count) * sizeof(T)), "host malloc");
  }
};

constexpr int kBlock = 256;
unsigned int Blocks(const std::size_t N) { return static_cast<unsigned int>((N + kBlock - 1) / kBlock); }

// The same four steps as MetalReweight's kernels and CpuReweight's loops, in double precision

__global__ void SigProducts(const double* __restrict__ Pool, const int* __restrict__ SigStart, const int* __restrict__ SigPool,
                            double* __restrict__ SigProd, const int NSig) {
  const int s = blockIdx.x * blockDim.x + threadIdx.x;
  if (s >= NSig) return;
  double p = 1.;
  for (int k = SigStart[s]; k < SigStart[s + 1]; ++k) p *= Pool[SigPool[k]];
  SigProd[s] = p;
}

// One block per bin, its threads taking every kBlock-th group, then a fixed order reduction.
// A bin holds thousands of groups (more with binned oscillation, which splits them by energy bin), so
// one thread per bin left the GPU almost idle: 3.5 ms per pass on an H100 for the DUNE FD oscillation fit
__global__ void BinSums(const int* __restrict__ BinStart, const int* __restrict__ GroupSig, const double* __restrict__ GroupConst,
                        const double* __restrict__ GroupConst2, const double* __restrict__ SigProd, double* __restrict__ MC,
                        double* __restrict__ W2, const int NBins) {
  __shared__ double Partial[kBlock];
  __shared__ double Partial2[kBlock];
  const int b = blockIdx.x;
  const int t = threadIdx.x;
  double Sum = 0., Sum2 = 0.;
  if (b < NBins) {
    for (int g = BinStart[b] + t; g < BinStart[b + 1]; g += kBlock) {
      const double P = SigProd[GroupSig[g]];
      const double W = GroupConst[g] * P;
      if (W <= 0.) continue;
      Sum += W;
      Sum2 += GroupConst2[g] * P * P;
    }
  }
  Partial[t] = Sum;
  Partial2[t] = Sum2;
  __syncthreads();
  for (int Half = kBlock / 2; Half > 0; Half /= 2) {
    if (t < Half) {
      Partial[t] += Partial[t + Half];
      Partial2[t] += Partial2[t + Half];
    }
    __syncthreads();
  }
  if (t == 0 && b < NBins) {
    MC[b] = Partial[0];
    W2[b] = Partial2[0];
  }
}

__global__ void SigShare(const int* __restrict__ SigGroupStart, const int* __restrict__ SigGroup, const int* __restrict__ GroupBin,
                         const double* __restrict__ GroupConst, const double* __restrict__ SigProd, const double* __restrict__ BinFactor,
                         double* __restrict__ Share, const int NSig) {
  const int s = blockIdx.x * blockDim.x + threadIdx.x;
  if (s >= NSig) return;
  double a = 0.;
  for (int k = SigGroupStart[s]; k < SigGroupStart[s + 1]; ++k) {
    const int g = SigGroup[k];
    const double W = GroupConst[g] * SigProd[s];
    if (W > 0.) a += BinFactor[GroupBin[g]] * W;
  }
  Share[s] = a;
}

/// Terms per block of ParamGrad. Each parameter's terms are cut into chunks of this many, so a
/// parameter every signature depends on (an oscillation parameter, say) is spread over many blocks
/// instead of keeping one busy while the rest of the GPU waits
constexpr int kChunk = 4096;

// One block per chunk of one parameter's terms, its threads taking every kBlock-th term, then a fixed
// order reduction into that chunk's partial sum
__global__ void ParamGradChunks(const int* __restrict__ ChunkStart, const int* __restrict__ ChunkEnd,
                                const int* __restrict__ TermSig, const int* __restrict__ TermDeriv,
                                const int* __restrict__ TermPool, const double* __restrict__ Deriv, const double* __restrict__ Pool,
                                const double* __restrict__ Share, double* __restrict__ ChunkSum, const int NChunks) {
  __shared__ double Partial[kBlock];
  const int c = blockIdx.x;
  const int t = threadIdx.x;
  double g = 0.;
  if (c < NChunks) {
    for (int k = ChunkStart[c] + t; k < ChunkEnd[c]; k += kBlock) {
      const double d = Deriv[TermDeriv[k]];
      if (d != 0.) g += Share[TermSig[k]] * d / Pool[TermPool[k]];
    }
  }
  Partial[t] = g;
  __syncthreads();
  for (int Half = kBlock / 2; Half > 0; Half /= 2) {
    if (t < Half) Partial[t] += Partial[t + Half];
    __syncthreads();
  }
  if (t == 0 && c < NChunks) ChunkSum[c] = Partial[0];
}

// Each parameter's chunks, added in order, so the result is the same from run to run
__global__ void ParamGradSum(const int* __restrict__ ParChunkStart, const double* __restrict__ ChunkSum,
                             double* __restrict__ Grad, const int NPar) {
  const int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= NPar) return;
  double g = 0.;
  for (int c = ParChunkStart[p]; c < ParChunkStart[p + 1]; ++c) g += ChunkSum[c];
  Grad[p] = g;
}


// ------------------------------------------------------------------------------------------------
// Full pass (the shifts moved): each live event's shifted variables from the linear shift table,
// its cuts and its bin, then the weights summed per bin in a fixed order
// ------------------------------------------------------------------------------------------------

/// Variable k of live event i. __dmul_rn and __dadd_rn are never fused into an FMA, so this is the
/// host's LinearShiftTable::Value bit for bit (the host build has no FMA), and every event lands in
/// the same bin as on the CPU
__device__ inline double SweepVariable(const int i, const int k, const unsigned int* __restrict__ VarStart,
                                       const double* __restrict__ Base, const unsigned int* __restrict__ TermStart,
                                       const unsigned short* __restrict__ TermShift, const double* __restrict__ TermCoef,
                                       const double* __restrict__ Delta) {
  const unsigned int Var = VarStart[i] + static_cast<unsigned int>(k);
  double v = Base[Var];
  for (unsigned int t = TermStart[Var]; t < TermStart[Var + 1]; ++t) v = __dadd_rn(v, __dmul_rn(Delta[TermShift[t]], TermCoef[t]));
  return v;
}

/// One thread per live event: its bin (-1 if cut or outside the binning), its weight, and the key it
/// is sorted by (its bin, or NBins for one that adds nothing)
__global__ void SweepBins(const int NLive, const int NBins, const int* __restrict__ Sample, const unsigned int* __restrict__ VarStart,
                          const double* __restrict__ Base, const unsigned int* __restrict__ TermStart,
                          const unsigned short* __restrict__ TermShift, const double* __restrict__ TermCoef,
                          const double* __restrict__ Delta, const int* __restrict__ CutStart, const double* __restrict__ CutLo,
                          const double* __restrict__ CutHi, const int* __restrict__ EdgeStart, const double* __restrict__ Edges,
                          const int* __restrict__ Stride, const int* __restrict__ Offset, const int* __restrict__ Sig,
                          const double* __restrict__ Const, const double* __restrict__ SigProd, int* __restrict__ EventBin,
                          unsigned int* __restrict__ Key, int* __restrict__ Index, double* __restrict__ Weight) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= NLive) return;
  const int s = Sample[i];
  int Bin = -1;
  bool Pass = true;
  for (int c = CutStart[s]; c < CutStart[s + 1]; ++c) {
    const double v = SweepVariable(i, 1 + c - CutStart[s], VarStart, Base, TermStart, TermShift, TermCoef, Delta);
    if (v < CutLo[c] || v >= CutHi[c]) { Pass = false; break; }
  }
  if (Pass) {
    const double v = SweepVariable(i, 0, VarStart, Base, TermStart, TermShift, TermCoef, Delta);
    const int Lo = EdgeStart[s], Hi = EdgeStart[s + 1] - 1;
    // SampleBinningInfo::FindBin: outside [first edge, last edge) is out, else the upper_bound search
    if (v >= Edges[Lo] && v < Edges[Hi]) {
      int a = Lo, b = Hi;   // first edge > v lies in (a, b]
      while (b - a > 1) {
        const int m = (a + b) / 2;
        if (Edges[m] > v) b = m; else a = m;
      }
      Bin = (a - Lo) * Stride[s] + Offset[s];
    }
  }
  EventBin[i] = Bin;
  const double W = Const[i] * SigProd[Sig[i]];
  Weight[i] = W;
  Key[i] = (Bin < 0 || !(W > 0.)) ? static_cast<unsigned int>(NBins) : static_cast<unsigned int>(Bin);
  Index[i] = i;
}

/// The weights and squared weights in sorted order
__global__ void GatherSorted(const int N, const int* __restrict__ SortedIndex, const double* __restrict__ Weight,
                             double* __restrict__ SortedW, double* __restrict__ SortedW2) {
  const int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= N) return;
  const double W = Weight[SortedIndex[j]];
  SortedW[j] = W;
  SortedW2[j] = W * W;
}

/// Where each bin's run of sorted keys starts (NBins + 1 entries; empty bins get an empty run)
__global__ void FindBinStarts(const int N, const int NBins, const unsigned int* __restrict__ SortedKey, int* __restrict__ BinStart) {
  const int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j > N) return;
  const int Prev = (j == 0) ? -1 : static_cast<int>(SortedKey[j - 1]);
  const int Cur = (j == N) ? NBins : static_cast<int>(SortedKey[j]);
  for (int b = Prev + 1; b <= Cur && b <= NBins; ++b) BinStart[b] = j;
}

} // namespace

struct CudaReweightDevice::Impl {
  int NBins = 0, NSig = 0, NPar = 0, NChunks = 0;
  std::size_t NPool = 0, NDeriv = 0;
  cudaStream_t Stream = nullptr;
  DeviceArray<int> SigStart, SigPool, BinStart, GroupSig, SigGroupStart, SigGroup, GroupBin, TermSig, TermDeriv, TermPool;
  /// Each chunk's first and one past last term, and where each parameter's chunks start
  DeviceArray<int> ChunkStart, ChunkEnd, ParChunkStart;
  DeviceArray<double> Pool, SigProd, GroupConst, GroupConst2, MC, W2, BinFactor, Share, Deriv, Grad, ChunkSum;
  PinnedArray<double> HostPool, HostDeriv, HostOut;

  /// @name Full pass, see SetSweep
  /// @{
  bool HaveSweep = false;
  int NLive = 0;
  int KeyBits = 1;
  DeviceArray<int> SwSample, SwCutStart, SwEdgeStart, SwStride, SwOffset, SwSig, EventBin, Index, SortedIndex, SwBinStart;
  DeviceArray<unsigned int> SwVarStart, SwTermStart, Key, SortedKey;
  DeviceArray<unsigned short> SwTermShift;
  DeviceArray<double> SwBase, SwTermCoef, SwCutLo, SwCutHi, SwEdges, SwConst, Delta, Weight, SortedW, SortedW2;
  DeviceArray<unsigned char> Temp;
  std::size_t TempBytes = 0;
  PinnedArray<double> HostDelta;
  PinnedArray<int> HostBins;
  /// @}
};

CudaReweightDevice::CudaReweightDevice(const std::vector<int>& SigStart, const std::vector<int>& SigPool, const std::size_t NPool,
                                       const int NBins)
    : D(new Impl()) {
  Check(cudaStreamCreate(&D->Stream), "stream");
  D->NBins = NBins;
  D->NSig = static_cast<int>(SigStart.size()) - 1;
  D->NPool = NPool;
  D->SigStart.Upload(SigStart);
  D->SigPool.Upload(SigPool);
  D->Pool.Allocate(NPool);
  D->SigProd.Allocate(static_cast<std::size_t>(D->NSig));
  D->MC.Allocate(static_cast<std::size_t>(NBins));
  D->W2.Allocate(static_cast<std::size_t>(NBins));
  D->BinFactor.Allocate(static_cast<std::size_t>(NBins));
  D->Share.Allocate(static_cast<std::size_t>(D->NSig));
  D->HostPool.Allocate(NPool);
  D->HostOut.Allocate(2 * static_cast<std::size_t>(NBins));
}

CudaReweightDevice::~CudaReweightDevice() {
  if (D && D->Stream) cudaStreamDestroy(D->Stream);
}

void CudaReweightDevice::SetGroups(const std::vector<int>& BinStart, const std::vector<int>& Sig, const std::vector<double>& Const,
                                   const std::vector<double>& Const2, const std::vector<int>& SigGroupStart,
                                   const std::vector<int>& SigGroup, const std::vector<int>& GroupBin) {
  D->BinStart.Upload(BinStart);
  D->GroupSig.Upload(Sig);
  D->GroupConst.Upload(Const);
  D->GroupConst2.Upload(Const2);
  D->SigGroupStart.Upload(SigGroupStart);
  D->SigGroup.Upload(SigGroup);
  D->GroupBin.Upload(GroupBin);
}

void CudaReweightDevice::SetTerms(const std::vector<int>& ParStart, const std::vector<int>& TermSig, const std::vector<int>& TermDeriv,
                                  const std::vector<int>& TermPool, const std::size_t NDeriv) {
  D->NPar = static_cast<int>(ParStart.size()) - 1;
  D->NDeriv = NDeriv;
  // Cut each parameter's terms into chunks of at most kChunk, keeping a parameter's chunks together
  std::vector<int> ChunkStart, ChunkEnd, ParChunkStart(1, 0);
  for (int p = 0; p < D->NPar; ++p) {
    const int First = ParStart[static_cast<std::size_t>(p)], Last = ParStart[static_cast<std::size_t>(p) + 1];
    for (int k = First; k < Last; k += kChunk) {
      ChunkStart.push_back(k);
      ChunkEnd.push_back(std::min(k + kChunk, Last));
    }
    ParChunkStart.push_back(static_cast<int>(ChunkStart.size()));
  }
  D->NChunks = static_cast<int>(ChunkStart.size());
  D->ChunkStart.Upload(ChunkStart);
  D->ChunkEnd.Upload(ChunkEnd);
  D->ParChunkStart.Upload(ParChunkStart);
  D->ChunkSum.Allocate(static_cast<std::size_t>(D->NChunks));
  D->TermSig.Upload(TermSig);
  D->TermDeriv.Upload(TermDeriv);
  D->TermPool.Upload(TermPool);
  D->Deriv.Allocate(NDeriv);
  D->Grad.Allocate(static_cast<std::size_t>(D->NPar));
  D->HostDeriv.Allocate(NDeriv);
}

double* CudaReweightDevice::PoolBuffer() { return D->HostPool.Ptr; }
double* CudaReweightDevice::DerivBuffer() { return D->HostDeriv.Ptr; }

void CudaReweightDevice::Fill(double* MC, double* W2) {
  const std::size_t NB = static_cast<std::size_t>(D->NBins);
  Check(cudaMemcpyAsync(D->Pool.Ptr, D->HostPool.Ptr, D->NPool * sizeof(double), cudaMemcpyHostToDevice, D->Stream), "pool");
  SigProducts<<<Blocks(static_cast<std::size_t>(D->NSig)), kBlock, 0, D->Stream>>>(D->Pool.Ptr, D->SigStart.Ptr, D->SigPool.Ptr, D->SigProd.Ptr, D->NSig);
  BinSums<<<static_cast<unsigned int>(std::max<std::size_t>(1, NB)), kBlock, 0, D->Stream>>>(D->BinStart.Ptr, D->GroupSig.Ptr, D->GroupConst.Ptr, D->GroupConst2.Ptr, D->SigProd.Ptr,
                                                D->MC.Ptr, D->W2.Ptr, D->NBins);
  // One copy back for both arrays, through pinned memory
  Check(cudaMemcpyAsync(D->HostOut.Ptr, D->MC.Ptr, NB * sizeof(double), cudaMemcpyDeviceToHost, D->Stream), "mc");
  Check(cudaMemcpyAsync(D->HostOut.Ptr + NB, D->W2.Ptr, NB * sizeof(double), cudaMemcpyDeviceToHost, D->Stream), "w2");
  Check(cudaStreamSynchronize(D->Stream), "fill");
  for (std::size_t b = 0; b < NB; ++b) { MC[b] = D->HostOut.Ptr[b]; W2[b] = D->HostOut.Ptr[NB + b]; }
}

void CudaReweightDevice::Gradient(const std::vector<double>& BinFactor, std::vector<double>& Grad) {
  Check(cudaMemcpyAsync(D->BinFactor.Ptr, BinFactor.data(), BinFactor.size() * sizeof(double), cudaMemcpyHostToDevice, D->Stream), "factor");
  Check(cudaMemcpyAsync(D->Deriv.Ptr, D->HostDeriv.Ptr, D->NDeriv * sizeof(double), cudaMemcpyHostToDevice, D->Stream), "deriv");
  SigShare<<<Blocks(static_cast<std::size_t>(D->NSig)), kBlock, 0, D->Stream>>>(D->SigGroupStart.Ptr, D->SigGroup.Ptr, D->GroupBin.Ptr,
                                                                                D->GroupConst.Ptr, D->SigProd.Ptr, D->BinFactor.Ptr,
                                                                                D->Share.Ptr, D->NSig);
  ParamGradChunks<<<static_cast<unsigned int>(std::max(1, D->NChunks)), kBlock, 0, D->Stream>>>(
      D->ChunkStart.Ptr, D->ChunkEnd.Ptr, D->TermSig.Ptr, D->TermDeriv.Ptr, D->TermPool.Ptr, D->Deriv.Ptr, D->Pool.Ptr,
      D->Share.Ptr, D->ChunkSum.Ptr, D->NChunks);
  ParamGradSum<<<Blocks(static_cast<std::size_t>(std::max(1, D->NPar))), kBlock, 0, D->Stream>>>(D->ParChunkStart.Ptr, D->ChunkSum.Ptr,
                                                                                                  D->Grad.Ptr, D->NPar);
  Grad.resize(static_cast<std::size_t>(D->NPar));
  Check(cudaMemcpyAsync(Grad.data(), D->Grad.Ptr, Grad.size() * sizeof(double), cudaMemcpyDeviceToHost, D->Stream), "grad");
  Check(cudaStreamSynchronize(D->Stream), "gradient");
}

void CudaReweightDevice::SetSweep(const SweepTables& T, const std::vector<int>& Sig, const std::vector<double>& Const) {
  if (static_cast<int>(Sig.size()) != T.NLive || static_cast<int>(Const.size()) != T.NLive) {
    throw std::runtime_error("CudaReweightDevice::SetSweep: tables and model disagree on the number of live events");
  }
  D->NLive = T.NLive;
  D->SwSample.Upload(T.Sample);
  D->SwVarStart.Upload(T.VarStart);
  D->SwBase.Upload(T.Base);
  D->SwTermStart.Upload(T.TermStart);
  D->SwTermShift.Upload(T.TermShift);
  D->SwTermCoef.Upload(T.TermCoef);
  D->SwCutStart.Upload(T.CutStart);
  D->SwCutLo.Upload(T.CutLo);
  D->SwCutHi.Upload(T.CutHi);
  D->SwEdgeStart.Upload(T.EdgeStart);
  D->SwEdges.Upload(T.Edges);
  D->SwStride.Upload(T.Stride);
  D->SwOffset.Upload(T.Offset);
  D->SwSig.Upload(Sig);
  D->SwConst.Upload(Const);
  const std::size_t N = static_cast<std::size_t>(T.NLive);
  D->Delta.Allocate(static_cast<std::size_t>(std::max(1, T.NShifts)));
  D->HostDelta.Allocate(static_cast<std::size_t>(std::max(1, T.NShifts)));
  D->EventBin.Allocate(N);
  D->Index.Allocate(N);
  D->SortedIndex.Allocate(N);
  D->Key.Allocate(N);
  D->SortedKey.Allocate(N);
  D->Weight.Allocate(N);
  D->SortedW.Allocate(N);
  D->SortedW2.Allocate(N);
  D->SwBinStart.Allocate(static_cast<std::size_t>(D->NBins) + 1);
  D->HostBins.Allocate(N);
  // Keys run up to NBins (the bin of an event that adds nothing)
  D->KeyBits = 1;
  while ((1u << D->KeyBits) <= static_cast<unsigned int>(D->NBins)) ++D->KeyBits;
  std::size_t SortBytes = 0, SumBytes = 0;
  Check(cub::DeviceRadixSort::SortPairs(nullptr, SortBytes, D->Key.Ptr, D->SortedKey.Ptr, D->Index.Ptr, D->SortedIndex.Ptr,
                                        T.NLive, 0, D->KeyBits), "sort size");
  Check(cub::DeviceSegmentedReduce::Sum(nullptr, SumBytes, D->SortedW.Ptr, D->MC.Ptr, D->NBins, D->SwBinStart.Ptr,
                                        D->SwBinStart.Ptr + 1), "reduce size");
  D->TempBytes = std::max<std::size_t>({SortBytes, SumBytes, 1});
  D->Temp.Allocate(D->TempBytes);
  D->HaveSweep = true;
}

bool CudaReweightDevice::HasSweep() const { return D->HaveSweep; }

void CudaReweightDevice::FullPass(const std::vector<double>& Delta, double* MC, double* W2) {
  if (!D->HaveSweep) throw std::runtime_error("CudaReweightDevice::FullPass before SetSweep");
  const std::size_t NB = static_cast<std::size_t>(D->NBins);
  const int N = D->NLive;
  std::copy(Delta.begin(), Delta.end(), D->HostDelta.Ptr);
  Check(cudaMemcpyAsync(D->Delta.Ptr, D->HostDelta.Ptr, Delta.size() * sizeof(double), cudaMemcpyHostToDevice, D->Stream), "delta");
  Check(cudaMemcpyAsync(D->Pool.Ptr, D->HostPool.Ptr, D->NPool * sizeof(double), cudaMemcpyHostToDevice, D->Stream), "pool");
  SigProducts<<<Blocks(static_cast<std::size_t>(D->NSig)), kBlock, 0, D->Stream>>>(D->Pool.Ptr, D->SigStart.Ptr, D->SigPool.Ptr, D->SigProd.Ptr, D->NSig);
  SweepBins<<<Blocks(static_cast<std::size_t>(N)), kBlock, 0, D->Stream>>>(
      N, D->NBins, D->SwSample.Ptr, D->SwVarStart.Ptr, D->SwBase.Ptr, D->SwTermStart.Ptr, D->SwTermShift.Ptr, D->SwTermCoef.Ptr,
      D->Delta.Ptr, D->SwCutStart.Ptr, D->SwCutLo.Ptr, D->SwCutHi.Ptr, D->SwEdgeStart.Ptr, D->SwEdges.Ptr, D->SwStride.Ptr,
      D->SwOffset.Ptr, D->SwSig.Ptr, D->SwConst.Ptr, D->SigProd.Ptr, D->EventBin.Ptr, D->Key.Ptr, D->Index.Ptr, D->Weight.Ptr);
  // A stable sort by bin keeps each bin's events in live order, so the sums below always add the
  // same numbers in the same order
  std::size_t Bytes = D->TempBytes;
  Check(cub::DeviceRadixSort::SortPairs(D->Temp.Ptr, Bytes, D->Key.Ptr, D->SortedKey.Ptr, D->Index.Ptr, D->SortedIndex.Ptr, N, 0,
                                        D->KeyBits, D->Stream), "sort");
  GatherSorted<<<Blocks(static_cast<std::size_t>(N)), kBlock, 0, D->Stream>>>(N, D->SortedIndex.Ptr, D->Weight.Ptr, D->SortedW.Ptr,
                                                                              D->SortedW2.Ptr);
  FindBinStarts<<<Blocks(static_cast<std::size_t>(N) + 1), kBlock, 0, D->Stream>>>(N, D->NBins, D->SortedKey.Ptr, D->SwBinStart.Ptr);
  Bytes = D->TempBytes;
  Check(cub::DeviceSegmentedReduce::Sum(D->Temp.Ptr, Bytes, D->SortedW.Ptr, D->MC.Ptr, D->NBins, D->SwBinStart.Ptr, D->SwBinStart.Ptr + 1,
                                        D->Stream), "sum");
  Bytes = D->TempBytes;
  Check(cub::DeviceSegmentedReduce::Sum(D->Temp.Ptr, Bytes, D->SortedW2.Ptr, D->W2.Ptr, D->NBins, D->SwBinStart.Ptr, D->SwBinStart.Ptr + 1,
                                        D->Stream), "sum w2");
  Check(cudaMemcpyAsync(D->HostOut.Ptr, D->MC.Ptr, NB * sizeof(double), cudaMemcpyDeviceToHost, D->Stream), "mc");
  Check(cudaMemcpyAsync(D->HostOut.Ptr + NB, D->W2.Ptr, NB * sizeof(double), cudaMemcpyDeviceToHost, D->Stream), "w2");
  Check(cudaStreamSynchronize(D->Stream), "full pass");
  for (std::size_t b = 0; b < NB; ++b) { MC[b] = D->HostOut.Ptr[b]; W2[b] = D->HostOut.Ptr[NB + b]; }
}

void CudaReweightDevice::GetEventBins(int* Out) {
  const std::size_t N = static_cast<std::size_t>(D->NLive);
  Check(cudaMemcpyAsync(D->HostBins.Ptr, D->EventBin.Ptr, N * sizeof(int), cudaMemcpyDeviceToHost, D->Stream), "bins");
  Check(cudaStreamSynchronize(D->Stream), "bins");
  std::copy(D->HostBins.Ptr, D->HostBins.Ptr + N, Out);
}

std::string CudaReweightDevice::DeviceName() {
  int Device = 0;
  cudaDeviceProp Props;
  if (cudaGetDevice(&Device) != cudaSuccess || cudaGetDeviceProperties(&Props, Device) != cudaSuccess) return "unknown";
  return Props.name;
}

bool CudaReweightDevice::Available() {
  int N = 0;
  return cudaGetDeviceCount(&N) == cudaSuccess && N > 0;
}

std::size_t CudaReweightDevice::Bytes() const {
  return D->SigStart.Bytes() + D->SigPool.Bytes() + D->BinStart.Bytes() + D->GroupSig.Bytes() + D->SigGroupStart.Bytes() +
         D->SigGroup.Bytes() + D->GroupBin.Bytes() + D->TermSig.Bytes() + D->TermDeriv.Bytes() +
         D->TermPool.Bytes() + D->Pool.Bytes() + D->SigProd.Bytes() + D->GroupConst.Bytes() + D->GroupConst2.Bytes() +
         D->MC.Bytes() + D->W2.Bytes() + D->BinFactor.Bytes() + D->Share.Bytes() + D->Deriv.Bytes() + D->Grad.Bytes() +
         D->ChunkStart.Bytes() + D->ChunkEnd.Bytes() + D->ParChunkStart.Bytes() + D->ChunkSum.Bytes() +
         D->SwSample.Bytes() + D->SwVarStart.Bytes() + D->SwBase.Bytes() + D->SwTermStart.Bytes() + D->SwTermShift.Bytes() +
         D->SwTermCoef.Bytes() + D->SwCutStart.Bytes() + D->SwCutLo.Bytes() + D->SwCutHi.Bytes() + D->SwEdgeStart.Bytes() +
         D->SwEdges.Bytes() + D->SwStride.Bytes() + D->SwOffset.Bytes() + D->SwSig.Bytes() + D->SwConst.Bytes() +
         D->Delta.Bytes() + D->EventBin.Bytes() + D->Index.Bytes() + D->SortedIndex.Bytes() + D->Key.Bytes() +
         D->SortedKey.Bytes() + D->Weight.Bytes() + D->SortedW.Bytes() + D->SortedW2.Bytes() + D->SwBinStart.Bytes() +
         D->Temp.Bytes();
}
