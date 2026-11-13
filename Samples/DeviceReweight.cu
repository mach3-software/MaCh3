// GPU version of SampleHandlerBase's cached, grouped reweight, see DeviceReweight.h
#include "Samples/DeviceReweight.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

#include <cub/block/block_reduce.cuh>
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_reduce.cuh>
#include <cub/device/device_segmented_reduce.cuh>
#include <cuda_runtime.h>

namespace {

void Check(const cudaError_t Err, const char* What, const int Line) {
  if (Err != cudaSuccess) {
    throw std::runtime_error(std::string("DeviceReweight: ") + What + " failed at line " + std::to_string(Line) +
                             ": " + cudaGetErrorString(Err));
  }
}
#define DR_CHECK(x) Check((x), #x, __LINE__)

constexpr int kThreads = 256;

/// One thread per spline: the cubic of BinnedSplineHandler::CalcSplineWeights, and its derivative
__global__ void SplineKernel(const int NSplines, const int MaxKnots, const short int* __restrict__ ParamNo,
                             const unsigned int* __restrict__ NKnots, const double* __restrict__ Many,
                             const double* __restrict__ X, const short int* __restrict__ Segments,
                             const int* __restrict__ GlobalPar, double* __restrict__ Pool, const int WOffset,
                             const int GOffset) {
  const int s = blockIdx.x * blockDim.x + threadIdx.x;
  if (s >= NSplines) return;
  const int Param = ParamNo[s];
  const int Segment = Segments[Param];
  const double* C = Many + static_cast<size_t>(NKnots[s] + Segment) * 4;
  const double y = C[0], b = C[1], c = C[2], d = C[3];
  const double dx = Pool[GlobalPar[Param]] - X[Param * MaxKnots + Segment];
  Pool[WOffset + s] = fma(dx, fma(dx, fma(dx, d, c), b), y);
  Pool[GOffset + s] = fma(dx, fma(dx, 3. * d, 2. * c), b);
}

/// Gradient of a full pass (only asked for at the start of a fit, the sweeps pause gradients): grid
/// stride over events, d(-lnL)/d(par) += (1 - data/MC) * event weight * d(w)/d(par) / w, summed per
/// block in shared memory first, so the global array only takes NPars atomics per block
__global__ void GradKernel(const int NGroups, const int NPars, const int* __restrict__ GroupBin,
                           const int* __restrict__ GroupSig, const double* __restrict__ GroupW,
                           const int* __restrict__ SigDStart, const int* __restrict__ SigDPar,
                           const int* __restrict__ SigDDeriv, const int* __restrict__ SigDWeight,
                           const double* __restrict__ Pool, const double* __restrict__ MC,
                           const double* __restrict__ Data, const double LowMC, double* __restrict__ Grad) {
  extern __shared__ double Local[];
  for (int p = threadIdx.x; p < NPars; p += blockDim.x) Local[p] = 0.;
  __syncthreads();
  for (int g = blockIdx.x * blockDim.x + threadIdx.x; g < NGroups; g += gridDim.x * blockDim.x) {
    const double W = GroupW[g];
    if (W == 0.) continue;
    const int Bin = GroupBin[g];
    const double Mu = MC[Bin];
    if (Mu < LowMC) continue;
    const double F = (1. - Data[Bin] / Mu) * W;
    const int Sig = GroupSig[g];
    for (int k = SigDStart[Sig]; k < SigDStart[Sig + 1]; ++k) {
      atomicAdd(&Local[SigDPar[k]], F * Pool[SigDDeriv[k]] / Pool[SigDWeight[k]]);
    }
  }
  __syncthreads();
  for (int p = threadIdx.x; p < NPars; p += blockDim.x) {
    if (Local[p] != 0.) atomicAdd(&Grad[p], Local[p]);
  }
}

/// One thread per live event, grid stride: the CPU full pass. Shifted variables are built in the
/// CPU's order with unfused multiply and add, the cuts and the bin search follow
/// IsEventSelectedFast and FindGlobalBinFast, and the bins are summed in shared memory per block
__global__ void FullPassKernel(const int NLive, const int NSlots, const int NBins, const int* __restrict__ Sample,
                               const int* __restrict__ Sig, const double* __restrict__ ConstW,
                               const double* __restrict__ Nominal, const int* __restrict__ CoefStart,
                               const int* __restrict__ CoefShift, const double* __restrict__ Coef,
                               const double* __restrict__ ShiftVal, const int* __restrict__ CutStart,
                               const int* __restrict__ CutSlot, const double* __restrict__ CutLo,
                               const double* __restrict__ CutHi, const int* __restrict__ EdgeStart,
                               const double* __restrict__ Edges, const int* __restrict__ Stride,
                               const int* __restrict__ Offset, const int* __restrict__ SigWStart,
                               const int* __restrict__ SigWIdx, const double* __restrict__ Pool,
                               const int UOF, const int NotSel, const bool FillW2, int* __restrict__ EventBin,
                               double* __restrict__ EventW, double* __restrict__ MC, double* __restrict__ W2) {
  extern __shared__ double Hist[];
  const int NH = FillW2 ? 2 * NBins : NBins;
  for (int b = threadIdx.x; b < NH; b += blockDim.x) Hist[b] = 0.;
  __syncthreads();
  constexpr int kMaxSlots = 8;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < NLive; i += gridDim.x * blockDim.x) {
    double V[kMaxSlots];
    for (int j = 0; j < NSlots; ++j) {
      const int Cell = i * NSlots + j;
      double x = Nominal[Cell];
      for (int k = CoefStart[Cell]; k < CoefStart[Cell + 1]; ++k) x = __dadd_rn(x, __dmul_rn(ShiftVal[CoefShift[k]], Coef[k]));
      V[j] = x;
    }
    const int S = Sample[i];
    bool Selected = true;
    for (int c = CutStart[S]; c < CutStart[S + 1]; ++c) {
      const double x = V[CutSlot[c]];
      if (x < CutLo[c] || x >= CutHi[c]) { Selected = false; break; }
    }
    int Bin = NotSel;
    double W = 0.;
    if (Selected) {
      const double* E = Edges + EdgeStart[S];
      const int NB = EdgeStart[S + 1] - EdgeStart[S] - 1;
      const double x = V[0];
      if (x < E[0] || x >= E[NB]) Bin = UOF;
      else {
        // upper_bound over the edges, minus one
        int Lo = 0, Hi = NB + 1;
        while (Lo < Hi) { const int Mid = (Lo + Hi) / 2; if (E[Mid] <= x) Lo = Mid + 1; else Hi = Mid; }
        Bin = (Lo - 1) * Stride[S] + Offset[S];
      }
      if (Bin > UOF) {
        const int G = Sig[i];
        double Shared = 1.;
        for (int k = SigWStart[G]; k < SigWStart[G + 1]; ++k) Shared *= Pool[SigWIdx[k]];
        W = ConstW[i] * Shared;
        if (W > 0.) {
          atomicAdd(&Hist[Bin], W);
          if (FillW2) atomicAdd(&Hist[NBins + Bin], W * W);
        } else W = 0.;
      }
    }
    EventBin[i] = Bin;
    EventW[i] = W;
  }
  __syncthreads();
  for (int b = threadIdx.x; b < NBins; b += blockDim.x) {
    if (Hist[b] != 0.) atomicAdd(&MC[b], Hist[b]);
    if (FillW2 && Hist[NBins + b] != 0.) atomicAdd(&W2[b], Hist[NBins + b]);
  }
}

/// Gradient, second half: the derivative entries, grouped by parameter, cut into chunks of at most
/// kChunk so a parameter every event depends on (a fiducial volume norm) doesn't leave one block
/// working alone. One block per chunk: the sum of SigSum * d(w)/d(par) / w over the chunk
constexpr int kChunk = 256;
__global__ void ParChunkKernel(const int* __restrict__ ChunkBegin, const int* __restrict__ ParSig,
                               const int* __restrict__ ParDeriv, const int* __restrict__ ParWeight,
                               const double* __restrict__ SigSum, const double* __restrict__ Pool,
                               double* __restrict__ Partial) {
  using Reduce = cub::BlockReduce<double, kThreads>;
  __shared__ typename Reduce::TempStorage Temp;
  const int c = blockIdx.x;
  double Acc = 0.;
  for (int k = ChunkBegin[c] + threadIdx.x; k < ChunkBegin[c + 1]; k += blockDim.x) {
    const double A = SigSum[ParSig[k]];
    if (A != 0.) Acc += A * Pool[ParDeriv[k]] / Pool[ParWeight[k]];
  }
  const double Total = Reduce(Temp).Sum(Acc);
  if (threadIdx.x == 0) Partial[c] = Total;
}

/// Gradient, last: one thread per parameter adds its chunks, in order. No atomics anywhere
__global__ void ParSumKernel(const int NPars, const int* __restrict__ ParChunkStart, const double* __restrict__ Partial,
                             double* __restrict__ Grad) {
  const int Par = blockIdx.x * blockDim.x + threadIdx.x;
  if (Par >= NPars) return;
  double Acc = 0.;
  for (int c = ParChunkStart[Par]; c < ParChunkStart[Par + 1]; ++c) Acc += Partial[c];
  Grad[Par] = Acc;
}

/// Groups from the last full pass, step 1: a sort key for every live event, (bin, signature, sign of its
/// constant weight), so a sort puts each group's events together in the CPU's order (bin, then
/// signature, positive before negative). Events that fill nothing get a key past every real one
__global__ void GroupKeyKernel(const int NLive, const int NSig, const int UOF, const int* __restrict__ EventBin,
                               const int* __restrict__ Sig, const double* __restrict__ ConstW,
                               const unsigned long long Invalid, unsigned long long* __restrict__ Keys,
                               int* __restrict__ Index, int* __restrict__ NValid) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= NLive) return;
  const int Bin = EventBin[i];
  const double C = ConstW[i];
  const bool Valid = (Bin > UOF) && (C != 0.);
  Keys[i] = Valid ? ((static_cast<unsigned long long>(Bin) * static_cast<unsigned long long>(NSig) + static_cast<unsigned long long>(Sig[i])) << 1) | (C < 0. ? 1ull : 0ull)
                  : Invalid;
  Index[i] = i;
  if (Valid) atomicAdd(NValid, 1);
}

/// Step 2: each sorted event's constant weight and its square
__global__ void GroupGatherKernel(const int N, const int* __restrict__ SortedIndex, const double* __restrict__ ConstW,
                                  double* __restrict__ C1, double* __restrict__ C2) {
  const int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= N) return;
  const double C = ConstW[SortedIndex[j]];
  C1[j] = C;
  C2[j] = C * C;
}

/// Step 4: each group's signature and bin back from its key
__global__ void GroupDecodeKernel(const int NGroups, const int NSig, const unsigned long long* __restrict__ Keys,
                                  int* __restrict__ GroupSig, int* __restrict__ GroupBin) {
  const int g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= NGroups) return;
  const unsigned long long k = Keys[g] >> 1;
  GroupSig[g] = static_cast<int>(k % static_cast<unsigned long long>(NSig));
  GroupBin[g] = static_cast<int>(k / static_cast<unsigned long long>(NSig));
}

/// Step 5: where each bin's groups start, by binary search over the (sorted) group bins
__global__ void GroupBinStartKernel(const int NBins, const int NGroups, const int* __restrict__ GroupBin, int* __restrict__ BinStart) {
  const int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b > NBins) return;
  int Lo = 0, Hi = NGroups;
  while (Lo < Hi) { const int Mid = (Lo + Hi) / 2; if (GroupBin[Mid] < b) Lo = Mid + 1; else Hi = Mid; }
  BinStart[b] = Lo;
}

/// Leapfrog pass, step 1: one thread per group, its weight and its share of W2. Bins hold very different
/// numbers of groups, so spreading groups rather than bins over the GPU keeps every thread busy
__global__ void GroupWeightKernel(const int NGroups, const int* __restrict__ GroupSig, const double* __restrict__ GroupConst,
                                  const double* __restrict__ GroupConst2, const int* __restrict__ SigWStart,
                                  const int* __restrict__ SigWIdx, const double* __restrict__ Pool,
                                  double* __restrict__ GroupW, double* __restrict__ GroupW2, const bool FillW2) {
  const int g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= NGroups) return;
  const int Sig = GroupSig[g];
  double Shared = 1.;
  for (int k = SigWStart[Sig]; k < SigWStart[Sig + 1]; ++k) Shared *= Pool[SigWIdx[k]];
  double W = GroupConst[g] * Shared;
  // As the CPU loop: a group whose weight isn't positive contributes nothing, gradient included
  double W2v = 0.;
  if (!(W > 0.)) W = 0.;
  else if (FillW2) W2v = GroupConst2[g] * Shared * Shared;
  GroupW[g] = W;
  if (FillW2) GroupW2[g] = W2v;
}

/// Step 3 (step 2 sums each bin with a segmented reduction): one thread per group, its
/// (1 - data/MC) * weight onto its signature, the weighting GetGradLikelihood gives the bin
__global__ void GroupSigSumKernel(const int NGroups, const int* __restrict__ GroupBin, const int* __restrict__ GroupSig,
                                  const double* __restrict__ GroupW, const double* __restrict__ MC,
                                  const double* __restrict__ Data, const double LowMC, double* __restrict__ SigSum) {
  const int g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= NGroups) return;
  const double W = GroupW[g];
  if (W == 0.) return;
  const int Bin = GroupBin[g];
  const double Mu = MC[Bin];
  if (Mu < LowMC) return;
  atomicAdd(&SigSum[GroupSig[g]], (1. - Data[Bin] / Mu) * W);
}

template <typename T>
T* Upload(const std::vector<T>& Host, std::size_t& Bytes) {
  T* Dev = nullptr;
  const std::size_t N = Host.empty() ? 1 : Host.size();
  DR_CHECK(cudaMalloc(&Dev, N * sizeof(T)));
  if (!Host.empty()) DR_CHECK(cudaMemcpy(Dev, Host.data(), Host.size() * sizeof(T), cudaMemcpyHostToDevice));
  Bytes += N * sizeof(T);
  return Dev;
}

}  // namespace

class DeviceReweightImpl {
 public:
  explicit DeviceReweightImpl(const DeviceReweightStatic& S) : St(S) {
    NPool = S.NPars + 2 * S.NSplines + S.NOther + 2;
    OtherOffset = S.NPars + 2 * S.NSplines;
    DR_CHECK(cudaMalloc(&Pool, static_cast<std::size_t>(NPool) * sizeof(double)));
    Bytes += static_cast<std::size_t>(NPool) * sizeof(double);
    const double Constants[2] = {1., 0.};
    DR_CHECK(cudaMemcpy(Pool + NPool - 2, Constants, 2 * sizeof(double), cudaMemcpyHostToDevice));
    if (S.NSplines > 0) {
      XCoeff = Upload(S.XCoeff, Bytes);
      Many = Upload(S.ManyCoeff, Bytes);
      NKnots = Upload(S.NKnots, Bytes);
      ParamNo = Upload(S.ParamNo, Bytes);
      GlobalPar = Upload(S.SplineGlobalPar, Bytes);
      DR_CHECK(cudaMalloc(&Segments, static_cast<std::size_t>(S.NSplineParams) * sizeof(short int)));
      Bytes += static_cast<std::size_t>(S.NSplineParams) * sizeof(short int);
    }
    // The derivative entries regrouped by parameter, for ParGradKernel
    NSig = static_cast<int>(S.SigWStart.size()) - 1;
    {
      std::vector<int> PStart(static_cast<std::size_t>(S.NPars) + 1, 0);
      for (const int Par : S.SigDPar) ++PStart[static_cast<std::size_t>(Par) + 1];
      for (int p = 0; p < S.NPars; ++p) PStart[static_cast<std::size_t>(p) + 1] += PStart[static_cast<std::size_t>(p)];
      std::vector<int> PSig(S.SigDPar.size()), PDeriv(S.SigDPar.size()), PWeight(S.SigDPar.size());
      std::vector<int> Next(PStart.begin(), PStart.end() - 1);
      for (int Sig = 0; Sig < NSig; ++Sig) {
        for (int k = S.SigDStart[static_cast<std::size_t>(Sig)]; k < S.SigDStart[static_cast<std::size_t>(Sig) + 1]; ++k) {
          const int At = Next[static_cast<std::size_t>(S.SigDPar[static_cast<std::size_t>(k)])]++;
          PSig[static_cast<std::size_t>(At)] = Sig;
          PDeriv[static_cast<std::size_t>(At)] = S.SigDDeriv[static_cast<std::size_t>(k)];
          PWeight[static_cast<std::size_t>(At)] = S.SigDWeight[static_cast<std::size_t>(k)];
        }
      }
      // Chunks: each parameter's entries in pieces of at most kChunk
      std::vector<int> ChunkBegin, PChunkStart(static_cast<std::size_t>(S.NPars) + 1, 0);
      for (int p = 0; p < S.NPars; ++p) {
        for (int k = PStart[static_cast<std::size_t>(p)]; k < PStart[static_cast<std::size_t>(p) + 1]; k += kChunk) ChunkBegin.push_back(k);
        PChunkStart[static_cast<std::size_t>(p) + 1] = static_cast<int>(ChunkBegin.size());
      }
      NChunks = static_cast<int>(ChunkBegin.size());
      // Chunk c covers ChunkBegin[c] to the next chunk's start, or to the end of its parameter
      std::vector<int> ChunkEdges(static_cast<std::size_t>(NChunks) + 1, 0);
      for (int p = 0; p < S.NPars; ++p) {
        for (int c = PChunkStart[static_cast<std::size_t>(p)]; c < PChunkStart[static_cast<std::size_t>(p) + 1]; ++c) {
          ChunkEdges[static_cast<std::size_t>(c)] = ChunkBegin[static_cast<std::size_t>(c)];
        }
      }
      ChunkEdges[static_cast<std::size_t>(NChunks)] = PStart[static_cast<std::size_t>(S.NPars)];
      // A chunk ends where the next begins: chunks are consecutive in entry order across parameters too
      ParSig = Upload(PSig, Bytes); ParDeriv = Upload(PDeriv, Bytes); ParWeight = Upload(PWeight, Bytes);
      ChunkBeginDev = Upload(ChunkEdges, Bytes); ParChunkStart = Upload(PChunkStart, Bytes);
      DR_CHECK(cudaMalloc(&Partial, static_cast<std::size_t>(std::max(1, NChunks)) * sizeof(double)));
      Bytes += static_cast<std::size_t>(std::max(1, NChunks)) * sizeof(double);
      DR_CHECK(cudaMalloc(&SigSum, static_cast<std::size_t>(std::max(1, NSig)) * sizeof(double)));
      Bytes += static_cast<std::size_t>(std::max(1, NSig)) * sizeof(double);
    }
    SigWStart = Upload(S.SigWStart, Bytes);
    SigWIdx = Upload(S.SigWIdx, Bytes);
    SigDStart = Upload(S.SigDStart, Bytes);
    SigDPar = Upload(S.SigDPar, Bytes);
    SigDDeriv = Upload(S.SigDDeriv, Bytes);
    SigDWeight = Upload(S.SigDWeight, Bytes);
    Data = Upload(S.Data, Bytes);
    //HW The outputs side by side, [MC | W2 | gradient], so one copy brings back whatever a reweight made.
    //Inputs and outputs are staged in pinned host memory, and everything goes on one stream with a
    //single synchronisation per reweight
    NOut = 2 * static_cast<std::size_t>(S.NBins) + static_cast<std::size_t>(S.NPars);
    DR_CHECK(cudaMalloc(&DevOut, NOut * sizeof(double)));
    MC = DevOut; W2 = DevOut + S.NBins; Grad = DevOut + 2 * S.NBins;
    Bytes += NOut * sizeof(double);
    DR_CHECK(cudaStreamCreate(&Stream));
    DR_CHECK(cudaMallocHost(&HostOut, NOut * sizeof(double)));
    DR_CHECK(cudaMallocHost(&HostPars, static_cast<std::size_t>(std::max(1, S.NPars)) * sizeof(double)));
    DR_CHECK(cudaMallocHost(&HostOther, static_cast<std::size_t>(std::max(1, S.NOther)) * sizeof(double)));
    DR_CHECK(cudaMallocHost(&HostSeg, static_cast<std::size_t>(std::max(1, S.NSplineParams)) * sizeof(short int)));
    DR_CHECK(cudaMallocHost(&HostShift, 64 * sizeof(double)));
    HostShiftCap = 64;
    GradShared = static_cast<std::size_t>(S.NPars) * sizeof(double);
    if (GradShared > 48 * 1024) {
      DR_CHECK(cudaFuncSetAttribute(GradKernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(GradShared)));
    }
    int Device = 0, SMs = 0;
    DR_CHECK(cudaGetDevice(&Device));
    DR_CHECK(cudaDeviceGetAttribute(&SMs, cudaDevAttrMultiProcessorCount, Device));
    GradBlocks = 4 * SMs;
  }

  ~DeviceReweightImpl() {
    for (void* P : {static_cast<void*>(HostOut), static_cast<void*>(HostPars), static_cast<void*>(HostOther),
                    static_cast<void*>(HostSeg), static_cast<void*>(HostShift)}) {
      if (P != nullptr) cudaFreeHost(P);
    }
    if (Stream != nullptr) cudaStreamDestroy(Stream);
    for (void* P : {static_cast<void*>(Pool), static_cast<void*>(XCoeff), static_cast<void*>(Many), static_cast<void*>(NKnots),
                    static_cast<void*>(ParamNo), static_cast<void*>(GlobalPar), static_cast<void*>(Segments),
                    static_cast<void*>(SigWStart), static_cast<void*>(SigWIdx), static_cast<void*>(SigDStart),
                    static_cast<void*>(SigDPar), static_cast<void*>(SigDDeriv), static_cast<void*>(SigDWeight),
                    static_cast<void*>(Data), static_cast<void*>(DevOut),
                    static_cast<void*>(ParSig), static_cast<void*>(ParDeriv), static_cast<void*>(ParWeight),
                    static_cast<void*>(SigSum), static_cast<void*>(ChunkBeginDev), static_cast<void*>(ParChunkStart),
                    static_cast<void*>(Partial)}) {
      if (P != nullptr) cudaFree(P);
    }
    FreeGroups();
    for (void* P : {static_cast<void*>(DKeys), static_cast<void*>(DKeysSorted), static_cast<void*>(DUnique), static_cast<void*>(DIdx),
                    static_cast<void*>(DIdxSorted), static_cast<void*>(DC1), static_cast<void*>(DC2), static_cast<void*>(DConst),
                    static_cast<void*>(DConst2), static_cast<void*>(DSig), static_cast<void*>(DBin), static_cast<void*>(DW),
                    static_cast<void*>(DBinStart), static_cast<void*>(DCounts), DTemp, static_cast<void*>(GW2Scratch), SegTemp}) {
      if (P != nullptr) cudaFree(P);
    }
    for (void* P : {static_cast<void*>(LSample), static_cast<void*>(LSig), static_cast<void*>(LCoefStart), static_cast<void*>(LCoefShift),
                    static_cast<void*>(LCutStart), static_cast<void*>(LCutSlot), static_cast<void*>(LEdgeStart), static_cast<void*>(LStride),
                    static_cast<void*>(LOffset), static_cast<void*>(LEventBin), static_cast<void*>(LConstW), static_cast<void*>(LNominal),
                    static_cast<void*>(LCoef), static_cast<void*>(LCutLo), static_cast<void*>(LCutHi), static_cast<void*>(LEdges),
                    static_cast<void*>(LShift), static_cast<void*>(LEventW)}) {
      if (P != nullptr) cudaFree(P);
    }
  }

  void FreeGroups() {
    for (void* P : {static_cast<void*>(GSig), static_cast<void*>(GConst), static_cast<void*>(GConst2),
                    static_cast<void*>(GBinStart), static_cast<void*>(GBin), static_cast<void*>(GW)}) {
      if (P != nullptr) cudaFree(P);
    }
    if (CSig == GSig) { CSig = nullptr; CConst = nullptr; CConst2 = nullptr; CBinStart = nullptr; CW = nullptr; CNGroups = 0; CBin = nullptr; }
    GSig = nullptr; GConst = nullptr; GConst2 = nullptr; GBinStart = nullptr; GBin = nullptr; GW = nullptr;
    Bytes -= GroupBytes; GroupBytes = 0; NGroups = 0;
  }

  void SetGroups(const DeviceReweightGroups& G, const std::vector<double>& NewData) {
    if (static_cast<int>(G.BinStart.size()) != St.NBins + 1) throw std::runtime_error("DeviceReweight: BinStart has the wrong length");
    if (static_cast<int>(NewData.size()) != St.NBins) throw std::runtime_error("DeviceReweight: data has the wrong length");
    FreeGroups();
    NGroups = static_cast<int>(G.Sig.size());
    std::vector<int> Bin(G.Sig.size());
    for (int b = 0; b < St.NBins; ++b) for (int g = G.BinStart[b]; g < G.BinStart[b + 1]; ++g) Bin[static_cast<std::size_t>(g)] = b;
    GSig = Upload(G.Sig, GroupBytes);
    GConst = Upload(G.Const, GroupBytes);
    GConst2 = Upload(G.Const2, GroupBytes);
    GBinStart = Upload(G.BinStart, GroupBytes);
    GBin = Upload(Bin, GroupBytes);
    DR_CHECK(cudaMalloc(&GW, (G.Sig.empty() ? 1 : G.Sig.size()) * sizeof(double)));
    GroupBytes += (G.Sig.empty() ? 1 : G.Sig.size()) * sizeof(double);
    Bytes += GroupBytes;
    DR_CHECK(cudaMemcpy(Data, NewData.data(), NewData.size() * sizeof(double), cudaMemcpyHostToDevice));
    CSig = GSig; CConst = GConst; CConst2 = GConst2; CBinStart = GBinStart; CW = GW; CNGroups = NGroups; CBin = GBin;
  }

  /// Rebuild the groups on the device from the last full pass's event bins, see GroupKeyKernel
  int BuildGroupsFromFullPass() {
    if (!HaveSweep) throw std::runtime_error("DeviceReweight: BuildGroupsFromFullPass before SetSweep");
    const int N = Sw.NLive;
    if (DKeys == nullptr) {
      const std::size_t n = static_cast<std::size_t>(std::max(1, N));
      DR_CHECK(cudaMalloc(&DKeys, n * sizeof(unsigned long long))); DR_CHECK(cudaMalloc(&DKeysSorted, n * sizeof(unsigned long long)));
      DR_CHECK(cudaMalloc(&DUnique, n * sizeof(unsigned long long)));
      DR_CHECK(cudaMalloc(&DIdx, n * sizeof(int))); DR_CHECK(cudaMalloc(&DIdxSorted, n * sizeof(int)));
      DR_CHECK(cudaMalloc(&DC1, n * sizeof(double))); DR_CHECK(cudaMalloc(&DC2, n * sizeof(double)));
      DR_CHECK(cudaMalloc(&DConst, n * sizeof(double))); DR_CHECK(cudaMalloc(&DConst2, n * sizeof(double)));
      DR_CHECK(cudaMalloc(&DSig, n * sizeof(int))); DR_CHECK(cudaMalloc(&DBin, n * sizeof(int))); DR_CHECK(cudaMalloc(&DW, n * sizeof(double)));
      DR_CHECK(cudaMalloc(&DBinStart, static_cast<std::size_t>(St.NBins + 1) * sizeof(int)));
      DR_CHECK(cudaMalloc(&DCounts, 2 * sizeof(int)));
      Bytes += n * (3 * sizeof(unsigned long long) + 4 * sizeof(int) + 5 * sizeof(double));
      // Keys only need enough bits for the largest real one plus the marker for events that fill nothing
      InvalidKey = (static_cast<unsigned long long>(St.NBins) * static_cast<unsigned long long>(NSig)) << 1;
      KeyBits = 1; while ((1ull << KeyBits) <= InvalidKey) ++KeyBits;
      std::size_t SortBytes = 0, ReduceBytes = 0;
      DR_CHECK(cub::DeviceRadixSort::SortPairs(nullptr, SortBytes, DKeys, DKeysSorted, DIdx, DIdxSorted, N, 0, KeyBits));
      DR_CHECK(cub::DeviceReduce::ReduceByKey(nullptr, ReduceBytes, DKeysSorted, DUnique, DC1, DConst, DCounts, cub::Sum(), N));
      TempBytes = std::max(SortBytes, ReduceBytes);
      DR_CHECK(cudaMalloc(&DTemp, std::max<std::size_t>(1, TempBytes)));
      Bytes += TempBytes;
    }
    DR_CHECK(cudaMemsetAsync(DCounts, 0, 2 * sizeof(int), Stream));
    GroupKeyKernel<<<(N + kThreads - 1) / kThreads, kThreads, 0, Stream>>>(N, NSig, Sw.UnderOverFlowBin, LEventBin, LSig, LConstW, InvalidKey, DKeys, DIdx, DCounts);
    DR_CHECK(cudaGetLastError());
    std::size_t Tb = TempBytes;
    DR_CHECK(cub::DeviceRadixSort::SortPairs(DTemp, Tb, DKeys, DKeysSorted, DIdx, DIdxSorted, N, 0, KeyBits, Stream));
    int NValid = 0;
    DR_CHECK(cudaMemcpyAsync(&NValid, DCounts, sizeof(int), cudaMemcpyDeviceToHost, Stream));
    DR_CHECK(cudaStreamSynchronize(Stream));
    int NG = 0;
    if (NValid > 0) {
      GroupGatherKernel<<<(NValid + kThreads - 1) / kThreads, kThreads, 0, Stream>>>(NValid, DIdxSorted, LConstW, DC1, DC2);
      DR_CHECK(cudaGetLastError());
      Tb = TempBytes;
      DR_CHECK(cub::DeviceReduce::ReduceByKey(DTemp, Tb, DKeysSorted, DUnique, DC1, DConst, DCounts + 1, cub::Sum(), NValid, Stream));
      Tb = TempBytes;
      DR_CHECK(cub::DeviceReduce::ReduceByKey(DTemp, Tb, DKeysSorted, DUnique, DC2, DConst2, DCounts + 1, cub::Sum(), NValid, Stream));
      DR_CHECK(cudaMemcpyAsync(&NG, DCounts + 1, sizeof(int), cudaMemcpyDeviceToHost, Stream));
      DR_CHECK(cudaStreamSynchronize(Stream));
      GroupDecodeKernel<<<(NG + kThreads - 1) / kThreads, kThreads, 0, Stream>>>(NG, NSig, DUnique, DSig, DBin);
      DR_CHECK(cudaGetLastError());
    }
    GroupBinStartKernel<<<(St.NBins + 1 + kThreads - 1) / kThreads, kThreads, 0, Stream>>>(St.NBins, NG, DBin, DBinStart);
    DR_CHECK(cudaGetLastError());
    CSig = DSig; CConst = DConst; CConst2 = DConst2; CBinStart = DBinStart; CW = DW; CNGroups = NG; CBin = DBin;
    return NG;
  }

  void SetData(const std::vector<double>& NewData) {
    if (static_cast<int>(NewData.size()) != St.NBins) throw std::runtime_error("DeviceReweight: data has the wrong length");
    DR_CHECK(cudaMemcpy(Data, NewData.data(), NewData.size() * sizeof(double), cudaMemcpyHostToDevice));
  }

  void FillPool(const double* ParValues, const short int* Seg, const double* Other) {
    std::memcpy(HostPars, ParValues, static_cast<std::size_t>(St.NPars) * sizeof(double));
    DR_CHECK(cudaMemcpyAsync(Pool, HostPars, static_cast<std::size_t>(St.NPars) * sizeof(double), cudaMemcpyHostToDevice, Stream));
    if (St.NOther > 0) {
      std::memcpy(HostOther, Other, static_cast<std::size_t>(St.NOther) * sizeof(double));
      DR_CHECK(cudaMemcpyAsync(Pool + OtherOffset, HostOther, static_cast<std::size_t>(St.NOther) * sizeof(double), cudaMemcpyHostToDevice, Stream));
    }
    if (St.NSplines > 0) {
      std::memcpy(HostSeg, Seg, static_cast<std::size_t>(St.NSplineParams) * sizeof(short int));
      DR_CHECK(cudaMemcpyAsync(Segments, HostSeg, static_cast<std::size_t>(St.NSplineParams) * sizeof(short int), cudaMemcpyHostToDevice, Stream));
      SplineKernel<<<(St.NSplines + kThreads - 1) / kThreads, kThreads, 0, Stream>>>(
          St.NSplines, St.MaxKnots, ParamNo, NKnots, Many, XCoeff, Segments, GlobalPar, Pool, St.NPars, St.NPars + St.NSplines);
      DR_CHECK(cudaGetLastError());
    }
  }

  void Run(const double* ParValues, const short int* Seg, const double* Other, const bool FillW2, const bool Gradient,
           double* MCOut, double* W2Out, double* GradOut) {
    if (CBinStart == nullptr) throw std::runtime_error("DeviceReweight: Run before any groups");
    FillPool(ParValues, Seg, Other);
    if (Gradient) DR_CHECK(cudaMemsetAsync(SigSum, 0, static_cast<std::size_t>(std::max(1, NSig)) * sizeof(double), Stream));
    EnsureGroupScratch();
    if (CNGroups > 0) {
      GroupWeightKernel<<<(CNGroups + kThreads - 1) / kThreads, kThreads, 0, Stream>>>(CNGroups, CSig, CConst, CConst2, SigWStart,
                                                                                         SigWIdx, Pool, CW, GW2Scratch, FillW2);
      DR_CHECK(cudaGetLastError());
    }
    std::size_t Tb = SegTempBytes;
    DR_CHECK(cub::DeviceSegmentedReduce::Sum(SegTemp, Tb, CW, MC, St.NBins, CBinStart, CBinStart + 1, Stream));
    if (FillW2) {
      Tb = SegTempBytes;
      DR_CHECK(cub::DeviceSegmentedReduce::Sum(SegTemp, Tb, GW2Scratch, W2, St.NBins, CBinStart, CBinStart + 1, Stream));
    }
    if (Gradient && CNGroups > 0) {
      GroupSigSumKernel<<<(CNGroups + kThreads - 1) / kThreads, kThreads, 0, Stream>>>(CNGroups, CBin, CSig, CW, MC, Data,
                                                                                         St.LowMCBound, SigSum);
      DR_CHECK(cudaGetLastError());
    }
    if (Gradient) {
      if (NChunks > 0) {
        ParChunkKernel<<<NChunks, kThreads, 0, Stream>>>(ChunkBeginDev, ParSig, ParDeriv, ParWeight, SigSum, Pool, Partial);
        DR_CHECK(cudaGetLastError());
      }
      ParSumKernel<<<(St.NPars + kThreads - 1) / kThreads, kThreads, 0, Stream>>>(St.NPars, ParChunkStart, Partial, Grad);
      DR_CHECK(cudaGetLastError());
    }
    CopyBack(FillW2, Gradient, MCOut, W2Out, GradOut);
  }

  /// Per group W2 scratch and the segmented reduction's working space, sized for the largest group count
  /// so far (at most the live event count)
  void EnsureGroupScratch() {
    const int Need = std::max(1, CNGroups);
    if (Need > GW2Cap) {
      if (GW2Scratch != nullptr) cudaFree(GW2Scratch);
      GW2Cap = std::max(Need, std::max(1, Sw.NLive));
      DR_CHECK(cudaMalloc(&GW2Scratch, static_cast<std::size_t>(GW2Cap) * sizeof(double)));
    }
    if (SegTemp == nullptr) {
      DR_CHECK(cub::DeviceSegmentedReduce::Sum(nullptr, SegTempBytes, CW, MC, St.NBins, CBinStart, CBinStart + 1, Stream));
      DR_CHECK(cudaMalloc(&SegTemp, std::max<std::size_t>(1, SegTempBytes)));
    }
  }

  /// One copy of whatever this reweight made, [MC | W2 | gradient] or less, then the only wait
  void CopyBack(const bool FillW2, const bool Gradient, double* MCOut, double* W2Out, double* GradOut) {
    const std::size_t NB = static_cast<std::size_t>(St.NBins);
    const std::size_t N = Gradient ? NOut : (FillW2 ? 2 * NB : NB);
    DR_CHECK(cudaMemcpyAsync(HostOut, DevOut, N * sizeof(double), cudaMemcpyDeviceToHost, Stream));
    DR_CHECK(cudaStreamSynchronize(Stream));
    std::memcpy(MCOut, HostOut, NB * sizeof(double));
    if (FillW2) std::memcpy(W2Out, HostOut + NB, NB * sizeof(double));
    if (Gradient) std::memcpy(GradOut, HostOut + 2 * NB, static_cast<std::size_t>(St.NPars) * sizeof(double));
  }

  void SetSweep(const DeviceSweepStatic& W) {
    if (W.NSlots > 8) throw std::runtime_error("DeviceReweight: more than 8 shifted variables per event");
    if (static_cast<int>(W.EdgeStart.size()) < 2) throw std::runtime_error("DeviceReweight: no binning for the sweep");
    Sw = W;
    Sw.Nominal.clear(); Sw.Coef.clear(); Sw.CoefShift.clear(); Sw.CoefStart.clear();  // only the counts are kept on the host
    LSample = Upload(W.Sample, Bytes); LSig = Upload(W.Sig, Bytes); LConstW = Upload(W.ConstW, Bytes);
    LNominal = Upload(W.Nominal, Bytes); LCoefStart = Upload(W.CoefStart, Bytes);
    LCoefShift = Upload(W.CoefShift, Bytes); LCoef = Upload(W.Coef, Bytes);
    LCutStart = Upload(W.CutStart, Bytes); LCutSlot = Upload(W.CutSlot, Bytes);
    LCutLo = Upload(W.CutLo, Bytes); LCutHi = Upload(W.CutHi, Bytes);
    LEdgeStart = Upload(W.EdgeStart, Bytes); LEdges = Upload(W.Edges, Bytes);
    LStride = Upload(W.Stride, Bytes); LOffset = Upload(W.Offset, Bytes);
    DR_CHECK(cudaMalloc(&LShift, static_cast<std::size_t>(std::max(1, W.NShifts)) * sizeof(double)));
    DR_CHECK(cudaMalloc(&LEventBin, static_cast<std::size_t>(std::max(1, W.NLive)) * sizeof(int)));
    DR_CHECK(cudaMalloc(&LEventW, static_cast<std::size_t>(std::max(1, W.NLive)) * sizeof(double)));
    Bytes += static_cast<std::size_t>(std::max(1, W.NLive)) * (sizeof(int) + sizeof(double)) + static_cast<std::size_t>(std::max(1, W.NShifts)) * sizeof(double);
    FullShared = static_cast<std::size_t>(2 * St.NBins) * sizeof(double);
    if (FullShared > 48 * 1024) {
      DR_CHECK(cudaFuncSetAttribute(FullPassKernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(FullShared)));
    }
    HaveSweep = true;
  }

  void RunFull(const double* ParValues, const short int* Seg, const double* Other, const double* ShiftValues,
               const bool FillW2, const bool Gradient, double* MCOut, double* W2Out, double* GradOut) {
    if (!HaveSweep) throw std::runtime_error("DeviceReweight: RunFull before SetSweep");
    FillPool(ParValues, Seg, Other);
    if (Sw.NShifts > 0) {
      if (Sw.NShifts > HostShiftCap) {
        cudaFreeHost(HostShift);
        HostShiftCap = Sw.NShifts;
        DR_CHECK(cudaMallocHost(&HostShift, static_cast<std::size_t>(HostShiftCap) * sizeof(double)));
      }
      std::memcpy(HostShift, ShiftValues, static_cast<std::size_t>(Sw.NShifts) * sizeof(double));
      DR_CHECK(cudaMemcpyAsync(LShift, HostShift, static_cast<std::size_t>(Sw.NShifts) * sizeof(double), cudaMemcpyHostToDevice, Stream));
    }
    DR_CHECK(cudaMemsetAsync(MC, 0, static_cast<std::size_t>(St.NBins) * sizeof(double), Stream));
    if (FillW2) DR_CHECK(cudaMemsetAsync(W2, 0, static_cast<std::size_t>(St.NBins) * sizeof(double), Stream));
    FullPassKernel<<<GradBlocks, kThreads, FullShared, Stream>>>(
        Sw.NLive, Sw.NSlots, St.NBins, LSample, LSig, LConstW, LNominal, LCoefStart, LCoefShift, LCoef, LShift,
        LCutStart, LCutSlot, LCutLo, LCutHi, LEdgeStart, LEdges, LStride, LOffset, SigWStart, SigWIdx, Pool,
        Sw.UnderOverFlowBin, Sw.NotSelected, FillW2, LEventBin, LEventW, MC, W2);
    DR_CHECK(cudaGetLastError());
    if (Gradient) {
      // Each event as a group of one: its bin, its weight (zero when it doesn't fill), its signature
      DR_CHECK(cudaMemsetAsync(Grad, 0, static_cast<std::size_t>(St.NPars) * sizeof(double), Stream));
      GradKernel<<<GradBlocks, kThreads, GradShared, Stream>>>(Sw.NLive, St.NPars, LEventBin, LSig, LEventW, SigDStart, SigDPar,
                                                               SigDDeriv, SigDWeight, Pool, MC, Data, St.LowMCBound, Grad);
      DR_CHECK(cudaGetLastError());
    }
    CopyBack(FillW2, Gradient, MCOut, W2Out, GradOut);
  }

  void GetEventBins(int* Out) const {
    if (!HaveSweep) throw std::runtime_error("DeviceReweight: GetEventBins before SetSweep");
    DR_CHECK(cudaStreamSynchronize(Stream));
    DR_CHECK(cudaMemcpy(Out, LEventBin, static_cast<std::size_t>(Sw.NLive) * sizeof(int), cudaMemcpyDeviceToHost));
  }

  bool HaveSweep = false;
  std::size_t Bytes = 0;

 private:
  DeviceReweightStatic St;
  int NPool = 0, OtherOffset = 0, NGroups = 0, GradBlocks = 0;
  std::size_t GradShared = 0, GroupBytes = 0;
  double* Pool = nullptr;
  double* XCoeff = nullptr;
  double* Many = nullptr;
  unsigned int* NKnots = nullptr;
  short int* ParamNo = nullptr;
  int* GlobalPar = nullptr;
  short int* Segments = nullptr;
  int *SigWStart = nullptr, *SigWIdx = nullptr, *SigDStart = nullptr, *SigDPar = nullptr, *SigDDeriv = nullptr, *SigDWeight = nullptr;
  double *Data = nullptr, *MC = nullptr, *W2 = nullptr, *Grad = nullptr, *DevOut = nullptr;
  std::size_t NOut = 0;
  cudaStream_t Stream = nullptr;
  double *HostOut = nullptr, *HostPars = nullptr, *HostOther = nullptr, *HostShift = nullptr;
  short int* HostSeg = nullptr;
  int HostShiftCap = 0;
  int *GSig = nullptr, *GBinStart = nullptr, *GBin = nullptr;
  // Groups the next Run uses: the host's (SetGroups) or the device's own (BuildGroupsFromFullPass)
  int *CSig = nullptr, *CBinStart = nullptr, *CBin = nullptr, CNGroups = 0;
  double* GW2Scratch = nullptr;
  int GW2Cap = 0;
  void* SegTemp = nullptr;
  std::size_t SegTempBytes = 0;
  double *CConst = nullptr, *CConst2 = nullptr, *CW = nullptr;
  // Device built groups and their working space
  unsigned long long *DKeys = nullptr, *DKeysSorted = nullptr, *DUnique = nullptr, InvalidKey = 0;
  int *DIdx = nullptr, *DIdxSorted = nullptr, *DSig = nullptr, *DBin = nullptr, *DBinStart = nullptr, *DCounts = nullptr;
  double *DC1 = nullptr, *DC2 = nullptr, *DConst = nullptr, *DConst2 = nullptr, *DW = nullptr;
  void* DTemp = nullptr;
  std::size_t TempBytes = 0;
  int KeyBits = 64;
  int NSig = 0;
  int *ParSig = nullptr, *ParDeriv = nullptr, *ParWeight = nullptr, *ChunkBeginDev = nullptr, *ParChunkStart = nullptr;
  int NChunks = 0;
  double* Partial = nullptr;
  double* SigSum = nullptr;
  double *GConst = nullptr, *GConst2 = nullptr, *GW = nullptr;
  // Full pass
  DeviceSweepStatic Sw;
  std::size_t FullShared = 0;
  int *LSample = nullptr, *LSig = nullptr, *LCoefStart = nullptr, *LCoefShift = nullptr, *LCutStart = nullptr, *LCutSlot = nullptr;
  int *LEdgeStart = nullptr, *LStride = nullptr, *LOffset = nullptr, *LEventBin = nullptr;
  double *LConstW = nullptr, *LNominal = nullptr, *LCoef = nullptr, *LCutLo = nullptr, *LCutHi = nullptr, *LEdges = nullptr;
  double *LShift = nullptr, *LEventW = nullptr;
};

DeviceReweight::DeviceReweight(const DeviceReweightStatic& Static) : Impl(new DeviceReweightImpl(Static)) {}
DeviceReweight::~DeviceReweight() = default;
void DeviceReweight::SetGroups(const DeviceReweightGroups& Groups, const std::vector<double>& Data) { Impl->SetGroups(Groups, Data); }
void DeviceReweight::Run(const double* ParValues, const short int* Segments, const double* Other, const bool FillW2,
                         const bool Gradient, double* MCOut, double* W2Out, double* GradOut) {
  Impl->Run(ParValues, Segments, Other, FillW2, Gradient, MCOut, W2Out, GradOut);
}
std::size_t DeviceReweight::DeviceBytes() const { return Impl->Bytes; }
void DeviceReweight::SetSweep(const DeviceSweepStatic& Sweep) { Impl->SetSweep(Sweep); }
void DeviceReweight::SetData(const std::vector<double>& NewData) { Impl->SetData(NewData); }
void DeviceReweight::RunFull(const double* ParValues, const short int* Segments, const double* Other, const double* ShiftValues,
                             const bool FillW2, const bool Gradient, double* MCOut, double* W2Out, double* GradOut) {
  Impl->RunFull(ParValues, Segments, Other, ShiftValues, FillW2, Gradient, MCOut, W2Out, GradOut);
}
void DeviceReweight::GetEventBins(int* Out) const { Impl->GetEventBins(Out); }
int DeviceReweight::BuildGroupsFromFullPass() { return Impl->BuildGroupsFromFullPass(); }
