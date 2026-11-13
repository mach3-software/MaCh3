// SampleHandlerBase's reweight on the GPU. See DeviceReweight.h for what is computed and why;
// this file is laid out as:
//   1. small helpers (error checks, launch sizes, an owning device array)
//   2. the kernels, in the order a reweight runs them
//   3. DeviceReweightImpl, which owns the device memory and runs the kernels
//   4. the thin DeviceReweight wrapper
//
// Every kernel does the same arithmetic as the CPU code it replaces, in double precision, so the GPU
// and CPU results agree to rounding error. Where the order of a sum differs (parallel reductions),
// that is the only difference.
#include "Samples/DeviceReweight.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

#include <cub/block/block_reduce.cuh>
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_reduce.cuh>
#include <cub/device/device_segmented_reduce.cuh>
#include <cuda_runtime.h>

namespace {

// ================================================================================================
// 1. Helpers
// ================================================================================================

/// Throw with the call, the line and CUDA's message if a CUDA call failed
void Check(const cudaError_t Err, const char* Call, const int Line) {
  if (Err != cudaSuccess) {
    throw std::runtime_error(std::string("DeviceReweight: ") + Call + " failed at line " + std::to_string(Line) + ": " +
                             cudaGetErrorString(Err));
  }
}
#define DR_CHECK(x) Check((x), #x, __LINE__)

/// Threads per block for every kernel
constexpr int kThreads = 256;

/// Blocks needed for one thread per item
int BlocksFor(const int NItems) {
  return (NItems + kThreads - 1) / kThreads;
}

/// An array in device memory that frees itself. Every allocation adds to a shared byte count, for
/// the "MB on the device" log line
template <typename T>
class DeviceArray {
 public:
  DeviceArray() = default;
  ~DeviceArray() { Free(); }
  DeviceArray(const DeviceArray&) = delete;
  DeviceArray& operator=(const DeviceArray&) = delete;

  /// Allocate room for Count elements (at least one, so the pointer is never null)
  void Allocate(const std::size_t Count, std::size_t& ByteCount) {
    Free();
    Size = std::max<std::size_t>(1, Count);
    DR_CHECK(cudaMalloc(&Ptr, Size * sizeof(T)));
    ByteCount += Size * sizeof(T);
  }

  /// Allocate and copy a host vector up
  void Upload(const std::vector<T>& Host, std::size_t& ByteCount) {
    Allocate(Host.size(), ByteCount);
    if (!Host.empty()) DR_CHECK(cudaMemcpy(Ptr, Host.data(), Host.size() * sizeof(T), cudaMemcpyHostToDevice));
  }

  void Free() {
    if (Ptr != nullptr) cudaFree(Ptr);
    Ptr = nullptr;
    Size = 0;
  }

  T* get() const { return Ptr; }
  std::size_t size() const { return Size; }
  bool empty() const { return Ptr == nullptr; }

 private:
  T* Ptr = nullptr;
  std::size_t Size = 0;
};

/// An array in pinned (page locked) host memory, which asynchronous copies need
template <typename T>
class PinnedArray {
 public:
  PinnedArray() = default;
  ~PinnedArray() { if (Ptr != nullptr) cudaFreeHost(Ptr); }
  PinnedArray(const PinnedArray&) = delete;
  PinnedArray& operator=(const PinnedArray&) = delete;

  /// Make room for at least Count elements, keeping the old allocation if it is big enough
  void Reserve(const std::size_t Count) {
    if (Count <= Size && Ptr != nullptr) return;
    if (Ptr != nullptr) cudaFreeHost(Ptr);
    Size = std::max<std::size_t>(1, Count);
    DR_CHECK(cudaMallocHost(&Ptr, Size * sizeof(T)));
  }

  T* get() const { return Ptr; }

 private:
  T* Ptr = nullptr;
  std::size_t Size = 0;
};

// ================================================================================================
// 2. Kernels
// ================================================================================================

// ------------------------------------------------------------------------------------------------
// Splines (both passes)
// ------------------------------------------------------------------------------------------------

/// One thread per spline: its weight and d(weight)/d(parameter) at the current parameter value,
/// written into the pool. Same cubic as BinnedSplineHandler::CalcSplineWeights:
///   weight = y + b*dx + c*dx^2 + d*dx^3,  derivative = b + 2c*dx + 3d*dx^2,  dx = x - knot
__global__ void EvaluateSplines(const int NSplines, const int MaxKnots,
                                const short int* __restrict__ SplineParam,        // spline parameter of each spline
                                const unsigned int* __restrict__ CoeffStart,       // first segment of each spline in Coeffs
                                const double* __restrict__ Coeffs,                 // y, b, c, d per segment
                                const double* __restrict__ KnotX,                  // knot positions, MaxKnots per spline parameter
                                const short int* __restrict__ Segment,             // current segment of each spline parameter
                                const int* __restrict__ GlobalParam,               // pool index of each spline parameter's value
                                double* __restrict__ Pool, const int WeightOffset, const int DerivOffset) {
  const int s = blockIdx.x * blockDim.x + threadIdx.x;
  if (s >= NSplines) return;

  const int Param = SplineParam[s];
  const int Seg = Segment[Param];
  const double* C = Coeffs + static_cast<size_t>(CoeffStart[s] + Seg) * 4;
  const double y = C[0];
  const double b = C[1];
  const double c = C[2];
  const double d = C[3];
  const double dx = Pool[GlobalParam[Param]] - KnotX[Param * MaxKnots + Seg];

  Pool[WeightOffset + s] = fma(dx, fma(dx, fma(dx, d, c), b), y);
  Pool[DerivOffset + s] = fma(dx, fma(dx, 3. * d, 2. * c), b);
}

// ------------------------------------------------------------------------------------------------
// Leapfrog pass (Run): groups -> bins -> gradient
// ------------------------------------------------------------------------------------------------

/// Step 1. One thread per group: its weight, which is the group's summed constant weights times the
/// product of its signature's varying weights, and the group's contribution to W2.
/// As in the CPU loop, a group whose weight isn't positive contributes nothing (gradient included).
__global__ void WeighGroups(const int NGroups, const int* __restrict__ GroupSig, const double* __restrict__ GroupConst,
                            const double* __restrict__ GroupConst2, const int* __restrict__ SigWeightStart,
                            const int* __restrict__ SigWeightIdx, const double* __restrict__ Pool,
                            double* __restrict__ GroupWeight, double* __restrict__ GroupW2, const bool FillW2) {
  const int g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= NGroups) return;

  const int Sig = GroupSig[g];
  double Shared = 1.;
  for (int k = SigWeightStart[Sig]; k < SigWeightStart[Sig + 1]; ++k) {
    Shared *= Pool[SigWeightIdx[k]];
  }

  double Weight = GroupConst[g] * Shared;
  double W2 = 0.;
  if (!(Weight > 0.)) {
    Weight = 0.;
  } else if (FillW2) {
    W2 = GroupConst2[g] * Shared * Shared;
  }
  GroupWeight[g] = Weight;
  if (FillW2) GroupW2[g] = W2;
}

// Step 2 sums each bin's group weights with CUB's segmented reduction (see Run).

/// Step 3, gradient only. One thread per group: add (1 - data/MC) * group weight to its signature's
/// total. That factor is d(-lnL)/d(MC) for the Poisson likelihood, the same pre-factor
/// SampleHandlerBase::GetGradLikelihood uses, and bins below the low-MC bound are skipped as there.
__global__ void SumBySignature(const int NGroups, const int* __restrict__ GroupBin, const int* __restrict__ GroupSig,
                               const double* __restrict__ GroupWeight, const double* __restrict__ MC,
                               const double* __restrict__ Data, const double LowMC, double* __restrict__ SigTotal) {
  const int g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= NGroups) return;

  const double Weight = GroupWeight[g];
  if (Weight == 0.) return;
  const int Bin = GroupBin[g];
  const double Mu = MC[Bin];
  if (Mu < LowMC) return;

  atomicAdd(&SigTotal[GroupSig[g]], (1. - Data[Bin] / Mu) * Weight);
}

/// The derivative entries of one parameter are cut into chunks of at most this many, so a
/// parameter that every event depends on (a fiducial volume norm) doesn't leave one block working alone
constexpr int kChunk = 256;

/// Step 4, gradient only. One block per chunk of a parameter's derivative entries: the sum over the
/// chunk of SigTotal * d(w)/d(par) / w, where w is the weight the entry differentiates.
/// Since the event weight is a product, d(weight)/d(par) = weight * sum over its factors of d(w)/d(par) / w
__global__ void SumChunks(const int* __restrict__ ChunkBegin, const int* __restrict__ EntrySig,
                          const int* __restrict__ EntryDeriv, const int* __restrict__ EntryWeight,
                          const double* __restrict__ SigTotal, const double* __restrict__ Pool,
                          double* __restrict__ ChunkTotal) {
  using Reduce = cub::BlockReduce<double, kThreads>;
  __shared__ typename Reduce::TempStorage Temp;

  const int c = blockIdx.x;
  double Sum = 0.;
  for (int k = ChunkBegin[c] + threadIdx.x; k < ChunkBegin[c + 1]; k += blockDim.x) {
    const double A = SigTotal[EntrySig[k]];
    if (A != 0.) Sum += A * Pool[EntryDeriv[k]] / Pool[EntryWeight[k]];
  }
  const double Total = Reduce(Temp).Sum(Sum);
  if (threadIdx.x == 0) ChunkTotal[c] = Total;
}

/// Step 5, gradient only. One thread per parameter: add its chunks, in order, into d(-lnL)/d(par)
__global__ void SumParameters(const int NPars, const int* __restrict__ ParChunkStart, const double* __restrict__ ChunkTotal,
                              double* __restrict__ Grad) {
  const int Par = blockIdx.x * blockDim.x + threadIdx.x;
  if (Par >= NPars) return;

  double Sum = 0.;
  for (int c = ParChunkStart[Par]; c < ParChunkStart[Par + 1]; ++c) Sum += ChunkTotal[c];
  Grad[Par] = Sum;
}

// ------------------------------------------------------------------------------------------------
// Full pass (RunFull): shifts, cuts, bins and weights of every live event
// ------------------------------------------------------------------------------------------------

/// Upper limit on shifted variables per event (the binning variable plus the shifted cuts)
constexpr int kMaxVariables = 8;

/// Global bin of value x in a uniform 1D binning given by Edges[0..NBins], or UnderOverFlow outside it.
/// Same answer as SampleHandlerBase::FindGlobalBinFast: upper_bound over the edges, minus one
__device__ int FindBin(const double x, const double* Edges, const int NBins, const int Stride, const int Offset,
                       const int UnderOverFlow) {
  if (x < Edges[0] || x >= Edges[NBins]) return UnderOverFlow;
  int Lo = 0;
  int Hi = NBins + 1;
  while (Lo < Hi) {
    const int Mid = (Lo + Hi) / 2;
    if (Edges[Mid] <= x) Lo = Mid + 1;
    else Hi = Mid;
  }
  return (Lo - 1) * Stride + Offset;
}

/// One thread per live event (grid stride), the CPU full pass:
///  - each shifted variable = nominal + sum over its shifts of (shift value * coefficient), added in the
///    CPU's order with separate multiply and add (no fused multiply-add, as the CPU does it)
///  - the shifted cuts of its sample (accepted if Lo <= x < Hi, as IsEventSelectedFast)
///  - its bin (FindBin) and weight (constant weights times its signature's varying weights)
/// Bin totals are summed per block in shared memory, then added to the global arrays once per block.
/// Each event's bin and weight are kept for the gradient and for building groups afterwards
__global__ void FullPass(const int NLive, const int NVars, const int NBins,
                         const int* __restrict__ Sample, const int* __restrict__ Sig, const double* __restrict__ ConstWeight,
                         const double* __restrict__ Nominal, const int* __restrict__ CoefStart,
                         const int* __restrict__ CoefShift, const double* __restrict__ Coef,
                         const double* __restrict__ ShiftValue,
                         const int* __restrict__ CutStart, const int* __restrict__ CutVar,
                         const double* __restrict__ CutLo, const double* __restrict__ CutHi,
                         const int* __restrict__ EdgeStart, const double* __restrict__ Edges,
                         const int* __restrict__ Stride, const int* __restrict__ Offset,
                         const int* __restrict__ SigWeightStart, const int* __restrict__ SigWeightIdx,
                         const double* __restrict__ Pool, const int UnderOverFlow, const int NotSelected,
                         const bool FillW2, int* __restrict__ EventBin, double* __restrict__ EventWeight,
                         double* __restrict__ MC, double* __restrict__ W2) {
  extern __shared__ double BlockHist[];  // [MC per bin | W2 per bin]
  const int NHist = FillW2 ? 2 * NBins : NBins;
  for (int b = threadIdx.x; b < NHist; b += blockDim.x) BlockHist[b] = 0.;
  __syncthreads();

  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < NLive; i += gridDim.x * blockDim.x) {
    // Shifted variables: 0 is the binning variable, 1.. the shifted cut variables
    double Var[kMaxVariables];
    for (int v = 0; v < NVars; ++v) {
      const int Cell = i * NVars + v;
      double x = Nominal[Cell];
      for (int k = CoefStart[Cell]; k < CoefStart[Cell + 1]; ++k) {
        x = __dadd_rn(x, __dmul_rn(ShiftValue[CoefShift[k]], Coef[k]));
      }
      Var[v] = x;
    }

    // Cuts
    const int S = Sample[i];
    bool Selected = true;
    for (int c = CutStart[S]; c < CutStart[S + 1]; ++c) {
      const double x = Var[CutVar[c]];
      if (x < CutLo[c] || x >= CutHi[c]) {
        Selected = false;
        break;
      }
    }

    // Bin and weight
    int Bin = NotSelected;
    double Weight = 0.;
    if (Selected) {
      const int NSampleBins = EdgeStart[S + 1] - EdgeStart[S] - 1;
      Bin = FindBin(Var[0], Edges + EdgeStart[S], NSampleBins, Stride[S], Offset[S], UnderOverFlow);
      if (Bin > UnderOverFlow) {
        const int G = Sig[i];
        double Shared = 1.;
        for (int k = SigWeightStart[G]; k < SigWeightStart[G + 1]; ++k) Shared *= Pool[SigWeightIdx[k]];
        Weight = ConstWeight[i] * Shared;
        if (Weight > 0.) {
          atomicAdd(&BlockHist[Bin], Weight);
          if (FillW2) atomicAdd(&BlockHist[NBins + Bin], Weight * Weight);
        } else {
          Weight = 0.;
        }
      }
    }
    EventBin[i] = Bin;
    EventWeight[i] = Weight;
  }

  __syncthreads();
  for (int b = threadIdx.x; b < NBins; b += blockDim.x) {
    if (BlockHist[b] != 0.) atomicAdd(&MC[b], BlockHist[b]);
    if (FillW2 && BlockHist[NBins + b] != 0.) atomicAdd(&W2[b], BlockHist[NBins + b]);
  }
}

/// Gradient of a full pass (only asked for when a fit starts; the random walk sweeps pause gradients).
/// Grid stride over events: d(-lnL)/d(par) += (1 - data/MC) * event weight * d(w)/d(par) / w.
/// Summed per block in shared memory first, so the global array only takes NPars atomics per block
__global__ void FullPassGradient(const int NEvents, const int NPars, const int* __restrict__ EventBin,
                                 const int* __restrict__ EventSig, const double* __restrict__ EventWeight,
                                 const int* __restrict__ SigDerivStart, const int* __restrict__ SigDerivPar,
                                 const int* __restrict__ SigDerivIdx, const int* __restrict__ SigDerivWeightIdx,
                                 const double* __restrict__ Pool, const double* __restrict__ MC,
                                 const double* __restrict__ Data, const double LowMC, double* __restrict__ Grad) {
  extern __shared__ double BlockGrad[];
  for (int p = threadIdx.x; p < NPars; p += blockDim.x) BlockGrad[p] = 0.;
  __syncthreads();

  for (int e = blockIdx.x * blockDim.x + threadIdx.x; e < NEvents; e += gridDim.x * blockDim.x) {
    const double Weight = EventWeight[e];
    if (Weight == 0.) continue;
    const int Bin = EventBin[e];
    const double Mu = MC[Bin];
    if (Mu < LowMC) continue;
    const double Factor = (1. - Data[Bin] / Mu) * Weight;
    const int Sig = EventSig[e];
    for (int k = SigDerivStart[Sig]; k < SigDerivStart[Sig + 1]; ++k) {
      atomicAdd(&BlockGrad[SigDerivPar[k]], Factor * Pool[SigDerivIdx[k]] / Pool[SigDerivWeightIdx[k]]);
    }
  }

  __syncthreads();
  for (int p = threadIdx.x; p < NPars; p += blockDim.x) {
    if (BlockGrad[p] != 0.) atomicAdd(&Grad[p], BlockGrad[p]);
  }
}

// ------------------------------------------------------------------------------------------------
// Building groups from the full pass's event bins (BuildGroupsFromFullPass)
// ------------------------------------------------------------------------------------------------
// The CPU groups the events of each bin by signature, positives and negatives of the constant weight
// apart, in the order (bin, signature, sign). Here: give every event that key, sort by it, then sum each
// run of equal keys. Events that fill nothing get a key larger than any real one, so they sort last.

/// Step 1. One thread per live event: its sort key, and its index to carry through the sort
__global__ void MakeGroupKeys(const int NLive, const int NSig, const int UnderOverFlow, const int* __restrict__ EventBin,
                              const int* __restrict__ Sig, const double* __restrict__ ConstWeight,
                              const unsigned long long NoGroupKey, unsigned long long* __restrict__ Key,
                              int* __restrict__ Index, int* __restrict__ NWithGroup) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= NLive) return;

  const int Bin = EventBin[i];
  const double C = ConstWeight[i];
  const bool InGroup = (Bin > UnderOverFlow) && (C != 0.);
  if (InGroup) {
    const unsigned long long BinSig = static_cast<unsigned long long>(Bin) * static_cast<unsigned long long>(NSig) +
                                      static_cast<unsigned long long>(Sig[i]);
    const unsigned long long Negative = (C < 0.) ? 1ull : 0ull;
    Key[i] = (BinSig << 1) | Negative;
    atomicAdd(NWithGroup, 1);
  } else {
    Key[i] = NoGroupKey;
  }
  Index[i] = i;
}

// Step 2 sorts the keys (CUB radix sort, which keeps equal keys in their original order).

/// Step 3. One thread per sorted event: its constant weight and the square, in sorted order
__global__ void GatherConstWeights(const int N, const int* __restrict__ SortedIndex, const double* __restrict__ ConstWeight,
                                   double* __restrict__ C, double* __restrict__ C2) {
  const int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= N) return;
  const double Value = ConstWeight[SortedIndex[j]];
  C[j] = Value;
  C2[j] = Value * Value;
}

// Step 4 sums each run of equal keys (CUB reduce by key), giving one group per run.

/// Step 5. One thread per group: its signature and bin back from its key
__global__ void DecodeGroupKeys(const int NGroups, const int NSig, const unsigned long long* __restrict__ Key,
                                int* __restrict__ GroupSig, int* __restrict__ GroupBin) {
  const int g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= NGroups) return;
  const unsigned long long BinSig = Key[g] >> 1;
  GroupSig[g] = static_cast<int>(BinSig % static_cast<unsigned long long>(NSig));
  GroupBin[g] = static_cast<int>(BinSig / static_cast<unsigned long long>(NSig));
}

/// Step 6. One thread per bin (and one past the last): where its groups start, by binary search over
/// the group bins, which the sort left in order
__global__ void FindBinStarts(const int NBins, const int NGroups, const int* __restrict__ GroupBin, int* __restrict__ BinStart) {
  const int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b > NBins) return;
  int Lo = 0;
  int Hi = NGroups;
  while (Lo < Hi) {
    const int Mid = (Lo + Hi) / 2;
    if (GroupBin[Mid] < b) Lo = Mid + 1;
    else Hi = Mid;
  }
  BinStart[b] = Lo;
}

}  // namespace

// ================================================================================================
// 3. DeviceReweightImpl
// ================================================================================================

class DeviceReweightImpl {
 public:
  // ----------------------------------------------------------------------------------------------
  // Setup
  // ----------------------------------------------------------------------------------------------

  explicit DeviceReweightImpl(const DeviceReweightStatic& Static) : Info(Static) {
    NSig = static_cast<int>(Static.SigWeightStart.size()) - 1;
    SetUpPool();
    if (Static.NSplines > 0) UploadSplines();
    UploadSignatures();
    BuildGradientEntries();

    Data.Upload(Static.Data, Bytes);

    // The outputs side by side, [MC | W2 | gradient], so one copy brings back whatever a reweight made
    NOut = 2 * static_cast<std::size_t>(Info.NBins) + static_cast<std::size_t>(Info.NPars);
    Outputs.Allocate(NOut, Bytes);
    MC = Outputs.get();
    W2 = Outputs.get() + Info.NBins;
    Grad = Outputs.get() + 2 * Info.NBins;

    // Everything runs on one stream, staged through pinned host memory, with one wait per reweight
    DR_CHECK(cudaStreamCreate(&Stream));
    HostOut.Reserve(NOut);
    HostPars.Reserve(static_cast<std::size_t>(Info.NPars));
    HostGathered.Reserve(static_cast<std::size_t>(Info.NGathered));
    HostSegments.Reserve(static_cast<std::size_t>(Info.NSplineParams));

    int Device = 0;
    int NMultiprocessors = 0;
    DR_CHECK(cudaGetDevice(&Device));
    DR_CHECK(cudaDeviceGetAttribute(&NMultiprocessors, cudaDevAttrMultiProcessorCount, Device));
    GridStrideBlocks = 4 * NMultiprocessors;

    FullGradShared = static_cast<std::size_t>(Info.NPars) * sizeof(double);
    if (FullGradShared > 48 * 1024) {
      DR_CHECK(cudaFuncSetAttribute(FullPassGradient, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(FullGradShared)));
    }
  }

  ~DeviceReweightImpl() {
    if (Stream != nullptr) cudaStreamDestroy(Stream);
  }

  /// Upload what full passes need, see DeviceSweepStatic
  void SetSweep(const DeviceSweepStatic& Sweep) {
    if (Sweep.NVariables > kMaxVariables) throw std::runtime_error("DeviceReweight: too many shifted variables per event");
    if (Sweep.EdgeStart.size() < 2) throw std::runtime_error("DeviceReweight: no binning for the full pass");

    NLive = Sweep.NLive;
    NVars = Sweep.NVariables;
    NShifts = Sweep.NShifts;
    UnderOverFlow = Sweep.UnderOverFlowBin;
    NotSelected = Sweep.NotSelected;

    Full.Sample.Upload(Sweep.Sample, Bytes);
    Full.Sig.Upload(Sweep.Sig, Bytes);
    Full.ConstWeight.Upload(Sweep.ConstWeight, Bytes);
    Full.Nominal.Upload(Sweep.Nominal, Bytes);
    Full.CoefStart.Upload(Sweep.CoefStart, Bytes);
    Full.CoefShift.Upload(Sweep.CoefShift, Bytes);
    Full.Coef.Upload(Sweep.Coef, Bytes);
    Full.CutStart.Upload(Sweep.CutStart, Bytes);
    Full.CutVar.Upload(Sweep.CutVariable, Bytes);
    Full.CutLo.Upload(Sweep.CutLo, Bytes);
    Full.CutHi.Upload(Sweep.CutHi, Bytes);
    Full.EdgeStart.Upload(Sweep.EdgeStart, Bytes);
    Full.Edges.Upload(Sweep.Edges, Bytes);
    Full.Stride.Upload(Sweep.Stride, Bytes);
    Full.Offset.Upload(Sweep.Offset, Bytes);
    Full.ShiftValue.Allocate(static_cast<std::size_t>(NShifts), Bytes);
    Full.EventBin.Allocate(static_cast<std::size_t>(NLive), Bytes);
    Full.EventWeight.Allocate(static_cast<std::size_t>(NLive), Bytes);
    HostShifts.Reserve(static_cast<std::size_t>(NShifts));

    FullHistShared = static_cast<std::size_t>(2 * Info.NBins) * sizeof(double);
    if (FullHistShared > 48 * 1024) {
      DR_CHECK(cudaFuncSetAttribute(FullPass, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(FullHistShared)));
    }
    HaveSweep = true;
  }

  /// Replace the data
  void SetData(const std::vector<double>& NewData) {
    if (static_cast<int>(NewData.size()) != Info.NBins) throw std::runtime_error("DeviceReweight: data has the wrong length");
    DR_CHECK(cudaMemcpy(Data.get(), NewData.data(), NewData.size() * sizeof(double), cudaMemcpyHostToDevice));
  }

  /// Use groups built on the host (after a CPU full pass), see DeviceReweightGroups
  void SetGroups(const DeviceReweightGroups& FromHost, const std::vector<double>& NewData) {
    if (static_cast<int>(FromHost.BinStart.size()) != Info.NBins + 1) throw std::runtime_error("DeviceReweight: BinStart has the wrong length");
    SetData(NewData);

    // The bin of every group, which the gradient needs, from the bins' group ranges
    std::vector<int> GroupBin(FromHost.Sig.size());
    for (int b = 0; b < Info.NBins; ++b) {
      for (int g = FromHost.BinStart[b]; g < FromHost.BinStart[b + 1]; ++g) GroupBin[static_cast<std::size_t>(g)] = b;
    }

    HostGroups.Sig.Upload(FromHost.Sig, Bytes);
    HostGroups.Const.Upload(FromHost.Const, Bytes);
    HostGroups.Const2.Upload(FromHost.Const2, Bytes);
    HostGroups.BinStart.Upload(FromHost.BinStart, Bytes);
    HostGroups.Bin.Upload(GroupBin, Bytes);
    HostGroups.Weight.Allocate(FromHost.Sig.size(), Bytes);
    HostGroups.Count = static_cast<int>(FromHost.Sig.size());
    Current = &HostGroups;
  }

  // ----------------------------------------------------------------------------------------------
  // Reweights
  // ----------------------------------------------------------------------------------------------

  /// Leapfrog reweight: splines, group weights, bin sums, gradient. See the kernels' steps 1 to 5
  void Run(const double* ParValues, const short int* Segments, const double* Gathered, const bool FillW2, const bool Gradient,
           double* MCOut, double* W2Out, double* GradOut) {
    if (Current == nullptr) throw std::runtime_error("DeviceReweight: Run before any groups");
    UploadInputs(ParValues, Segments, Gathered);
    EnsureLeapfrogScratch();
    const Groups& G = *Current;

    // Step 1: group weights
    if (G.Count > 0) {
      WeighGroups<<<BlocksFor(G.Count), kThreads, 0, Stream>>>(G.Count, G.Sig.get(), G.Const.get(), G.Const2.get(),
                                                               Sig.WeightStart.get(), Sig.WeightIdx.get(), Pool.get(),
                                                               G.Weight.get(), GroupW2.get(), FillW2);
      DR_CHECK(cudaGetLastError());
    }

    // Step 2: each bin's total, a segmented sum over its range of groups
    std::size_t TempSize = SegmentedSumTempBytes;
    DR_CHECK(cub::DeviceSegmentedReduce::Sum(SegmentedSumTemp.get(), TempSize, G.Weight.get(), MC, Info.NBins,
                                             G.BinStart.get(), G.BinStart.get() + 1, Stream));
    if (FillW2) {
      TempSize = SegmentedSumTempBytes;
      DR_CHECK(cub::DeviceSegmentedReduce::Sum(SegmentedSumTemp.get(), TempSize, GroupW2.get(), W2, Info.NBins,
                                               G.BinStart.get(), G.BinStart.get() + 1, Stream));
    }

    // Steps 3 to 5: the gradient
    if (Gradient) {
      DR_CHECK(cudaMemsetAsync(SigTotal.get(), 0, SigTotal.size() * sizeof(double), Stream));
      if (G.Count > 0) {
        SumBySignature<<<BlocksFor(G.Count), kThreads, 0, Stream>>>(G.Count, G.Bin.get(), G.Sig.get(), G.Weight.get(), MC,
                                                                    Data.get(), Info.LowMCBound, SigTotal.get());
        DR_CHECK(cudaGetLastError());
      }
      if (Entries.NChunks > 0) {
        SumChunks<<<Entries.NChunks, kThreads, 0, Stream>>>(Entries.ChunkBegin.get(), Entries.Sig.get(), Entries.DerivIdx.get(),
                                                            Entries.WeightIdx.get(), SigTotal.get(), Pool.get(),
                                                            Entries.ChunkTotal.get());
        DR_CHECK(cudaGetLastError());
      }
      SumParameters<<<BlocksFor(Info.NPars), kThreads, 0, Stream>>>(Info.NPars, Entries.ParChunkStart.get(),
                                                                    Entries.ChunkTotal.get(), Grad);
      DR_CHECK(cudaGetLastError());
    }

    CopyOutputsBack(FillW2, Gradient, MCOut, W2Out, GradOut);
  }

  /// Full pass: every live event shifted, cut, binned and weighed. Its bins stay on the device
  void RunFull(const double* ParValues, const short int* Segments, const double* Gathered, const double* ShiftValues,
               const bool FillW2, const bool Gradient, double* MCOut, double* W2Out, double* GradOut) {
    if (!HaveSweep) throw std::runtime_error("DeviceReweight: RunFull before SetSweep");
    UploadInputs(ParValues, Segments, Gathered);
    if (NShifts > 0) {
      std::memcpy(HostShifts.get(), ShiftValues, static_cast<std::size_t>(NShifts) * sizeof(double));
      DR_CHECK(cudaMemcpyAsync(Full.ShiftValue.get(), HostShifts.get(), static_cast<std::size_t>(NShifts) * sizeof(double),
                               cudaMemcpyHostToDevice, Stream));
    }

    // The kernel adds into the totals, so they start at zero
    DR_CHECK(cudaMemsetAsync(MC, 0, static_cast<std::size_t>(Info.NBins) * sizeof(double), Stream));
    if (FillW2) DR_CHECK(cudaMemsetAsync(W2, 0, static_cast<std::size_t>(Info.NBins) * sizeof(double), Stream));

    FullPass<<<GridStrideBlocks, kThreads, FullHistShared, Stream>>>(
        NLive, NVars, Info.NBins, Full.Sample.get(), Full.Sig.get(), Full.ConstWeight.get(), Full.Nominal.get(),
        Full.CoefStart.get(), Full.CoefShift.get(), Full.Coef.get(), Full.ShiftValue.get(), Full.CutStart.get(),
        Full.CutVar.get(), Full.CutLo.get(), Full.CutHi.get(), Full.EdgeStart.get(), Full.Edges.get(), Full.Stride.get(),
        Full.Offset.get(), Sig.WeightStart.get(), Sig.WeightIdx.get(), Pool.get(), UnderOverFlow, NotSelected, FillW2,
        Full.EventBin.get(), Full.EventWeight.get(), MC, W2);
    DR_CHECK(cudaGetLastError());

    if (Gradient) {
      DR_CHECK(cudaMemsetAsync(Grad, 0, static_cast<std::size_t>(Info.NPars) * sizeof(double), Stream));
      FullPassGradient<<<GridStrideBlocks, kThreads, FullGradShared, Stream>>>(
          NLive, Info.NPars, Full.EventBin.get(), Full.Sig.get(), Full.EventWeight.get(), Sig.DerivStart.get(),
          Sig.DerivPar.get(), Sig.DerivIdx.get(), Sig.DerivWeightIdx.get(), Pool.get(), MC, Data.get(), Info.LowMCBound, Grad);
      DR_CHECK(cudaGetLastError());
    }

    CopyOutputsBack(FillW2, Gradient, MCOut, W2Out, GradOut);
  }

  /// The bin (or code) of every live event from the last RunFull
  void GetEventBins(int* Out) const {
    if (!HaveSweep) throw std::runtime_error("DeviceReweight: GetEventBins before SetSweep");
    DR_CHECK(cudaStreamSynchronize(Stream));
    DR_CHECK(cudaMemcpy(Out, Full.EventBin.get(), static_cast<std::size_t>(NLive) * sizeof(int), cudaMemcpyDeviceToHost));
  }

  /// Groups for the leapfrogs, built here from the last RunFull's event bins. See steps 1 to 6 of the
  /// group building kernels
  int BuildGroupsFromFullPass() {
    if (!HaveSweep) throw std::runtime_error("DeviceReweight: BuildGroupsFromFullPass before SetSweep");
    EnsureGroupBuildingScratch();
    GroupScratch& S = Scratch;

    // Step 1: keys
    DR_CHECK(cudaMemsetAsync(S.Counts.get(), 0, 2 * sizeof(int), Stream));
    MakeGroupKeys<<<BlocksFor(NLive), kThreads, 0, Stream>>>(NLive, NSig, UnderOverFlow, Full.EventBin.get(), Full.Sig.get(),
                                                             Full.ConstWeight.get(), S.NoGroupKey, S.Keys.get(),
                                                             S.Index.get(), S.Counts.get());
    DR_CHECK(cudaGetLastError());

    // Step 2: sort
    std::size_t TempSize = S.TempBytes;
    DR_CHECK(cub::DeviceRadixSort::SortPairs(S.Temp.get(), TempSize, S.Keys.get(), S.SortedKeys.get(), S.Index.get(),
                                             S.SortedIndex.get(), NLive, 0, S.KeyBits, Stream));
    const int NWithGroup = ReadCount(0);

    int NGroups = 0;
    if (NWithGroup > 0) {
      // Step 3: constant weights in sorted order
      GatherConstWeights<<<BlocksFor(NWithGroup), kThreads, 0, Stream>>>(NWithGroup, S.SortedIndex.get(),
                                                                          Full.ConstWeight.get(), S.C.get(), S.C2.get());
      DR_CHECK(cudaGetLastError());

      // Step 4: one group per run of equal keys, summing C and C^2
      TempSize = S.TempBytes;
      DR_CHECK(cub::DeviceReduce::ReduceByKey(S.Temp.get(), TempSize, S.SortedKeys.get(), S.GroupKeys.get(), S.C.get(),
                                              DeviceGroups.Const.get(), S.Counts.get() + 1, cub::Sum(), NWithGroup, Stream));
      TempSize = S.TempBytes;
      DR_CHECK(cub::DeviceReduce::ReduceByKey(S.Temp.get(), TempSize, S.SortedKeys.get(), S.GroupKeys.get(), S.C2.get(),
                                              DeviceGroups.Const2.get(), S.Counts.get() + 1, cub::Sum(), NWithGroup, Stream));
      NGroups = ReadCount(1);

      // Step 5: signature and bin of each group
      DecodeGroupKeys<<<BlocksFor(NGroups), kThreads, 0, Stream>>>(NGroups, NSig, S.GroupKeys.get(), DeviceGroups.Sig.get(),
                                                                   DeviceGroups.Bin.get());
      DR_CHECK(cudaGetLastError());
    }

    // Step 6: where each bin's groups start
    FindBinStarts<<<BlocksFor(Info.NBins + 1), kThreads, 0, Stream>>>(Info.NBins, NGroups, DeviceGroups.Bin.get(),
                                                                      DeviceGroups.BinStart.get());
    DR_CHECK(cudaGetLastError());

    DeviceGroups.Count = NGroups;
    Current = &DeviceGroups;
    return NGroups;
  }

  /// Device memory allocated so far, in bytes
  std::size_t Bytes = 0;

 private:
  // ----------------------------------------------------------------------------------------------
  // Data
  // ----------------------------------------------------------------------------------------------

  /// One set of groups, see DeviceReweightGroups. Count of them, sorted by bin
  struct Groups {
    DeviceArray<int> Sig;       // signature of each group
    DeviceArray<double> Const;  // summed constant weights of each group
    DeviceArray<double> Const2; // summed squared constant weights
    DeviceArray<int> BinStart;  // bin b's groups are BinStart[b] to BinStart[b+1]-1
    DeviceArray<int> Bin;       // bin of each group
    DeviceArray<double> Weight; // each group's weight in the last reweight
    int Count = 0;
  };

  /// Signatures: their varying weights and derivatives, as pool indices, see DeviceReweightStatic
  struct Signatures {
    DeviceArray<int> WeightStart, WeightIdx;
    DeviceArray<int> DerivStart, DerivPar, DerivIdx, DerivWeightIdx;
  };

  /// The derivative entries regrouped by parameter and cut into chunks, for steps 4 and 5
  struct GradientEntries {
    DeviceArray<int> Sig;           // signature of each entry
    DeviceArray<int> DerivIdx;      // pool index of the entry's derivative
    DeviceArray<int> WeightIdx;     // pool index of the weight it differentiates
    DeviceArray<int> ChunkBegin;    // chunk c is entries ChunkBegin[c] to ChunkBegin[c+1]-1
    DeviceArray<int> ParChunkStart; // parameter p's chunks are ParChunkStart[p] to ParChunkStart[p+1]-1
    DeviceArray<double> ChunkTotal; // each chunk's sum
    int NChunks = 0;
  };

  /// Splines, as BinnedSplineHandler lays them out
  struct Splines {
    DeviceArray<double> KnotX, Coeffs;
    DeviceArray<unsigned int> CoeffStart;
    DeviceArray<short int> SplineParam, Segment;
    DeviceArray<int> GlobalParam;
  };

  /// What a full pass needs, see DeviceSweepStatic
  struct FullPassData {
    DeviceArray<int> Sample, Sig;
    DeviceArray<double> ConstWeight;
    DeviceArray<double> Nominal, Coef;
    DeviceArray<int> CoefStart, CoefShift;
    DeviceArray<int> CutStart, CutVar;
    DeviceArray<double> CutLo, CutHi;
    DeviceArray<int> EdgeStart, Stride, Offset;
    DeviceArray<double> Edges;
    DeviceArray<double> ShiftValue;
    DeviceArray<int> EventBin;        // output: each live event's bin
    DeviceArray<double> EventWeight;  // output: each live event's weight
  };

  /// Working space for building groups on the device
  struct GroupScratch {
    DeviceArray<unsigned long long> Keys, SortedKeys, GroupKeys;
    DeviceArray<int> Index, SortedIndex, Counts;
    DeviceArray<double> C, C2;
    DeviceArray<char> Temp;
    std::size_t TempBytes = 0;
    unsigned long long NoGroupKey = 0;
    int KeyBits = 64;
  };

  DeviceReweightStatic Info;  // sizes and settings (the big vectors are released after upload)
  int NSig = 0;               // number of signatures
  int NPool = 0;              // pool size
  int GatheredOffset = 0;        // where the host gathered values start in the pool
  int GridStrideBlocks = 0;   // blocks for the grid stride kernels

  DeviceArray<double> Pool;
  Splines Spl;
  Signatures Sig;
  GradientEntries Entries;
  DeviceArray<double> SigTotal;  // step 3's per signature totals
  DeviceArray<double> Data;

  DeviceArray<double> Outputs;   // [MC | W2 | gradient]
  double* MC = nullptr;
  double* W2 = nullptr;
  double* Grad = nullptr;
  std::size_t NOut = 0;

  Groups HostGroups;             // from SetGroups
  Groups DeviceGroups;           // from BuildGroupsFromFullPass
  Groups* Current = nullptr;     // the ones the next Run uses
  DeviceArray<double> GroupW2;   // per group W2, scratch
  DeviceArray<char> SegmentedSumTemp;
  std::size_t SegmentedSumTempBytes = 0;

  bool HaveSweep = false;
  int NLive = 0, NVars = 0, NShifts = 0, UnderOverFlow = -1, NotSelected = -2;
  FullPassData Full;
  GroupScratch Scratch;
  std::size_t FullHistShared = 0;  // shared memory of the full pass kernel
  std::size_t FullGradShared = 0;  // shared memory of the full pass gradient kernel

  cudaStream_t Stream = nullptr;
  PinnedArray<double> HostOut, HostPars, HostGathered, HostShifts;
  PinnedArray<short int> HostSegments;

  // ----------------------------------------------------------------------------------------------
  // Setup helpers
  // ----------------------------------------------------------------------------------------------

  /// The pool: parameters, spline weights, spline derivatives, gathered values, then the constants 1 and 0
  void SetUpPool() {
    NPool = Info.NPars + 2 * Info.NSplines + Info.NGathered + 2;
    GatheredOffset = Info.NPars + 2 * Info.NSplines;
    Pool.Allocate(static_cast<std::size_t>(NPool), Bytes);
    const double Constants[2] = {1., 0.};
    DR_CHECK(cudaMemcpy(Pool.get() + NPool - 2, Constants, 2 * sizeof(double), cudaMemcpyHostToDevice));
  }

  void UploadSplines() {
    Spl.KnotX.Upload(Info.KnotX, Bytes);
    Spl.Coeffs.Upload(Info.SegmentCoeffs, Bytes);
    Spl.CoeffStart.Upload(Info.SplineCoeffStart, Bytes);
    Spl.SplineParam.Upload(Info.SplineParam, Bytes);
    Spl.GlobalParam.Upload(Info.SplineParamIndex, Bytes);
    Spl.Segment.Allocate(static_cast<std::size_t>(Info.NSplineParams), Bytes);
    // Only the sizes are needed from here on
    Info.KnotX.clear();
    Info.SegmentCoeffs.clear();
    Info.SplineCoeffStart.clear();
    Info.SplineParam.clear();
    Info.SplineParamIndex.clear();
  }

  void UploadSignatures() {
    Sig.WeightStart.Upload(Info.SigWeightStart, Bytes);
    Sig.WeightIdx.Upload(Info.SigWeightIdx, Bytes);
    Sig.DerivStart.Upload(Info.SigDerivStart, Bytes);
    Sig.DerivPar.Upload(Info.SigDerivPar, Bytes);
    Sig.DerivIdx.Upload(Info.SigDerivIdx, Bytes);
    Sig.DerivWeightIdx.Upload(Info.SigDerivWeightIdx, Bytes);
    SigTotal.Allocate(static_cast<std::size_t>(NSig), Bytes);
  }

  /// Regroup the signatures' derivative entries by parameter (a counting sort), then cut each
  /// parameter's entries into chunks of at most kChunk
  void BuildGradientEntries() {
    const std::size_t NPars = static_cast<std::size_t>(Info.NPars);
    const std::size_t NEntries = Info.SigDerivPar.size();

    // Where each parameter's entries start
    std::vector<int> ParStart(NPars + 1, 0);
    for (const int Par : Info.SigDerivPar) ++ParStart[static_cast<std::size_t>(Par) + 1];
    for (std::size_t p = 0; p < NPars; ++p) ParStart[p + 1] += ParStart[p];

    // Place every entry under its parameter, keeping signature order within a parameter
    std::vector<int> EntrySig(NEntries), EntryDeriv(NEntries), EntryWeight(NEntries);
    std::vector<int> Next(ParStart.begin(), ParStart.end() - 1);
    for (int s = 0; s < NSig; ++s) {
      for (int k = Info.SigDerivStart[static_cast<std::size_t>(s)]; k < Info.SigDerivStart[static_cast<std::size_t>(s) + 1]; ++k) {
        const std::size_t Entry = static_cast<std::size_t>(k);
        const std::size_t At = static_cast<std::size_t>(Next[static_cast<std::size_t>(Info.SigDerivPar[Entry])]++);
        EntrySig[At] = s;
        EntryDeriv[At] = Info.SigDerivIdx[Entry];
        EntryWeight[At] = Info.SigDerivWeightIdx[Entry];
      }
    }

    // Chunks: a parameter's entries in pieces of at most kChunk. Chunks are consecutive in entry
    // order, so a chunk ends where the next one starts, and the last one ends at NEntries
    std::vector<int> ChunkBegin;
    std::vector<int> ParChunkStart(NPars + 1, 0);
    for (std::size_t p = 0; p < NPars; ++p) {
      for (int k = ParStart[p]; k < ParStart[p + 1]; k += kChunk) ChunkBegin.push_back(k);
      ParChunkStart[p + 1] = static_cast<int>(ChunkBegin.size());
    }
    Entries.NChunks = static_cast<int>(ChunkBegin.size());
    ChunkBegin.push_back(static_cast<int>(NEntries));

    Entries.Sig.Upload(EntrySig, Bytes);
    Entries.DerivIdx.Upload(EntryDeriv, Bytes);
    Entries.WeightIdx.Upload(EntryWeight, Bytes);
    Entries.ChunkBegin.Upload(ChunkBegin, Bytes);
    Entries.ParChunkStart.Upload(ParChunkStart, Bytes);
    Entries.ChunkTotal.Allocate(static_cast<std::size_t>(Entries.NChunks), Bytes);
  }

  /// Per group W2 scratch and the segmented sum's working space, allocated once at the most groups
  /// there can be (one per live event, or the host groups if there are more)
  void EnsureLeapfrogScratch() {
    const std::size_t MostGroups = static_cast<std::size_t>(std::max({1, NLive, Current->Count}));
    if (GroupW2.size() < MostGroups) GroupW2.Allocate(MostGroups, Bytes);
    if (SegmentedSumTemp.empty()) {
      DR_CHECK(cub::DeviceSegmentedReduce::Sum(nullptr, SegmentedSumTempBytes, Current->Weight.get(), MC, Info.NBins,
                                               Current->BinStart.get(), Current->BinStart.get() + 1, Stream));
      SegmentedSumTemp.Allocate(SegmentedSumTempBytes, Bytes);
    }
  }

  /// Arrays for building groups on the device, allocated on first use at one group per live event
  void EnsureGroupBuildingScratch() {
    if (!Scratch.Keys.empty()) return;
    const std::size_t N = static_cast<std::size_t>(std::max(1, NLive));
    Scratch.Keys.Allocate(N, Bytes);
    Scratch.SortedKeys.Allocate(N, Bytes);
    Scratch.GroupKeys.Allocate(N, Bytes);
    Scratch.Index.Allocate(N, Bytes);
    Scratch.SortedIndex.Allocate(N, Bytes);
    Scratch.C.Allocate(N, Bytes);
    Scratch.C2.Allocate(N, Bytes);
    Scratch.Counts.Allocate(2, Bytes);
    DeviceGroups.Sig.Allocate(N, Bytes);
    DeviceGroups.Const.Allocate(N, Bytes);
    DeviceGroups.Const2.Allocate(N, Bytes);
    DeviceGroups.Bin.Allocate(N, Bytes);
    DeviceGroups.Weight.Allocate(N, Bytes);
    DeviceGroups.BinStart.Allocate(static_cast<std::size_t>(Info.NBins) + 1, Bytes);

    // Keys need only enough bits for the largest real key plus the "no group" marker above it,
    // which keeps the radix sort to as few passes as possible
    Scratch.NoGroupKey = (static_cast<unsigned long long>(Info.NBins) * static_cast<unsigned long long>(NSig)) << 1;
    Scratch.KeyBits = 1;
    while ((1ull << Scratch.KeyBits) <= Scratch.NoGroupKey) ++Scratch.KeyBits;

    // One working space big enough for both the sort and the reduce by key
    std::size_t SortBytes = 0;
    std::size_t ReduceBytes = 0;
    DR_CHECK(cub::DeviceRadixSort::SortPairs(nullptr, SortBytes, Scratch.Keys.get(), Scratch.SortedKeys.get(),
                                             Scratch.Index.get(), Scratch.SortedIndex.get(), NLive, 0, Scratch.KeyBits));
    DR_CHECK(cub::DeviceReduce::ReduceByKey(nullptr, ReduceBytes, Scratch.SortedKeys.get(), Scratch.GroupKeys.get(),
                                            Scratch.C.get(), DeviceGroups.Const.get(), Scratch.Counts.get(), cub::Sum(), NLive));
    Scratch.TempBytes = std::max(SortBytes, ReduceBytes);
    Scratch.Temp.Allocate(Scratch.TempBytes, Bytes);
  }

  // ----------------------------------------------------------------------------------------------
  // Per reweight helpers
  // ----------------------------------------------------------------------------------------------

  /// Parameter values and gathered values into the pool, then the splines evaluated at them
  void UploadInputs(const double* ParValues, const short int* Segments, const double* Gathered) {
    const std::size_t NPars = static_cast<std::size_t>(Info.NPars);
    std::memcpy(HostPars.get(), ParValues, NPars * sizeof(double));
    DR_CHECK(cudaMemcpyAsync(Pool.get(), HostPars.get(), NPars * sizeof(double), cudaMemcpyHostToDevice, Stream));

    if (Info.NGathered > 0) {
      const std::size_t NGathered = static_cast<std::size_t>(Info.NGathered);
      std::memcpy(HostGathered.get(), Gathered, NGathered * sizeof(double));
      DR_CHECK(cudaMemcpyAsync(Pool.get() + GatheredOffset, HostGathered.get(), NGathered * sizeof(double), cudaMemcpyHostToDevice, Stream));
    }

    if (Info.NSplines > 0) {
      const std::size_t NSegments = static_cast<std::size_t>(Info.NSplineParams);
      std::memcpy(HostSegments.get(), Segments, NSegments * sizeof(short int));
      DR_CHECK(cudaMemcpyAsync(Spl.Segment.get(), HostSegments.get(), NSegments * sizeof(short int), cudaMemcpyHostToDevice, Stream));
      EvaluateSplines<<<BlocksFor(Info.NSplines), kThreads, 0, Stream>>>(
          Info.NSplines, Info.MaxKnots, Spl.SplineParam.get(), Spl.CoeffStart.get(), Spl.Coeffs.get(), Spl.KnotX.get(),
          Spl.Segment.get(), Spl.GlobalParam.get(), Pool.get(), Info.NPars, Info.NPars + Info.NSplines);
      DR_CHECK(cudaGetLastError());
    }
  }

  /// One copy of whatever this reweight produced ([MC], [MC | W2] or [MC | W2 | gradient]), the
  /// reweight's only wait, then out of the pinned buffer into the caller's arrays
  void CopyOutputsBack(const bool FillW2, const bool Gradient, double* MCOut, double* W2Out, double* GradOut) {
    const std::size_t NBins = static_cast<std::size_t>(Info.NBins);
    const std::size_t NCopy = Gradient ? NOut : (FillW2 ? 2 * NBins : NBins);
    DR_CHECK(cudaMemcpyAsync(HostOut.get(), Outputs.get(), NCopy * sizeof(double), cudaMemcpyDeviceToHost, Stream));
    DR_CHECK(cudaStreamSynchronize(Stream));

    std::memcpy(MCOut, HostOut.get(), NBins * sizeof(double));
    if (FillW2) std::memcpy(W2Out, HostOut.get() + NBins, NBins * sizeof(double));
    if (Gradient) std::memcpy(GradOut, HostOut.get() + 2 * NBins, static_cast<std::size_t>(Info.NPars) * sizeof(double));
  }

  /// One of the group building counters, after waiting for the stream
  int ReadCount(const int Which) {
    int Value = 0;
    DR_CHECK(cudaMemcpyAsync(&Value, Scratch.Counts.get() + Which, sizeof(int), cudaMemcpyDeviceToHost, Stream));
    DR_CHECK(cudaStreamSynchronize(Stream));
    return Value;
  }
};

// ================================================================================================
// 4. DeviceReweight: forwards to the implementation
// ================================================================================================

DeviceReweight::DeviceReweight(const DeviceReweightStatic& Static) : Impl(new DeviceReweightImpl(Static)) {}
DeviceReweight::~DeviceReweight() = default;

void DeviceReweight::SetGroups(const DeviceReweightGroups& FromHost, const std::vector<double>& Data) {
  Impl->SetGroups(FromHost, Data);
}

void DeviceReweight::SetData(const std::vector<double>& NewData) {
  Impl->SetData(NewData);
}

void DeviceReweight::Run(const double* ParValues, const short int* Segments, const double* Gathered, const bool FillW2,
                         const bool Gradient, double* MCOut, double* W2Out, double* GradOut) {
  Impl->Run(ParValues, Segments, Gathered, FillW2, Gradient, MCOut, W2Out, GradOut);
}

std::size_t DeviceReweight::DeviceBytes() const {
  return Impl->Bytes;
}

void DeviceReweight::SetSweep(const DeviceSweepStatic& Sweep) {
  Impl->SetSweep(Sweep);
}

void DeviceReweight::RunFull(const double* ParValues, const short int* Segments, const double* Gathered, const double* ShiftValues,
                             const bool FillW2, const bool Gradient, double* MCOut, double* W2Out, double* GradOut) {
  Impl->RunFull(ParValues, Segments, Gathered, ShiftValues, FillW2, Gradient, MCOut, W2Out, GradOut);
}

void DeviceReweight::GetEventBins(int* Out) const {
  Impl->GetEventBins(Out);
}

int DeviceReweight::BuildGroupsFromFullPass() {
  return Impl->BuildGroupsFromFullPass();
}
