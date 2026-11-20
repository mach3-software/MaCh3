// OscillationGradientDevice: OscillationGradient's derivatives on an NVIDIA GPU, see OscillationGradientDevice.h
#include "Samples/Oscillation/OscillationGradientDevice.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace {

void Check(const cudaError_t Status, const char* What) {
  if (Status != cudaSuccess) throw std::runtime_error(std::string("CUDA ") + What + ": " + cudaGetErrorString(Status));
}

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
  template <typename V>
  void Upload(const std::vector<V>& Values) {
    static_assert(sizeof(V) == sizeof(T), "element size");
    Allocate(Values.size());
    if (!Values.empty()) Check(cudaMemcpy(Ptr, Values.data(), Values.size() * sizeof(T), cudaMemcpyHostToDevice), "upload");
  }
};

/// Small blocks: a thread's dual numbers take a lot of registers
constexpr int kBlock = 128;

/// Thread (w, blockIdx.y): derivative of weight w with respect to tangent slot blockIdx.y
__global__ void __launch_bounds__(kBlock) WeightDerivatives(const int NWeights, const int NFree,
                                                           const M3::Osc::Constants<M3::Osc::OscDual1>* __restrict__ C,
                                                           const double* __restrict__ Energy, const char* __restrict__ Anti,
                                                           const signed char* __restrict__ From,
                                                           const signed char* __restrict__ To, double* __restrict__ Out) {
  const int w = blockIdx.x * blockDim.x + threadIdx.x;
  if (w >= NWeights) return;
  const int j = blockIdx.y;
  const M3::Osc::OscDual1 P = M3::Osc::Probability(C[j], Energy[w], Anti[w] != 0, From[w], To[w]);
  Out[static_cast<std::size_t>(w) * NFree + j] = P.d[0];
}

}  // namespace

struct OscillationGradientDevice::Impl {
  int N = 0;
  DeviceArray<double> Energy;
  DeviceArray<char> Anti;
  DeviceArray<signed char> From;
  DeviceArray<signed char> To;
  DeviceArray<M3::Osc::Constants<M3::Osc::OscDual1>> Tangents;
  DeviceArray<double> Out;
  /// Pinned staging for the sliced constants and the results, so both copies are fast
  M3::Osc::Constants<M3::Osc::OscDual1>* HostTangents = nullptr;
  double* HostOut = nullptr;
  std::size_t HostOutSize = 0;
  cudaStream_t Stream = nullptr;
  ~Impl() {
    if (HostTangents) cudaFreeHost(HostTangents);
    if (HostOut) cudaFreeHost(HostOut);
    if (Stream) cudaStreamDestroy(Stream);
  }
};

// ************************************************
OscillationGradientDevice::OscillationGradientDevice(const std::vector<double>& Energy, const std::vector<char>& Anti,
                                                     const std::vector<signed char>& From, const std::vector<signed char>& To)
    : D(std::make_unique<Impl>()) {
// ************************************************
  D->N = static_cast<int>(Energy.size());
  D->Energy.Upload(Energy);
  D->Anti.Upload(Anti);
  D->From.Upload(From);
  D->To.Upload(To);
  D->Tangents.Allocate(M3::Osc::kNPars);
  Check(cudaMallocHost(&D->HostTangents, M3::Osc::kNPars * sizeof(M3::Osc::Constants<M3::Osc::OscDual1>)), "host malloc");
  Check(cudaStreamCreateWithFlags(&D->Stream, cudaStreamNonBlocking), "stream");
}

// ************************************************
OscillationGradientDevice::~OscillationGradientDevice() = default;
// ************************************************

// ************************************************
bool OscillationGradientDevice::Available() {
// ************************************************
  int Count = 0;
  return cudaGetDeviceCount(&Count) == cudaSuccess && Count > 0;
}

// ************************************************
void OscillationGradientDevice::Evaluate(const M3::Osc::Constants<M3::Osc::OscDual>& C, const std::vector<int>& FreeSlot,
                                         double* Out) {
// ************************************************
  const int NFree = static_cast<int>(FreeSlot.size());
  if (NFree == 0 || D->N == 0) return;
  const std::size_t NOut = static_cast<std::size_t>(D->N) * static_cast<std::size_t>(NFree);
  if (D->Out.N != NOut) {
    D->Out.Allocate(NOut);
    if (D->HostOut) cudaFreeHost(D->HostOut);
    Check(cudaMallocHost(&D->HostOut, NOut * sizeof(double)), "host malloc");
    D->HostOutSize = NOut;
  }
  for (int j = 0; j < NFree; ++j) D->HostTangents[j] = M3::Osc::SliceTangent(C, FreeSlot[static_cast<std::size_t>(j)]);
  Check(cudaMemcpyAsync(D->Tangents.Ptr, D->HostTangents, static_cast<std::size_t>(NFree) * sizeof(M3::Osc::Constants<M3::Osc::OscDual1>),
                        cudaMemcpyHostToDevice, D->Stream), "copy constants");
  const dim3 Grid(static_cast<unsigned int>((D->N + kBlock - 1) / kBlock), static_cast<unsigned int>(NFree));
  WeightDerivatives<<<Grid, kBlock, 0, D->Stream>>>(D->N, NFree, D->Tangents.Ptr, D->Energy.Ptr, D->Anti.Ptr, D->From.Ptr,
                                                    D->To.Ptr, D->Out.Ptr);
  Check(cudaGetLastError(), "launch");
  Check(cudaMemcpyAsync(D->HostOut, D->Out.Ptr, NOut * sizeof(double), cudaMemcpyDeviceToHost, D->Stream), "copy back");
  Check(cudaStreamSynchronize(D->Stream), "synchronize");
  std::copy(D->HostOut, D->HostOut + NOut, Out);
}
