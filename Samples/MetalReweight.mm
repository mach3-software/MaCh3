// MetalReweight: CpuReweight's cached pass and gradient on an Apple GPU, see MetalReweight.h
#include "Samples/MetalReweight.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <sstream>
#include <stdexcept>

namespace {

/// The kernels, compiled when the backend is made: Metal needs no offline compiler for that
const char* const kKernels = R"(
#include <metal_stdlib>
using namespace metal;

// Product of every signature's pool entries
kernel void sig_products(device const float* pool [[buffer(0)]], device const int* sig_start [[buffer(1)]],
                         device const int* sig_pool [[buffer(2)]], device float* sig_prod [[buffer(3)]],
                         constant uint& n_sig [[buffer(4)]], uint s [[thread_position_in_grid]]) {
  if (s >= n_sig) return;
  float p = 1.0f;
  for (int k = sig_start[s]; k < sig_start[s + 1]; ++k) p *= pool[sig_pool[k]];
  sig_prod[s] = p;
}

// Each bin's groups: the MC prediction and the sum of squared weights. One threadgroup per bin, its
// threads taking every 256th group, then a fixed order reduction: a bin holds thousands of groups,
// and one thread per bin left the GPU almost idle (1.2 ms per pass on an M5 Pro)
constant constexpr uint kBinWidth = 256;
kernel void bin_sums(device const int* bin_start [[buffer(0)]], device const int* group_sig [[buffer(1)]],
                     device const float* group_const [[buffer(2)]], device const float* group_const2 [[buffer(3)]],
                     device const float* sig_prod [[buffer(4)]], device float* mc [[buffer(5)]],
                     device float* w2 [[buffer(6)]], constant uint& n_bins [[buffer(7)]], uint b [[threadgroup_position_in_grid]],
                     uint t [[thread_index_in_threadgroup]]) {
  threadgroup float partial[kBinWidth];
  threadgroup float partial2[kBinWidth];
  float sum = 0.0f, sum2 = 0.0f;
  if (b < n_bins) {
    for (int g = bin_start[b] + int(t); g < bin_start[b + 1]; g += int(kBinWidth)) {
      const float p = sig_prod[group_sig[g]];
      const float w = group_const[g] * p;
      if (w <= 0.0f) continue;
      sum += w;
      sum2 += group_const2[g] * p * p;
    }
  }
  partial[t] = sum;
  partial2[t] = sum2;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  for (uint half_width = kBinWidth / 2; half_width > 0; half_width /= 2) {
    if (t < half_width) {
      partial[t] += partial[t + half_width];
      partial2[t] += partial2[t + half_width];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }
  if (t == 0 && b < n_bins) {
    mc[b] = partial[0];
    w2[b] = partial2[0];
  }
}

// Each signature's share of the bin factors, over its groups that contribute
kernel void sig_share(device const int* sig_group_start [[buffer(0)]], device const int* sig_group [[buffer(1)]],
                      device const int* group_bin [[buffer(2)]], device const float* group_const [[buffer(3)]],
                      device const float* sig_prod [[buffer(4)]], device const float* bin_factor [[buffer(5)]],
                      device float* share [[buffer(6)]], constant uint& n_sig [[buffer(7)]], uint s [[thread_position_in_grid]]) {
  if (s >= n_sig) return;
  float a = 0.0f;
  for (int k = sig_group_start[s]; k < sig_group_start[s + 1]; ++k) {
    const int g = sig_group[k];
    const float w = group_const[g] * sig_prod[s];
    if (w > 0.0f) a += bin_factor[group_bin[g]] * w;
  }
  share[s] = a;
}

// Each parameter's gradient: share times d(weight)/d(parameter) over the weight, over its terms.
// Each parameter's terms are cut into chunks of at most kChunk; one threadgroup per chunk, its threads
// taking every 256th term, then a fixed order reduction into the chunk's sum. A parameter every
// signature depends on has hundreds of thousands of terms, which one threadgroup took 0.27 ms over
constant constexpr uint kGroupWidth = 256;
kernel void param_grad_chunks(device const int* chunk_start [[buffer(0)]], device const int* chunk_end [[buffer(1)]],
                              device const int* term_sig [[buffer(2)]], device const int* term_deriv [[buffer(3)]],
                              device const int* term_pool [[buffer(4)]], device const float* deriv [[buffer(5)]],
                              device const float* pool [[buffer(6)]], device const float* share [[buffer(7)]],
                              device float* chunk_sum [[buffer(8)]], constant uint& n_chunks [[buffer(9)]],
                              uint c [[threadgroup_position_in_grid]], uint t [[thread_index_in_threadgroup]]) {
  threadgroup float partial[kGroupWidth];
  float g = 0.0f;
  if (c < n_chunks) {
    for (int k = chunk_start[c] + int(t); k < chunk_end[c]; k += int(kGroupWidth)) {
      const float d = deriv[term_deriv[k]];
      if (d != 0.0f) g += share[term_sig[k]] * d / pool[term_pool[k]];
    }
  }
  partial[t] = g;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  for (uint half_width = kGroupWidth / 2; half_width > 0; half_width /= 2) {
    if (t < half_width) partial[t] += partial[t + half_width];
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }
  if (t == 0 && c < n_chunks) chunk_sum[c] = partial[0];
}

// Each parameter's chunks, added in order, so the result is the same from run to run
kernel void param_grad_sum(device const int* par_chunk_start [[buffer(0)]], device const float* chunk_sum [[buffer(1)]],
                           device float* grad [[buffer(2)]], constant uint& n_par [[buffer(3)]], uint p [[thread_position_in_grid]]) {
  if (p >= n_par) return;
  float g = 0.0f;
  for (int c = par_chunk_start[p]; c < par_chunk_start[p + 1]; ++c) g += chunk_sum[c];
  grad[p] = g;
}
)";

/// @brief A shared buffer holding Count values of T
template <typename T>
id<MTLBuffer> MakeBuffer(id<MTLDevice> Device, const std::size_t Count) {
  return [Device newBufferWithLength:std::max<std::size_t>(1, Count) * sizeof(T) options:MTLResourceStorageModeShared];
}

/// @brief A shared buffer holding a copy of Values
template <typename T, typename U>
id<MTLBuffer> UploadBuffer(id<MTLDevice> Device, const std::vector<U>& Values) {
  id<MTLBuffer> Buffer = MakeBuffer<T>(Device, Values.size());
  T* Data = static_cast<T*>([Buffer contents]);
  for (std::size_t i = 0; i < Values.size(); ++i) Data[i] = static_cast<T>(Values[i]);
  return Buffer;
}

} // namespace

/// Everything Metal
struct MetalReweight::Device {
  id<MTLDevice> Gpu = nil;
  id<MTLCommandQueue> Queue = nil;
  id<MTLComputePipelineState> SigProducts = nil, BinSums = nil, SigShare = nil, ParamGradChunks = nil, ParamGradSum = nil;
  // The model, uploaded once
  id<MTLBuffer> Pool = nil, SigStart = nil, SigPool = nil, SigProd = nil;
  // The groups, uploaded whenever they are rebuilt
  unsigned int GroupsVersion = 0;
  id<MTLBuffer> BinStart = nil, GroupSig = nil, GroupConst = nil, GroupConst2 = nil, GroupBin = nil;
  id<MTLBuffer> SigGroupStart = nil, SigGroup = nil;
  // Outputs and gradient inputs
  id<MTLBuffer> MC = nil, W2 = nil, BinFactor = nil, Share = nil, Grad = nil;
  // The gradient terms, grouped by parameter, uploaded at AttachDerivatives
  id<MTLBuffer> ChunkStart = nil, ChunkEnd = nil, ParChunkStart = nil, ChunkSum = nil, TermSig = nil, TermDeriv = nil, TermPool = nil, Deriv = nil;
  std::size_t NDeriv = 0;
  unsigned int NChunks = 0;

  id<MTLComputePipelineState> Pipeline(id<MTLLibrary> Library, NSString* Name) {
    NSError* Error = nil;
    id<MTLFunction> Function = [Library newFunctionWithName:Name];
    id<MTLComputePipelineState> State = [Gpu newComputePipelineStateWithFunction:Function error:&Error];
    if (State == nil) throw std::runtime_error("Metal pipeline " + std::string([Name UTF8String]) + " failed");
    return State;
  }

  /// @brief One kernel launch: the kernel, its buffers, how many threads, and whether it runs one
  /// threadgroup of 256 per item rather than one thread
  struct Launch {
    id<MTLComputePipelineState> Kernel;
    std::vector<id<MTLBuffer>> Buffers;
    unsigned int N;
    bool GroupPerItem = false;
  };

  /// @brief Run kernels one after the other in one command buffer, and wait for the last
  void Run(const std::vector<Launch>& Launches) {
    id<MTLCommandBuffer> Commands = [Queue commandBuffer];
    id<MTLComputeCommandEncoder> Encoder = [Commands computeCommandEncoder];
    for (const Launch& L : Launches) {
      [Encoder setComputePipelineState:L.Kernel];
      NSUInteger Index = 0;
      for (id<MTLBuffer> Buffer : L.Buffers) [Encoder setBuffer:Buffer offset:0 atIndex:Index++];
      [Encoder setBytes:&L.N length:sizeof(L.N) atIndex:Index];
      if (L.GroupPerItem) {
        [Encoder dispatchThreadgroups:MTLSizeMake(std::max(1u, L.N), 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
      } else {
        const NSUInteger Width = std::min<NSUInteger>(L.Kernel.maxTotalThreadsPerThreadgroup, 256);
        [Encoder dispatchThreads:MTLSizeMake(std::max(1u, L.N), 1, 1) threadsPerThreadgroup:MTLSizeMake(Width, 1, 1)];
      }
      // Each kernel reads what the one before wrote
      [Encoder memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    [Encoder endEncoding];
    [Commands commit];
    [Commands waitUntilCompleted];
  }
};

// ************************************************
bool MetalReweight::Available() {
// ************************************************
  @autoreleasepool {
    return MTLCreateSystemDefaultDevice() != nil;
  }
}

// ************************************************
MetalReweight::MetalReweight(ReweightModel Model_, ReweightEventAccess& Access_)
    : CpuReweight(std::move(Model_), Access_), GPU(std::make_unique<Device>()) {
// ************************************************
  @autoreleasepool {
    GPU->Gpu = MTLCreateSystemDefaultDevice();
    if (GPU->Gpu == nil) throw std::runtime_error("No Metal device");
    GPU->Queue = [GPU->Gpu newCommandQueue];
    NSError* Error = nil;
    id<MTLLibrary> Library = [GPU->Gpu newLibraryWithSource:[NSString stringWithUTF8String:kKernels] options:nil error:&Error];
    if (Library == nil) {
      throw std::runtime_error(std::string("Metal kernels failed to compile: ") + [[Error localizedDescription] UTF8String]);
    }
    GPU->SigProducts = GPU->Pipeline(Library, @"sig_products");
    GPU->BinSums = GPU->Pipeline(Library, @"bin_sums");
    GPU->SigShare = GPU->Pipeline(Library, @"sig_share");
    GPU->ParamGradChunks = GPU->Pipeline(Library, @"param_grad_chunks");
    GPU->ParamGradSum = GPU->Pipeline(Library, @"param_grad_sum");

    GPU->Pool = MakeBuffer<float>(GPU->Gpu, Model.PoolPointer.size());
    GPU->SigStart = UploadBuffer<int>(GPU->Gpu, Model.SigStart);
    GPU->SigPool = UploadBuffer<int>(GPU->Gpu, Model.SigPool);
    GPU->SigProd = MakeBuffer<float>(GPU->Gpu, static_cast<std::size_t>(Model.NSignatures()));
    GPU->MC = MakeBuffer<float>(GPU->Gpu, static_cast<std::size_t>(Model.NBins));
    GPU->W2 = MakeBuffer<float>(GPU->Gpu, static_cast<std::size_t>(Model.NBins));
    GPU->BinFactor = MakeBuffer<float>(GPU->Gpu, static_cast<std::size_t>(Model.NBins));
    GPU->Share = MakeBuffer<float>(GPU->Gpu, static_cast<std::size_t>(Model.NSignatures()));
  }
}

// ************************************************
MetalReweight::~MetalReweight() = default;
// ************************************************

// ************************************************
std::string MetalReweight::Summary() const {
// ************************************************
  std::ostringstream Out;
  Out << CpuReweight::Summary() << "; cached passes on " << [[GPU->Gpu name] UTF8String] << " in single precision";
  return Out.str();
}

// ************************************************
void MetalReweight::Fill(const std::vector<double>& ShiftKey, const bool FillW2, double* MC, double* W2) {
// ************************************************
  // The shifts moved: a full pass, on the CPU, through the experiment's own code
  if (!BinsCurrent(ShiftKey)) {
    LastFillOnGPU = false;
    CpuReweight::Fill(ShiftKey, FillW2, MC, W2);
    return;
  }
  @autoreleasepool {
    EnsureGroups();
    if (GPU->GroupsVersion != GroupsVersion) {
      // Upload the groups, with each group's bin and every signature's groups for the gradient
      std::vector<int> GroupBin, SigGroupStart, SigGroup;
      IndexGroups(GroupBin, SigGroupStart, SigGroup);
      GPU->BinStart = UploadBuffer<int>(GPU->Gpu, Groups.BinStart);
      GPU->GroupSig = UploadBuffer<int>(GPU->Gpu, Groups.Sig);
      GPU->GroupConst = UploadBuffer<float>(GPU->Gpu, Groups.Const);
      GPU->GroupConst2 = UploadBuffer<float>(GPU->Gpu, Groups.Const2);
      GPU->GroupBin = UploadBuffer<int>(GPU->Gpu, GroupBin);
      GPU->SigGroupStart = UploadBuffer<int>(GPU->Gpu, SigGroupStart);
      GPU->SigGroup = UploadBuffer<int>(GPU->Gpu, SigGroup);
      GPU->GroupsVersion = GroupsVersion;
    }
    // Read the pool straight into the GPU's buffer
    float* Pool = static_cast<float*>([GPU->Pool contents]);
    const int NPool = static_cast<int>(Model.PoolPointer.size());
    const M3::float_t* const* Pointer = Model.PoolPointer.data();
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int k = 0; k < NPool; ++k) {
      PoolValue[static_cast<std::size_t>(k)] = static_cast<double>(*Pointer[k]);
      Pool[k] = static_cast<float>(PoolValue[static_cast<std::size_t>(k)]);
    }
    GPU->Run({{GPU->SigProducts, {GPU->Pool, GPU->SigStart, GPU->SigPool, GPU->SigProd}, static_cast<unsigned int>(Model.NSignatures())},
              {GPU->BinSums, {GPU->BinStart, GPU->GroupSig, GPU->GroupConst, GPU->GroupConst2, GPU->SigProd, GPU->MC, GPU->W2},
               static_cast<unsigned int>(Model.NBins), true}});
    const float* GpuMC = static_cast<const float*>([GPU->MC contents]);
    const float* GpuW2 = static_cast<const float*>([GPU->W2 contents]);
    for (int b = 0; b < Model.NBins; ++b) {
      MC[b] += static_cast<double>(GpuMC[b]);
      if (FillW2) W2[b] += static_cast<double>(GpuW2[b]);
    }
  }
  LastFillOnGPU = true;
}

// ************************************************
void MetalReweight::AttachDerivatives(const std::vector<WeightDerivativeSource>& Sources, const int NPars_) {
// ************************************************
  CpuReweight::AttachDerivatives(Sources, NPars_);
  @autoreleasepool {
    // Every term of the gradient, grouped by its parameter: (signature, derivative, pool entry)
    std::vector<int> ParStart, TermSig, TermDeriv, TermPool;
    ParameterTerms(ParStart, TermSig, TermDeriv, TermPool);
    // Cut each parameter's terms into chunks, keeping a parameter's chunks together
    constexpr int kChunk = 4096;
    std::vector<int> ChunkStart, ChunkEnd, ParChunkStart(1, 0);
    for (int p = 0; p < NPars; ++p) {
      const int First = ParStart[static_cast<std::size_t>(p)], Last = ParStart[static_cast<std::size_t>(p) + 1];
      for (int k = First; k < Last; k += kChunk) {
        ChunkStart.push_back(k);
        ChunkEnd.push_back(std::min(k + kChunk, Last));
      }
      ParChunkStart.push_back(static_cast<int>(ChunkStart.size()));
    }
    GPU->NChunks = static_cast<unsigned int>(ChunkStart.size());
    GPU->ChunkStart = UploadBuffer<int>(GPU->Gpu, ChunkStart);
    GPU->ChunkEnd = UploadBuffer<int>(GPU->Gpu, ChunkEnd);
    GPU->ParChunkStart = UploadBuffer<int>(GPU->Gpu, ParChunkStart);
    GPU->ChunkSum = MakeBuffer<float>(GPU->Gpu, ChunkStart.size());
    GPU->TermSig = UploadBuffer<int>(GPU->Gpu, TermSig);
    GPU->TermDeriv = UploadBuffer<int>(GPU->Gpu, TermDeriv);
    GPU->TermPool = UploadBuffer<int>(GPU->Gpu, TermPool);
    GPU->NDeriv = PoolDerivative.size();
    GPU->Deriv = MakeBuffer<float>(GPU->Gpu, GPU->NDeriv);
    GPU->Grad = MakeBuffer<float>(GPU->Gpu, static_cast<std::size_t>(NPars));
  }
}

// ************************************************
std::vector<double> MetalReweight::Gradient(const std::vector<double>& BinFactor) {
// ************************************************
  if (!LastFillOnGPU) return CpuReweight::Gradient(BinFactor);
  @autoreleasepool {
    float* Factor = static_cast<float*>([GPU->BinFactor contents]);
    for (std::size_t b = 0; b < BinFactor.size(); ++b) Factor[b] = static_cast<float>(BinFactor[b]);
    float* Deriv = static_cast<float*>([GPU->Deriv contents]);
    for (std::size_t j = 0; j < GPU->NDeriv; ++j) Deriv[j] = (PoolDerivative[j] == nullptr) ? 1.0f : static_cast<float>(*PoolDerivative[j]);
    GPU->Run({{GPU->SigShare, {GPU->SigGroupStart, GPU->SigGroup, GPU->GroupBin, GPU->GroupConst, GPU->SigProd, GPU->BinFactor, GPU->Share},
               static_cast<unsigned int>(Model.NSignatures())},
              {GPU->ParamGradChunks, {GPU->ChunkStart, GPU->ChunkEnd, GPU->TermSig, GPU->TermDeriv, GPU->TermPool, GPU->Deriv, GPU->Pool,
                                      GPU->Share, GPU->ChunkSum}, GPU->NChunks, true},
              {GPU->ParamGradSum, {GPU->ParChunkStart, GPU->ChunkSum, GPU->Grad}, static_cast<unsigned int>(NPars)}});
    const float* Grad = static_cast<const float*>([GPU->Grad contents]);
    return std::vector<double>(Grad, Grad + NPars);
  }
}
