// MaCh3 includes
#include "manager/gpuUtils.cuh"

// **************************************************
// Check for a safe call on GPU
void __cudaSafeCall( cudaError err, const char *file, const int line ) {
// **************************************************
#ifdef CUDA_ERROR_CHECK
  if (cudaSuccess != err) {
    fprintf(stderr, "cudaSafeCall() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
    exit(-1);
  }
#endif
  return;
}

// **************************************************
// Check if there's been an error
void __cudaCheckError( const char *file, const int line ) {
// **************************************************
#ifdef CUDA_ERROR_CHECK
  cudaError err = cudaGetLastError();
  if (cudaSuccess != err) {
    fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
    exit(-1);
  }

  // More careful checking. However, this will affect performance.
  // Comment away if needed.
  err = cudaDeviceSynchronize();
  if (cudaSuccess != err) {
    fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
    exit(-1);
  }
#endif
  return;
}

// *******************************************
//              Utils
// *******************************************

// *******************************************
// KS: Get some fancy info about VRAM usage
void checkGpuMem() {
// *******************************************

  float free_m, total_m,used_m;
  size_t free_t, total_t;

  cudaMemGetInfo(&free_t, &total_t);
  CudaCheckError();

  free_m = (uint)free_t/1048576.0;
  total_m = (uint)total_t/1048576.0;
  used_m = total_m - free_m;

  printf("  Memory free %f MB, total memory %f MB, memory used %f MB\n", free_m, total_m, used_m);
}

// *******************************************
// KS: Get some fancy info about GPU
void PrintNdevices() {
// *******************************************

  int nDevices;
  cudaGetDeviceCount(&nDevices);
  CudaCheckError();

  if (nDevices == 0) {
    printf("No CUDA devices found");
    throw;
  }

  printf("  Found %i GPUs, currently I only support one GPU\n", nDevices);
}


// *******************************************
// KS: Completely clean GPU, this is time consuming and may lead to unexpected behaviour.
void ResetDevice() {
// *******************************************

  cudaDeviceReset();
  CudaCheckError();
}


// *******************************************
/// @brief Only useful if using multiple GPU
void SetDevice(const int deviceId) {
// *******************************************

  // Check if the device ID is valid
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceId < 0 || deviceId >= deviceCount) {
    printf("Invalid device ID: %i \n", deviceId);
    throw;
  }

  cudaSetDevice(deviceId);
  CudaCheckError();
  printf("GPU device set to ID: %i \n", deviceId);

}

// *******************************************
// Get number of GPU threads for currently used GPU
int GetNumGPUThreads(const int Device) {
// *******************************************

  int deviceCount;
  cudaGetDeviceCount(&deviceCount);

  if (deviceCount == 0) {
    printf("No CUDA devices found");
    throw;
  }

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, Device);

  // Define the number of threads per block
  int nThreadsBlocks = (deviceProp.multiProcessorCount * deviceProp.maxThreadsPerMultiProcessor);

  return nThreadsBlocks;
}
