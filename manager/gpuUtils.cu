// C i/o  for printf and others
#include <stdio.h>
#include <vector>

// CUDA specifics

#include <cuda_runtime.h>

#ifdef CUDA_ERROR_CHECK
#include <helper_functions.h>
#include <helper_cuda.h>
#endif

// Define the macros
#define CudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define CudaCheckError()  __cudaCheckError(__FILE__, __LINE__)

/// KS: Need it for shared memory, there is way to use dynamic shared memory but I am lazy right now
#define __BlockSize__ 1024

// CUDA_ERROR_CHECK is now defined in the makefile instead
//#define CUDA_ERROR_CHECK

// **************************************************
//             ERROR CHECKING ROUTINES
// Also exist in helper_cuda.h
// **************************************************

// **************************************************
// Check for a safe call on GPU
inline void __cudaSafeCall( cudaError err, const char *file, const int line ) {
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
/// Check if there's been an error
inline void __cudaCheckError( const char *file, const int line ) {
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
/// KS: Get some fancy info about VRAM usage
inline void checkGpuMem() {
// *******************************************

  float free_m, total_m,used_m;
  size_t free_t, total_t;

  cudaMemGetInfo(&free_t, &total_t);

  free_m = (uint)free_t/1048576.0;
  total_m = (uint)total_t/1048576.0;
  used_m = total_m - free_m;

  printf("  Memory free %f MB, total memory %f MB, memory used %f MB\n", free_m, total_m, used_m);
}

// *******************************************
/// KS: Get some fancy info about GPU
inline void PrintNdevices() {
// *******************************************

  int nDevices;
  cudaGetDeviceCount(&nDevices);

  printf("  Found %i GPUs, currently I only support one GPU\n", nDevices);
}
