#pragma once

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
#define _BlockSize_ 1024

//KS: TODO
// There is plenty of useful stuff here https://github.com/NVIDIA/cuda-samples/blob/master/Samples/1_Utilities/deviceQuery/deviceQuery.cpp
// We might want to port some of these utilities, for example having bool if there is unified memory etc.

// CUDA_ERROR_CHECK is now defined in the makefile instead
//#define CUDA_ERROR_CHECK

// **************************************************
//             ERROR CHECKING ROUTINES
// Also exist in helper_cuda.h
// **************************************************

/// @brief Check for a safe call on GPU
void __cudaSafeCall( cudaError err, const char *file, const int line );

/// @brief Check if there's been an error
void __cudaCheckError( const char *file, const int line );

// *******************************************
//              Utils
// *******************************************

// *******************************************
/// @brief KS: Get some fancy info about VRAM usage
void checkGpuMem();

/// @brief KS: Get some fancy info about GPU
void PrintNdevices();

/// @brief KS: Completely clean GPU, this is time consuming and may lead to unexpected behaviour.
void ResetDevice();

/// @brief KS: Only useful if using multiple GPU
void SetDevice(const int deviceId);

/// @brief KS: Get number of GPU threads for currently used GPU
int GetNumGPUThreads(const int Device = 0);