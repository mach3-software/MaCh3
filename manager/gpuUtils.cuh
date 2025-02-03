#pragma once


// C i/o  for printf and others
#include <stdio.h>
#include <vector>

// CUDA specifics
// Because CUDA is cuda, need to make sure we don't check C-style floats...
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <cuda_runtime.h>
#pragma GCC diagnostic pop

#ifdef CUDA_ERROR_CHECK
#include <helper_functions.h>
#include <helper_cuda.h>
#endif

// Define the macros
#define CudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define CudaCheckError()  __cudaCheckError(__FILE__, __LINE__)

/// KS: Need it for shared memory, there is way to use dynamic shared memory but I am lazy right now
#define _BlockSize_ 1024


/// @file gpuUtils.cuh
/// @brief Common CUDA utilities and definitions for shared GPU functionality.
/// @author Richard Calland
/// @author Kamil Skwarczynski

/// @todo KS: There is plenty of useful stuff here https://github.com/NVIDIA/cuda-samples/blob/master/Samples/1_Utilities/deviceQuery/deviceQuery.cpp
/// @todo KS: We might want to port some of these utilities, for example having bool if there is unified memory etc.


// **************************************************
//             ERROR CHECKING ROUTINES
// Also exist in helper_cuda.h
// **************************************************

/// @brief Check for a safe call on GPU.
/// @details This function checks the error status returned by CUDA runtime API functions and reports any errors.
/// @param err The CUDA error code.
/// @param file The file name where the error occurred.
/// @param line The line number where the error occurred.
void __cudaSafeCall( cudaError err, const char *file, const int line );

/// @brief Check if there's been an error.
/// @details This function checks if there has been any CUDA runtime API error and reports it.
/// @param file The file name where the error occurred.
/// @param line The line number where the error occurred.
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
/// @param deviceId The ID of the device to be set as active.
void SetDevice(const int deviceId);

/// @brief KS: Get number of GPU threads for currently used GPU
/// @param Device The ID of the device. Defaults to 0.
/// @return The number of GPU threads.
int GetNumGPUThreads(const int Device = 0);
