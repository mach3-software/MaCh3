#pragma once
// MaCh3 utils for processing/diagnostic MCMC
// Written by Kamil Skwarczynski
//
// Contains code to run on CUDA GPUs. Right now only can calculate autocorrelations
// Potential extensions:
// -Covariance matrix calculations and other matrix operations
// -Effective Sample Size evaluation

#include "manager/gpuUtils.cuh"

// *******************************************
//              INITIALISE GPU
// *******************************************

/// @brief KS: Initialiser, here we allocate memory for variables and copy constants
__host__ void InitGPU_AutoCorr(
                          float **ParStep_gpu,
                          float **NumeratorSum_gpu,
                          float **ParamSums_gpu,
                          float **DenomSum_gpu,

                          int n_Entries,
                          int n_Pars,
                          const int n_Lags);

/// @brief KS: Copy necessary variables from CPU to GPU
__host__ void CopyToGPU_AutoCorr(
                            float *ParStep_cpu,
                            float *NumeratorSum_cpu,
                            float *ParamSums_cpu,
                            float *DenomSum_cpu,

                            float *ParStep_gpu,
                            float *NumeratorSum_gpu,
                            float *ParamSums_gpu,
                            float *DenomSum_gpu);


/// @brief Eval autocorrelations based on Box and Jenkins
__global__ void EvalOnGPU_AutoCorr(
    const float* __restrict__ ParStep_gpu,
    const float* __restrict__ ParamSums_gpu,
    float*  NumeratorSum_gpu,
    float*  DenomSum_gpu);

/// @brief KS: This call the main kernel responsible for calculating LagL and later copy results back to CPU
__host__ void RunGPU_AutoCorr(
    float*  ParStep_gpu,
    float*  ParamSums_gpu,
    float*  NumeratorSum_gpu,
    float*  DenomSum_gpu,
    float*  NumeratorSum_cpu,
    float*  DenomSum_cpu);

/// @brief KS: free memory on gpu
__host__ void CleanupGPU_AutoCorr(
    float *ParStep_gpu,
    float *NumeratorSum_gpu,
    float *ParamSums_gpu,
    float *DenomSum_gpu);
