#pragma once
// MaCh3 utils for processing/diagnostic MCMC
// Written by Kamil Skwarczynski
//
// Contains code to run on CUDA GPUs. Right now only can calculate autocorrelations
// Potential extensions:
// -Covariance matrix calculations and other matrix operations
// -Effective Sample Size evaluation

#include "Manager/gpuUtils.cuh"

/// @file gpuMCMCProcessorUtils.cuh
/// @author Kamil Skwarczynski

/// @brief Eval autocorrelations based on Box and Jenkins
/// @param ParStep_gpu Parameter value at each step
/// @param NumeratorSum_gpu Sum used for nominator of autocorrelation calculations
/// @param ParamSums_gpu Overall sum for each parameter over all steps
/// @param DenomSum_gpu Sum used for denominator of autocorrelation calculations
__global__ void EvalOnGPU_AutoCorr(
  const float* __restrict__ ParStep_gpu,
  const float* __restrict__ ParamSums_gpu,
  float*  NumeratorSum_gpu,
  float*  DenomSum_gpu);


/// @brief Class responsible for performing MCMC Processing with GPU
/// @author Kamil Skwarczynski
class MCMCProcessorGPU
{
 public:
  /// @brief constructor
  MCMCProcessorGPU();
  /// @brief destructor
  virtual ~MCMCProcessorGPU();

  /// @brief KS: Initialiser, here we allocate memory for variables and copy constants
  /// @param n_Entries Total number of entries in mcmc chain
  /// @param n_Pars Number of relevant parameters
  /// @param n_Lags Value of Lag in autocreation calculation
  __host__ void InitGPU_AutoCorr(
                            int n_Entries,
                            int n_Pars,
                            const int n_Lags);

  /// @brief KS: Copy necessary variables from CPU to GPU
  /// @param ParStep_cpu Parameter value at each step
  /// @param NumeratorSum_cpu Sum used for nominator of autocorrelation calculations
  /// @param ParamSums_cpu Overall sum for each parameter over all steps
  /// @param DenomSum_cpu Sum used for denominator of autocorrelation calculations
  __host__ void CopyToGPU_AutoCorr(
                              float *ParStep_cpu,
                              float *NumeratorSum_cpu,
                              float *ParamSums_cpu,
                              float *DenomSum_cpu);

  /// @brief KS: This call the main kernel responsible for calculating LagL and later copy results back to CPU
  /// @param NumeratorSum_cpu Sum used for nominator of autocorrelation calculations
  /// @param DenomSum_cpu Sum used for denominator of autocorrelation calculations
  __host__ void RunGPU_AutoCorr(float*  NumeratorSum_cpu,
                                float*  DenomSum_cpu);

  /// @brief KS: free memory on gpu
  __host__ void CleanupGPU_AutoCorr();


 private:
  /// Value of each param at GPU
  float* ParStep_gpu;
  float* NumeratorSum_gpu;
  float* ParamSums_gpu;
  float* DenomSum_gpu;

  // h_NAME declares HOST constants (live on CPU)
  int h_nLag;
  int h_nDraws;
  int h_nEntries;
};
