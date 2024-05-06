// MaCh3 utils for processing/diagnostic MCMC
// Written by Kamil Skwarczynski
//
// Contains code to run on CUDA GPUs. Right now only can calculate autocorrelations
// Potential extensions:
// -Covariance matrix calculations and other matrix operations
// -Effective Sample Size evaluation

#include "manager/gpuUtils.cu"

// ******************************************
// CONSTANTS
// ******************************************

// d_NAME declares DEVICE constants (live on GPU)
__device__ __constant__ int d_nLag;
__device__ __constant__ int d_nDraws;
__device__ __constant__ int d_nEntries;

// h_NAME declares HOST constants (live on CPU)
static int h_nLag     = -1;
static int h_nDraws   = -1;
static int h_nEntries = -1;

// *******************************************
//              INITIALISE GPU
// *******************************************

// *******************************************
/// KS: Initialiser, here we allocate memory for variables and copy constants
__host__ void InitGPU_AutoCorr(
// *******************************************
                          float **ParStep_gpu,
                          float **NumeratorSum_gpu,
                          float **ParamSums_gpu,
                          float **DenomSum_gpu,

                          int n_Entries,
                          int n_Pars,
                          const int n_Lags) {

  // Write to the global statics (h_* denotes host stored variable)
  h_nDraws = n_Pars;
  h_nLag = n_Lags;
  h_nEntries = n_Entries;

  // Copy the constants
  cudaMemcpyToSymbol(d_nLag,   &h_nLag,   sizeof(h_nLag));
  CudaCheckError();

  cudaMemcpyToSymbol(d_nDraws, &h_nDraws, sizeof(h_nDraws));
  CudaCheckError();

  cudaMemcpyToSymbol(d_nEntries, &h_nEntries, sizeof(h_nEntries));
  CudaCheckError();

  // Allocate chunks of memory to GPU
  //Numerator which is directly used for calculating LagL
  cudaMalloc((void **) NumeratorSum_gpu, h_nLag*h_nDraws*sizeof(float));
  CudaCheckError();

  //Denominator which is directly used for calculating LagL
  cudaMalloc((void **) DenomSum_gpu, h_nLag*h_nDraws*sizeof(float));
  CudaCheckError();

  //Mean value for a given parameter
  cudaMalloc((void **) ParamSums_gpu, h_nDraws*sizeof(float));
  CudaCheckError();

  //store value of paramter for each step
  cudaMalloc((void **) ParStep_gpu, h_nDraws*h_nEntries*sizeof(float*));
  CudaCheckError();

  printf(" Allocated in total %f MB for autocorrelations calculations on GPU\n", double(sizeof(float)*(h_nLag*h_nDraws+h_nLag*h_nDraws+h_nDraws+h_nDraws*h_nEntries))/1.E6);

}

// ******************************************************
//                START COPY TO GPU
// ******************************************************

// ******************************************************
/// KS: Copy necessary variables from CPU to GPU
__host__ void CopyToGPU_AutoCorr(
// ******************************************************
                            float *ParStep_cpu,
                            float *NumeratorSum_cpu,
                            float *ParamSums_cpu,
                            float *DenomSum_cpu,

                            float *ParStep_gpu,
                            float *NumeratorSum_gpu,
                            float *ParamSums_gpu,
                            float *DenomSum_gpu) {

  //store value of parameter for each step
  cudaMemcpy(ParStep_gpu, ParStep_cpu, h_nDraws*h_nEntries*sizeof(float), cudaMemcpyHostToDevice);
  CudaCheckError();

  //Mean value for a given parameter
  cudaMemcpy(ParamSums_gpu, ParamSums_cpu, h_nDraws*sizeof(float), cudaMemcpyHostToDevice);
  CudaCheckError();

  //Numerator which is directly used for calculating LagL
  cudaMemcpy(NumeratorSum_gpu, NumeratorSum_cpu, h_nLag*h_nDraws*sizeof(float), cudaMemcpyHostToDevice);
  CudaCheckError();

  //Denominator which is directly used for calculating LagL
  cudaMemcpy(DenomSum_gpu, DenomSum_cpu, h_nLag*h_nDraws*sizeof(float), cudaMemcpyHostToDevice);
  CudaCheckError();
}


// ********************************************************
//                  START GPU KERNELS
//*********************************************************

//*********************************************************
/// Eval autocorrelations based on Box and Jenkins
__global__ void EvalOnGPU_AutoCorr(
    const float* __restrict__ ParStep_gpu,
    const float* __restrict__ ParamSums_gpu,
    float*  NumeratorSum_gpu,
    float*  DenomSum_gpu) {
//*********************************************************

  const unsigned int CurrentLagNum = (blockIdx.x * blockDim.x + threadIdx.x);

  //KS: Accessing shared memory is much much faster than global memory hence we use shared memory for calculation and then write to global memory
  __shared__ float shared_NumeratorSum[_BlockSize_];
  __shared__ float shared_DenomSum[_BlockSize_];

  // this is the stopping condition!
  if (CurrentLagNum < d_nLag*d_nDraws)
  {
      shared_NumeratorSum[threadIdx.x] = 0;
      shared_DenomSum[threadIdx.x] = 0;

      //KS: Might consider caching this information, which MIGHT be faster too lazy right now
      const int Param = int(CurrentLagNum/d_nLag);
      const int nLag = CurrentLagNum - Param*d_nLag;
      // Loop over the number of entries
      for (int i = 0; i < d_nEntries; ++i)
      {
        //KS: Use fmaf to have it tiny bit faster, for something easier to read: Param*d_nEntries + i
        int CurrParStep = fmaf(Param, d_nEntries, i);
        const float Diff = ParStep_gpu[CurrParStep]-ParamSums_gpu[Param];
        // Only sum the numerator up to i = N-k
        if (i < d_nEntries-nLag)
        {
          //KS: Use fmaf to have it tiny bit faster, for something easier to read: Param*d_nEntries + (i + nLag)
          CurrParStep = fmaf(Param, d_nEntries, i + nLag);
          const float LagTerm = ParStep_gpu[CurrParStep]-ParamSums_gpu[Param];
          const float Product = Diff*LagTerm;
          shared_NumeratorSum[threadIdx.x] += Product;
        }
        // Square the difference to form the denominator
        const float Denom = Diff*Diff;
        shared_DenomSum[threadIdx.x] += Denom;
      }

      //KS: Make sure threads are synchronised before moving to global memory
      __syncthreads();
      NumeratorSum_gpu[CurrentLagNum] = shared_NumeratorSum[threadIdx.x];
      DenomSum_gpu[CurrentLagNum]     = shared_DenomSum[threadIdx.x];
  }
}

// *****************************************
/// KS: This call the main kernel responsible for calculating LagL and later copy results back to CPU
__host__ void RunGPU_AutoCorr(
    float*  ParStep_gpu,
    float*  ParamSums_gpu,
    float*  NumeratorSum_gpu,
    float*  DenomSum_gpu,
    float*  NumeratorSum_cpu,
    float*  DenomSum_cpu) {
// *****************************************

  dim3 block_size;
  dim3 grid_size;

  block_size.x = _BlockSize_;
  grid_size.x = (h_nLag*h_nDraws / block_size.x) + 1;

  EvalOnGPU_AutoCorr<<<grid_size, block_size>>>(
      ParStep_gpu,
      ParamSums_gpu,
      NumeratorSum_gpu,
      DenomSum_gpu);
  CudaCheckError();

  printf(" Finished calculating now copying results back to CPU \n");

  //KS: Finally copy paste memory from GPU to CPU
  cudaMemcpy(NumeratorSum_cpu, NumeratorSum_gpu, h_nLag*h_nDraws*sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();

  cudaMemcpy(DenomSum_cpu, DenomSum_gpu, h_nLag*h_nDraws*sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();
}

// *********************************
// CLEANING
// *********************************

// *********************************
/// KS: free memory on gpu
__host__ void CleanupGPU_AutoCorr(
    float *ParStep_gpu,
    float *NumeratorSum_gpu,
    float *ParamSums_gpu,
    float *DenomSum_gpu) {
// *********************************
  cudaFree(ParStep_gpu);
  cudaFree(NumeratorSum_gpu);
  cudaFree(ParamSums_gpu);
  cudaFree(DenomSum_gpu);

  printf(" Cleared memory at GPU, I am free \n");
  return;
}
