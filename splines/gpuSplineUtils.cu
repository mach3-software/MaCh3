// MaCh3 event-by-event cross-section spline code
// Written by Richard Calland, Asher Kaboth, Clarence Wret, Kamil Skwarczynski
//
// Contains code to run on CUDA GPUs. Essentially we load up stripped TSpline3 objects to the GPU and do the equivalent of TSpline3->Eval(double) for all events
// Now also supports TF1 evals
// Called from samplePDFN/samplePDFND.cpp -> splines/SplineMonolith.cpp -> splines/gpuSplineUtils.cu

// C i/o  for printf and others
#include <stdio.h>

// CUDA specifics
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

// Define the macros
#define CudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define CudaCheckError()  __cudaCheckError(__FILE__, __LINE__)

// Hard code the number of splines
// Not entirely necessary: only used for val_gpu and segment_gpu being device constants. Could move them to not being device constants
#define __N_SPLINES__ 48

// CUDA_ERROR_CHECK is now defined in the makefile instead
//#define CUDA_ERROR_CHECK
//#define Weight_On_SplineBySpline_Basis

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
// Check if there's been an error
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

// ******************************************
// CONSTANTS
// ******************************************

// d_NAME declares DEVICE constants (live on GPU)
__device__ __constant__ int d_n_splines;
__device__ __constant__ int d_spline_size;
#ifndef Weight_On_SplineBySpline_Basis
__device__ __constant__ int d_n_params;
__device__ __constant__ int d_n_events;
#endif
// Constant memory needs to be hard-coded on compile time
// Could make this texture memory instead, but don't care enough right now...
__device__ __constant__ float val_gpu[__N_SPLINES__];
__device__ __constant__ int segment_gpu[__N_SPLINES__];

// h_NAME declares HOST constants (live on CPU)
static int h_n_splines    = -1;
static int h_spline_size  = -1;
static int h_n_params     = -1;
#ifndef Weight_On_SplineBySpline_Basis
static int h_n_events     = -1;
#endif
// *******************************************
//              INITIALISE GPU
// *******************************************

// *******************************************
// Initaliser when using the x array and combined y,b,c,d array
__host__ void InitGPU_SepMany(
// *******************************************
                          float **gpu_x_array,
                          float **gpu_many_array,
                          float **gpu_weights,

                          int** gpu_paramNo_arr,
                 #ifndef Weight_On_SplineBySpline_Basis
                          float **gpu_total_weights,
                          int n_events,

                          int** index_gpu,
                          unsigned int NSplines_total_large,
                  #endif
                          int sizeof_array,
                          int n_splines) {

  // Allocate chunks of memory to GPU
  cudaMalloc((void **) gpu_paramNo_arr, n_splines*sizeof(int));
  CudaCheckError();

  cudaMalloc((void **) gpu_x_array, sizeof_array*sizeof(float));
  CudaCheckError();

  cudaMalloc((void **) gpu_many_array, 4*sizeof_array*sizeof(float));
  CudaCheckError();

  // Allocate memory for the array of weights to be returned to CPU
  cudaMalloc((void **) gpu_weights, n_splines*sizeof(float));
  CudaCheckError();
#ifndef Weight_On_SplineBySpline_Basis
  //KS: Allocate memory for the array of total weights to be returned to CPU
  cudaMalloc((void **) gpu_total_weights, n_events*sizeof(float));
  CudaCheckError();

  //KS: Allocate memory for the index of GPU splines to global spline number
  cudaMalloc((void **) index_gpu, NSplines_total_large*sizeof(int));
  CudaCheckError();
#endif

  // Print allocation info to user
  printf("  Allocated %i entries for paramNo arrays, size = %f MB\n", n_splines, double(sizeof(int)*n_splines)/1.E6);
  printf("  Allocated %i entries for x coeff arrays, size = %f MB\n", sizeof_array, double(sizeof(float)*sizeof_array)/1.E6);
  printf("  Allocated %i entries for {ybcd} coeff arrays, size = %f MB\n", 4*sizeof_array, double(sizeof(float)*4*sizeof_array)/1.E6);
}

// *******************************************
// Initaliser when using the x array and combined y,b,c,d array
__host__ void InitGPU_TF1(
// *******************************************
                          float **gpu_coeffs,
                          int** gpu_paramNo_arr,
                          int** gpu_nPoints_arr,
                          float **gpu_weights,

                    #ifndef Weight_On_SplineBySpline_Basis
                          float **gpu_total_weights,
                          int n_events,

                          int** index_gpu,
                          unsigned int NSplines_total_large,
                    #endif
                          int n_splines) {

  // Holds the parameter number
  cudaMalloc((void **) gpu_paramNo_arr, n_splines*sizeof(int));
  CudaCheckError();

  // Holds the number of points
  cudaMalloc((void **) gpu_nPoints_arr, n_splines*sizeof(int));
  CudaCheckError();

  // Holds the coefficients (5th order polynomial and constant term == 1) -> 5
  cudaMalloc((void **) gpu_coeffs, 5*n_splines*sizeof(float));
  CudaCheckError();

  // Allocate memory for the array of weights to be returned to CPU
  cudaMalloc((void **) gpu_weights, n_splines*sizeof(float));
  CudaCheckError();

#ifndef Weight_On_SplineBySpline_Basis
  //KS: Allocate memory for the array of total weights to be returned to CPU
  cudaMalloc((void **) gpu_total_weights, n_events*sizeof(float));
  CudaCheckError();

  //KS: Allocate memory for the index of GPU splines to global spline number
  cudaMalloc((void **) index_gpu, NSplines_total_large*sizeof(int));
  CudaCheckError();
#endif

  // Print allocation info to user
  printf("  Allocated %i entries for paramNo and nPoints arrays, size = %f MB\n", n_splines, double(2.0*sizeof(int)*n_splines)/1.E6);
  printf("  Allocated %i entries for coefficient arrays, size = %f MB\n", 5*n_splines, double(sizeof(float)*5*n_splines)/1.E6);

}


// ******************************************************
//                START COPY TO GPU
// ******************************************************

// ******************************************************
// Copy to GPU for x array and separate ybcd array
__host__ void CopyToGPU_SepMany(
// ******************************************************
                            int *gpu_paramNo_arr,
                            float *gpu_x_array,
                            float *gpu_many_array,

                            int *paramNo_arr,
                            float *cpu_x_array,
                            float *cpu_many_array,

                    #ifndef Weight_On_SplineBySpline_Basis
                            int *index_cpu,
                            int *index_gpu,
                            unsigned int NSplines_total_large,
                            int n_events,
                    #endif
                            int n_params,
                            int n_splines,
                            int spline_size) {
  if (n_params != __N_SPLINES__) {
    printf("Number of splines not equal to %i, GPU code for event-by-event splines will fail\n", __N_SPLINES__);
    printf("n_params = %i\n", n_params);
    printf("%s : %i\n", __FILE__, __LINE__);
    exit(-1);
  }

  // Write to the global statics (h_* denotes host stored variable)
  h_n_params = n_params;
  h_n_splines = n_splines;
  h_spline_size = spline_size;
#ifndef Weight_On_SplineBySpline_Basis
  h_n_events    = n_events;
#endif
  // Copy the constants
  // Total number of valid splines for all loaded events
  cudaMemcpyToSymbol(d_n_splines,   &h_n_splines,   sizeof(h_n_splines));
  CudaCheckError();
  // Total spline size per spline; i.e. just the number of points or knots in the spline
  cudaMemcpyToSymbol(d_spline_size, &h_spline_size, sizeof(h_spline_size));
  CudaCheckError();
#ifndef Weight_On_SplineBySpline_Basis
  // Number of dials
  cudaMemcpyToSymbol(d_n_params, &n_params, sizeof(n_params));
  CudaCheckError();

  // Number of events
  cudaMemcpyToSymbol(d_n_events, &h_n_events, sizeof(h_n_events));
  CudaCheckError();
#endif
  // Copy the coefficient arrays to the GPU; this only happens once per entire Markov Chain so is OK to do multiple extensive memory copies
  size_t coeff_size = sizeof(float)*(n_splines * spline_size);

  cudaMemcpy(gpu_many_array, cpu_many_array, coeff_size*4, cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_x_array, cpu_x_array, coeff_size, cudaMemcpyHostToDevice);
  CudaCheckError();

  // Also copy the parameter number for each spline onto the GPU; i.e. what spline parameter are we calculating right now
  cudaMemcpy(gpu_paramNo_arr, paramNo_arr, h_n_splines*sizeof(int), cudaMemcpyHostToDevice);
  CudaCheckError();
  #ifndef Weight_On_SplineBySpline_Basis
  //KS: Lastly copy the index to global spline number
  cudaMemcpy(index_gpu, index_cpu, NSplines_total_large*sizeof(int), cudaMemcpyHostToDevice);
  CudaCheckError();
  #endif
}

// ******************************************************
// Copy to GPU for x array and separate ybcd array
__host__ void CopyToGPU_TF1(
// ******************************************************
                            float *gpu_coeffs,
                            int *gpu_paramNo_arr,
                            int *gpu_nPoints_arr,

                            float *cpu_coeffs,
                            int *paramNo_arr,
                            int *nPoints_arr,

                  #ifndef Weight_On_SplineBySpline_Basis
                            int *index_cpu,
                            int *index_gpu,
                            unsigned int NSplines_total_large,
                            int n_events,
                  #endif

                            int n_params,
                            int n_splines,
                            int _max_knots) {

  if (n_params != __N_SPLINES__) {
    printf("Number of splines not equal to %i, GPU code for event-by-event splines will fail\n", __N_SPLINES__);
    printf("n_params = %i\n", n_params);
    printf("%s : %i\n", __FILE__, __LINE__);
    exit(-1);
  }

  // Write to the global statics (h_* denotes host stored variable)
  h_n_params = n_params;
  h_n_splines = n_splines;
  h_spline_size = _max_knots;
#ifndef Weight_On_SplineBySpline_Basis
  h_n_events    = n_events;
#endif
  // Copy the constants
  // Total number of valid splines for all loaded events
  cudaMemcpyToSymbol(d_n_splines,   &h_n_splines,   sizeof(h_n_splines));
  CudaCheckError();
  // Total spline size per spline; i.e. just the number of points or knots in the spline
  cudaMemcpyToSymbol(d_spline_size, &h_spline_size, sizeof(h_spline_size));
  CudaCheckError();

#ifndef Weight_On_SplineBySpline_Basis
  // Number of dials
  cudaMemcpyToSymbol(d_n_params, &n_params, sizeof(n_params));
  CudaCheckError();

  // Number of events
  cudaMemcpyToSymbol(d_n_events, &h_n_events, sizeof(h_n_events));
  CudaCheckError();
#endif
  // Move the coefficients
  cudaMemcpy(gpu_coeffs, cpu_coeffs, n_splines*5*sizeof(float), cudaMemcpyHostToDevice);
  CudaCheckError();

  // Also copy the parameter number for each spline onto the GPU; i.e. what spline parameter are we calculating right now
  cudaMemcpy(gpu_paramNo_arr, paramNo_arr, h_n_splines*sizeof(int), cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_nPoints_arr, nPoints_arr, h_n_splines*sizeof(int), cudaMemcpyHostToDevice);
  CudaCheckError();

  #ifndef Weight_On_SplineBySpline_Basis
  //KS: Lastly copy the index to global spline number
  cudaMemcpy(index_gpu, index_cpu, NSplines_total_large*sizeof(int), cudaMemcpyHostToDevice);
  CudaCheckError();
  #endif
}

// ********************************************************
//                  START GPU KERNELS
//*********************************************************
// All the GPU kernels have similar tasks but different implementations
// Essentially they perform a binary search to find which TSpline3 point is nearest to our parameter variation
// Once it knows this, we simply extract the pre-computed coefficients for that spline point and multiply together to get a weight

//*********************************************************
// Evaluate the spline on the GPU
// Using one {y,b,c,d} array
// And one {x} array
// Should be most efficient at cache hitting and memory coalescense
__global__ void EvalOnGPU_SepMany(const int* __restrict__ gpu_paramNo_arr,
//*********************************************************
    const float* __restrict__ gpu_coeff_x,
    const float* __restrict__ gpu_coeff_many,
    float *gpu_weights) {
//*********************************************************

  // points per spline is the offset to skip in the index to move between splines
  const unsigned int splineNum = (blockIdx.x * blockDim.x + threadIdx.x);

  // Note, the new arrays are arranged as:
  //       gpu_nPoints_arr has length = spln_counter (total number of valid splines for all events in the loop)
  //       gpu_coeff_x has length = spln_counter * spline_size
  //       gpu_coeff_many has length = spln_counter * spline_size * 4
  //       ...
  //       gpu_weights has length = spln_counter * spline_size

  // this is the stopping condition!
  if (splineNum < d_n_splines) {

    // This is the parameter number for this particular splineNum; 0 = MACCQE, 1 = pFC, 2 = EBC, etc (see SplineMonolith constructor for this numbering scheme; should probably be made into an enum!)
    //const unsigned int paramNo = gpu_paramNo_arr[splineNum];
    //const float paramVal = val_gpu[paramNo];
    const float paramVal = val_gpu[gpu_paramNo_arr[splineNum]];

    // Number of points in TSpline3 already saved in d_spline_size

    // Seg (segment) is the part of the TSpline3 that needs to be
    // evaluated. Determining the segment determines the subset of
    // parameters used to evaluate the polynomial.

    // val[splineNum] is the variation in x for this spline
    // gpu_coeff_x[splineNum] is the first x-value; if variation is beneath first point, choose zeroth segment
    unsigned int klow = 0;

    if (paramVal <= gpu_coeff_x[splineNum*d_spline_size]) {
      klow = 0;
      // gpu_coeff_x[splineNum+d_spline_size-1] is the last x-value; if the variation is above this point, choose high segment
    } else if (paramVal >= gpu_coeff_x[splineNum*d_spline_size+d_spline_size-1]) {
      // This is DEFINITELY correct; TSpline3 put d_spline_size-1 but then changes to d_spline_size-2, see http://savannah.web.cern.ch/savannah/HEP_Applications/savroot/bugs/71651.html
      klow = d_spline_size - 2;
    } else {
      unsigned int khig = d_spline_size-1;
      unsigned int khalf = 0;

      // Try to not unroll here because of register pressure
//#pragma unroll 4
      while (khig - klow > 1) {
        khalf = (klow + khig)/2;
        if (paramVal > gpu_coeff_x[splineNum*d_spline_size + khalf]) {
          klow = khalf;
        } else {
          khig = khalf;
        }
      }
    }

    // After the segment is determined (klow) grab the appropriate
    // polynomial parameters from the monolithic splineMonolith
    const float fY = gpu_coeff_many[splineNum*d_spline_size*4+klow*4];
    const float fB = gpu_coeff_many[splineNum*d_spline_size*4+klow*4+1];
    const float fC = gpu_coeff_many[splineNum*d_spline_size*4+klow*4+2];
    const float fD = gpu_coeff_many[splineNum*d_spline_size*4+klow*4+3];
    const float dx = paramVal - gpu_coeff_x[splineNum*d_spline_size+klow];

    // Wooow, let's use some fancy intrinsics and pull down the processing time by <1% from normal multiplication! HURRAY
    gpu_weights[splineNum] = fmaf(dx, fmaf(dx, fmaf(dx, fD, fC), fB), fY);
    // Or for the more "easy to read" version: gpu_weights[splineNum] = (fY+dx*(fB+dx*(fC+dx*fD)));
#ifdef DEBUG
    printf("splineNum = %i/%i, paramNo = %i, paramval = %f, fX = %f, fX+1 = %f, dx = %f, klow = %i, d_n_splines = %i, d_spline_size = %i\n", splineNum, d_n_splines, gpu_paramNo_arr[splineNum], paramVal, gpu_coeff_x[splineNum*d_spline_size+klow], gpu_coeff_x[splineNum*d_spline_size+klow+1], dx, klow, d_n_splines, d_spline_size);
#endif

  }
}

//*********************************************************
// Evaluate the spline on the GPU
// Using one {y,b,c,d} array
// And one {x} array
// Should be most efficient at cache hitting and memory coalescense
// But using spline segments rather than the parameter value: avoids doing binary search on GPU
__global__ void EvalOnGPU_SepMany_seg( const int* __restrict__ gpu_paramNo_arr,
    const float* __restrict__ gpu_coeff_x,
    const float* __restrict__ gpu_coeff_many,
    float *gpu_weights) {
//*********************************************************

  // points per spline is the offset to skip in the index to move between splines
  const unsigned int splineNum = (blockIdx.x * blockDim.x + threadIdx.x);

  // Note, the new arrays are arranged as:
  //       gpu_paramNo_arr has length = spln_counter (keeps track of which parameter we're using on this thread)
  //       gpu_coeff_x has length = spln_counter * spline_size
  //       gpu_coeff_many has length = spln_counter * spline_size * 4
  //       ...
  //       gpu_weights has length = spln_counter * spline_size

  // this is the stopping condition!
  if (splineNum < d_n_splines) {
    // This is the segment we want for this parameter variation
    // for this particular splineNum; 0 = MACCQE, 1 = pFC, 2 = EBC, etc
    // Avoids doing costly binary search on GPU
    const int segment = segment_gpu[gpu_paramNo_arr[splineNum]];

    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float variation = val_gpu[gpu_paramNo_arr[splineNum]];

    // We've read the segment straight from CPU and is saved in segment_gpu
    // polynomial parameters from the monolithic splineMonolith
    const float fY = gpu_coeff_many[splineNum*d_spline_size*4+segment*4];
    const float fB = gpu_coeff_many[splineNum*d_spline_size*4+segment*4+1];
    const float fC = gpu_coeff_many[splineNum*d_spline_size*4+segment*4+2];
    const float fD = gpu_coeff_many[splineNum*d_spline_size*4+segment*4+3];
    const float dx = variation - gpu_coeff_x[splineNum*d_spline_size+segment];

    // Wooow, let's use some fancy intrinsics and pull down the processing time by <1% from normal multiplication! HURRAY
    gpu_weights[splineNum] = fmaf(dx, fmaf(dx, fmaf(dx, fD, fC), fB), fY);
    // Or for the more "easy to read" version:
    //gpu_weights[splineNum] = (fY+dx*(fB+dx*(fC+dx*fD)));

#ifdef DEBUG
  printf("splineNum = %i/%i, paramNo = %i, variation = %f, segment = %i, fX = %f, fX+1 = %f, dx = %f, d_n_splines = %i, d_spline_size = %i\n", splineNum, d_n_splines, gpu_paramNo_arr[splineNum], variation, segment, gpu_coeff_x[splineNum*d_spline_size+segment], gpu_coeff_x[splineNum*d_spline_size+segment+1], dx, d_n_splines, d_spline_size);
#endif
  }
}

//*********************************************************
// Evaluate the TF1 on the GPU
// Using 5th order polynomial
__global__ void EvalOnGPU_TF1(
    const float* __restrict__ gpu_coeffs,
    const int* __restrict__ gpu_paramNo_arr,
    const int* __restrict__ gpu_nPoints_arr,
    float *gpu_weights) {
//*********************************************************

  // points per spline is the offset to skip in the index to move between splines
  const unsigned int splineNum = (blockIdx.x * blockDim.x + threadIdx.x);

  // Note, the new arrays are arranged as:
  //       gpu_paramNo_arr has length = spln_counter (keeps track of which parameter we're using on this thread)
  //       gpu_coeff_x has length = spln_counter * spline_size
  //       gpu_coeff_many has length = spln_counter * spline_size * 4
  //       ...
  //       gpu_weights has length = spln_counter * spline_size

  if (splineNum < d_n_splines) {
    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float x = val_gpu[gpu_paramNo_arr[splineNum]];

    // Read the coefficients
    const float a = gpu_coeffs[splineNum*d_spline_size];
    const float b = gpu_coeffs[splineNum*d_spline_size+1];
    const float c = gpu_coeffs[splineNum*d_spline_size+2];
    const float d = gpu_coeffs[splineNum*d_spline_size+3];
    const float e = gpu_coeffs[splineNum*d_spline_size+4];

    // Match these with form in SetSplines
    // Might not be great to have this if statement: maybe split two kernels?
    if (gpu_nPoints_arr[splineNum] == 5) {
      gpu_weights[splineNum] = 1 + a*x + b*x*x + c*x*x*x + d*x*x*x*x + e*x*x*x*x*x;
    } else if (gpu_nPoints_arr[splineNum] == 2) {
      gpu_weights[splineNum] = (x<=0)*(1+a*x) + (x>0)*(1+b*x);
    } else {
      printf("Big problems, I found a nPoints array which is not 5 or 2 on GPU!\n");
    }

#ifdef DEBUG
    //if (splineNum < 200) {
    if (gpu_nPoints_arr[splineNum] == 2) {
      printf("splineNum = %i, spline_size=%i, paramNo = %i, variation = %f, a = %f, b = %f, c = %f, d = %f, e = %f, weight = %f\n", splineNum, d_spline_size, gpu_paramNo_arr[splineNum], x, a, b, c, d, e, gpu_weights[splineNum] );
    }
#endif
  }

}

#ifndef Weight_On_SplineBySpline_Basis
//*********************************************************
//KS: Evaluate the total spline event weight on the GPU, as in most cases GPU is faster, even more this significant reduce memory transfer from GPU to CPU
__global__ void EvalOnGPU_TotWeight(
   const int* __restrict__ index_gpu,
   const float* __restrict__ gpu_weights,
   float *gpu_total_weights) {
//*********************************************************
    const unsigned int EventNum = (blockIdx.x * blockDim.x + threadIdx.x);
    if(EventNum < d_n_events) //stopping condition
    {
        gpu_total_weights[EventNum] = 1.;
        for (int id = 0; id < d_n_params; id++)
        {
            const unsigned int Global_Spline_Num = EventNum * d_n_params + id;
            if (index_gpu[Global_Spline_Num] >= 0) gpu_total_weights[EventNum] *= gpu_weights[index_gpu[Global_Spline_Num]];
            else gpu_total_weights[EventNum] *=  1.;

            #ifdef DEBUG
            printf("Event = %i, Spline_Num = %i, index_gpu = %i, gpu_weights = %f, total gpu_total_weights = %f \n",
                   EventNum, Global_Spline_Num, index_gpu[Global_Spline_Num], gpu_weights[index_gpu[Global_Spline_Num]],  gpu_total_weights[EventNum]);
            #endif
        }
    }
}
#endif
// *****************************************
// Run the GPU code for the separate many arrays
// As in separate {x}, {y,b,c,d} arrays
// Passing the parameter variations (val) so binary search on GPU
__host__ void RunGPU_SepMany(
    int* gpu_paramNo_arr,

    float *gpu_coeff_x,
    float *gpu_coeff_many,

    float* gpu_weights,
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    int* index_gpu,
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif
    // Holds the changes in parameters, only nParams long now so nParams*sizeof(float) occupied on GPU and transferred
    float* val ) {
  // *****************************************

  dim3 block_size;
  dim3 grid_size;

  block_size.x = 1024;
  grid_size.x = (h_n_splines / block_size.x) + 1;

  // Copy the variation values (val) to the GPU (val_gpu), which is h_n_params long
  cudaMemcpyToSymbol(val_gpu, val, h_n_params*sizeof(float));
  CudaCheckError();

#ifdef DEBUG
  printf("\n***********************\nGPU DEBUGGING ENABLED\n***********************\n");
  printf("block_size.x = %i, grid_size.x = %i \n", block_size.x, grid_size.x);
  printf("RunGPU_SepMany parameter variations\n");
  for (int i = 0; i < h_n_params; i++) {
    printf("val[%i] = %f \n", i, val[i]);
  }
  printf("nParams = %i, n_splines = %i", h_n_params, h_n_splines);
  printf("\n***********************\nAM NOW CALLING KERNEL\n***********************\n");
#endif

  // Set the cache config to prefer L1 for the kernel
  //cudaFuncSetCacheConfig(EvalOnGPU_SepMany, cudaFuncCachePreferL1);

  EvalOnGPU_SepMany<<<grid_size, block_size>>>(
      gpu_paramNo_arr,

      gpu_coeff_x,
      gpu_coeff_many,

      gpu_weights
      );
  CudaCheckError();

#ifdef DEBUG
  printf("Evaluated kernel with SUCCESS (drink beer)\n");
#endif

//KS: We can either copy gpu_weight and calculate total weight in reweighting loop, or not copy and calcualate total weight stall at GPU, which means less memory transfer
#ifdef Weight_On_SplineBySpline_Basis
  // Here we have to make a somewhat large GPU->CPU transfer because it's all the splines' response
  cudaMemcpy(cpu_weights, gpu_weights, h_n_splines*sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();
  #ifdef DEBUG
  printf("Copied GPU weights to CPU with SUCCESS (drink moar beer)\n");
  printf("Released calculated response from GPU with SUCCESS (drink most beer)\n");
  #endif

//KS: Else calculate Total Weight
#else
  grid_size.x = (h_n_events / block_size.x) + 1;

#ifdef DEBUG
  printf("\n***********************\nGPU DEBUGGING ENABLED\n***********************\n");
  printf("block_size.x = %i, grid_size.x = %i \n", block_size.x, grid_size.x);
  printf("RunGPU_TotWeight\n");

  printf("nEvents = %i, n_splines = %i, d_n_params", h_n_events, h_n_splines, h_n_params);
  printf("\n***********************\nI AM NOW CALLING KERNEL\n***********************\n");
#endif
  EvalOnGPU_TotWeight<<<grid_size, block_size>>>(
      index_gpu,
      gpu_weights,
      gpu_total_weights
      );
  #ifdef DEBUG
  CudaCheckError();
  printf("Evaluated kernel with SUCCESS (drink tea)\n");
  #endif
  //KS: Here we have to make a somewhat large GPU->CPU transfer because it is proportional to number of events
  cudaMemcpy(cpu_total_weights, gpu_total_weights, h_n_events*sizeof(float), cudaMemcpyDeviceToHost);
  #ifdef DEBUG
  CudaCheckError();
  printf("Copied GPU total weights to CPU with SUCCESS (drink moar tea)\n");
  printf("Released calculated response from GPU with SUCCESS (drink most tea)\n");
  #endif
#endif
}

// *****************************************
// Run the GPU code for the separate many arrays
// As in separate {x}, {y,b,c,d} arrays
// Pass the segment and the parameter values
// (binary search already performed in samplePDFND::FindSplineSegment()
__host__ void RunGPU_SepMany_seg(
    int* gpu_paramNo_arr,

    float *gpu_coeff_x,
    float *gpu_coeff_many,

    float* gpu_weights,
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    int* index_gpu,
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif
    // Holds the changes in parameters
    float *vals,
    // Holds the segments for parameters
    int *segment ) {
  // *****************************************

  dim3 block_size;
  dim3 grid_size;

  block_size.x = 1024;
  grid_size.x = (h_n_splines / block_size.x) + 1;

  // Copy the segment values to the GPU (segment_gpu), which is h_n_params long
  cudaMemcpyToSymbol(segment_gpu, segment, h_n_params*sizeof(int));
  CudaCheckError();

  // Copy the parameter values values to the GPU (vals_gpu), which is h_n_params long
  cudaMemcpyToSymbol(val_gpu, vals, h_n_params*sizeof(float));
  CudaCheckError();

#ifdef DEBUG
  printf("\n***********************\nGPU DEBUGGING ENABLED\n***********************\n");
  printf("block_size.x = %i, grid_size.x = %i \n", block_size.x, grid_size.x);
  printf("RunGPU_SepMany_seg segments\n");
  for (int i = 0; i < h_n_params; i++) {
    printf("val[%i] = %f in segment %i\n", i, vals[i], segment[i]);
  }
  printf("nParams = %i, n_splines = %i", h_n_params, h_n_splines);
  printf("\n***********************\nAM NOW CALLING KERNEL\n***********************\n");
#endif

  // Set the cache config to prefer L1 for the kernel
  //cudaFuncSetCacheConfig(EvalOnGPU_SepMany, cudaFuncCachePreferL1);
  EvalOnGPU_SepMany_seg<<<grid_size, block_size>>>(
      gpu_paramNo_arr,

      gpu_coeff_x,
      gpu_coeff_many,

      gpu_weights
      );
  CudaCheckError();

#ifdef DEBUG
  printf("Evaluated kernel with SUCCESS (drink beer)\n");
#endif

//KS: We can either copy gpu_weight and calculate total weight in reweighting loop, or not copy and calcualate total weight stall at GPU, which means less memory transfer
#ifdef Weight_On_SplineBySpline_Basis
  // Here we have to make a somewhat large GPU->CPU transfer because it's all the splines' response
  cudaMemcpy(cpu_weights, gpu_weights, h_n_splines*sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();
  #ifdef DEBUG
  printf("Copied GPU weights to CPU with SUCCESS (drink moar beer)\n");
  printf("Released calculated response from GPU with SUCCESS (drink most beer)\n");
  #endif

//KS: Else calculate Total Weight
#else
  grid_size.x = (h_n_events / block_size.x) + 1;

  #ifdef DEBUG
  printf("\n***********************\nGPU DEBUGGING ENABLED\n***********************\n");
  printf("block_size.x = %i, grid_size.x = %i \n", block_size.x, grid_size.x);
  printf("RunGPU_TotWeight\n");

  printf("nEvents = %i, n_splines = %i, d_n_params", h_n_events, h_n_splines, h_n_params);
  printf("\n***********************\nI AM NOW CALLING KERNEL\n***********************\n");
  #endif
  EvalOnGPU_TotWeight<<<grid_size, block_size>>>(
      index_gpu,
      gpu_weights,
      gpu_total_weights
      );
  #ifdef DEBUG
  CudaCheckError();
  printf("Evaluated kernel with SUCCESS (drink tea)\n");
  #endif
  //KS: Here we have to make a somewhat large GPU->CPU transfer because it is proportional to number of events
  cudaMemcpy(cpu_total_weights, gpu_total_weights, h_n_events*sizeof(float), cudaMemcpyDeviceToHost);
  #ifdef DEBUG
  CudaCheckError();
  printf("Copied GPU total weights to CPU with SUCCESS (drink moar tea)\n");
  printf("Released calculated response from GPU with SUCCESS (drink most tea)\n");
  #endif
#endif
}

// *****************************************
// Run the GPU code for the TF1
__host__ void RunGPU_TF1(
    float *gpu_coeffs,
    int* gpu_paramNo_arr,
    int* gpu_nPoints_arr,

    float* gpu_weights,
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    int* index_gpu,
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif

    // Holds the changes in parameters
    float *vals) {
  // *****************************************

  dim3 block_size;
  dim3 grid_size;

  block_size.x = 1024;
  grid_size.x = (h_n_splines / block_size.x) + 1;

  // Copy the parameter values values to the GPU (vals_gpu), which is h_n_params long
  cudaMemcpyToSymbol(val_gpu, vals, h_n_params*sizeof(float));
  CudaCheckError();

#ifdef DEBUG
  printf("\n***********************\nGPU DEBUGGING ENABLED\n***********************\n");
  printf("block_size.x = %i, grid_size.x = %i \n", block_size.x, grid_size.x);
  printf("RunGPU_TF1 segments\n");
  for (int i = 0; i < h_n_params; i++) {
    printf("val[%i] = %f \n", i, vals[i]);
  }
  printf("nParams = %i, n_splines = %i", h_n_params, h_n_splines);
  printf("\n***********************\nAM NOW CALLING KERNEL\n***********************\n");
#endif

  // Set the cache config to prefer L1 for the kernel
  //cudaFuncSetCacheConfig(EvalOnGPU_SepMany, cudaFuncCachePreferL1);
  EvalOnGPU_TF1<<<grid_size, block_size>>>(
      gpu_coeffs,
      gpu_paramNo_arr,
      gpu_nPoints_arr,

      gpu_weights
      );
  CudaCheckError();

#ifdef DEBUG
  printf("Evaluated TF1 kernel with SUCCESS (drink beer)\n");
#endif

//KS: We can either copy gpu_weight and calculate total weight in reweighting loop, or not copy and calcualate total weight stall at GPU, which means less memory transfer
#ifdef Weight_On_SplineBySpline_Basis
  // Here we have to make a somewhat large GPU->CPU transfer because it's all the splines' response
  cudaMemcpy(cpu_weights, gpu_weights, h_n_splines*sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();

  #ifdef DEBUG
  printf("Copied TF1 GPU weights to CPU with SUCCESS (drink moar beer)\n");
  printf("Released TF1 calculated response from GPU with SUCCESS (drink most beer)\n");
  #endif

//KS: Else calculate Total Weight
#else
  grid_size.x = (h_n_events / block_size.x) + 1;

  #ifdef DEBUG
  printf("\n***********************\nGPU DEBUGGING ENABLED\n***********************\n");
  printf("block_size.x = %i, grid_size.x = %i \n", block_size.x, grid_size.x);
  printf("RunGPU_TotWeight\n");

  printf("nEvents = %i, n_splines = %i, d_n_params", h_n_events, h_n_splines, h_n_params);
  printf("\n***********************\nI AM NOW CALLING KERNEL\n***********************\n");
  #endif
  EvalOnGPU_TotWeight<<<grid_size, block_size>>>(
      index_gpu,
      gpu_weights,
      gpu_total_weights
      );
  #ifdef DEBUG
  CudaCheckError();
  printf("Evaluated kernel with SUCCESS (drink tea)\n");
  #endif
  //KS: Here we have to make a somewhat large GPU->CPU transfer because it is proportional to number of events
  cudaMemcpy(cpu_total_weights, gpu_total_weights, h_n_events*sizeof(float), cudaMemcpyDeviceToHost);
  #ifdef DEBUG
  CudaCheckError();
  printf("Copied GPU total weights to CPU with SUCCESS (drink moar tea)\n");
  printf("Released calculated response from GPU with SUCCESS (drink most tea)\n");
  #endif
#endif
}

// *********************************
// CLEANING
// *********************************

// *********************************
// Clean up the {x},{ybcd} arrays
__host__ void CleanupGPU_SepMany( int *gpu_paramNo_arr,
    // *********************************
    float *gpu_x_array,
    float *gpu_many_array,
#ifndef Weight_On_SplineBySpline_Basis
    float *gpu_total_weights,
    int *index_gpu,
#endif
    float *gpu_weights) {

  cudaFree(gpu_paramNo_arr);

  // free the coefficient arrays
  cudaFree(gpu_x_array);
  cudaFree(gpu_many_array);

  // free weights on the gpu
  cudaFree(gpu_weights);
#ifndef Weight_On_SplineBySpline_Basis
  cudaFree(gpu_total_weights);
  cudaFree(index_gpu);
#endif
  return;
}

// *********************************
// Clean up the TF1 arrays
__host__ void CleanupGPU_TF1(
    // *********************************
    float *gpu_coeffs,
    int *gpu_paramNo_arr,
    int *gpu_nPoints_arr,

#ifndef Weight_On_SplineBySpline_Basis
    int *index_gpu,
    float *gpu_total_weights,
#endif
    float *gpu_weights) {

  cudaFree(gpu_coeffs);
  cudaFree(gpu_paramNo_arr);
  cudaFree(gpu_nPoints_arr);
  cudaFree(gpu_weights);
#ifndef Weight_On_SplineBySpline_Basis
  cudaFree(gpu_total_weights);
  cudaFree(index_gpu);
#endif

  return;
}

