// MaCh3 event-by-event cross-section spline code
// Written by Richard Calland, Asher Kaboth, Clarence Wret, Kamil Skwarczynski
// 
// Contains code to run on CUDA GPUs. Essentially we load up stripped TSpline3 objects to the GPU and do the equivalent of TSpline3->Eval(double) for all events
// Now also supports TF1 evals
// Called from samplePDF/samplePDFND.cpp -> splines/SplineMonolith.cpp -> splines/gpuSplineUtils.cu

//MaCh3 included
#include "manager/gpuUtils.cu"
#include "splines/SplineCommon.h"

// Hard code the number of splines
// Not entirely necessary: only used for val_gpu and segment_gpu being device constants. Could move them to not being device constants
// EM: for OA2022:
#ifdef NSPLINES_ND280
#pragma message("using User Specified N splines")
#define _N_SPLINES_ NSPLINES_ND280
// EM: for OA2024:
#else
#define _N_SPLINES_ 160
#pragma message("using default N splines")
#endif


// KS: Forgive me father, for I have sinned.
#if defined(__CUDA_ARCH__)
  #if __CUDA_ARCH__ >= 1200
    #pragma message("Compiling with CUDA Architecture: 12.x")
  #elif __CUDA_ARCH__ >= 1100
    #pragma message("Compiling with CUDA Architecture: 11.x")
  #elif __CUDA_ARCH__ >= 1000
    #pragma message("Compiling with CUDA Architecture: 10.x")
  #elif __CUDA_ARCH__ >= 900
    #pragma message("Compiling with CUDA Architecture: 9.x")
  #elif __CUDA_ARCH__ >= 800
    #pragma message("Compiling with CUDA Architecture: 8.x")
  #elif __CUDA_ARCH__ >= 750
    #pragma message("Compiling with CUDA Architecture: 7.5")
  #elif __CUDA_ARCH__ >= 730
    #pragma message("Compiling with CUDA Architecture: 7.3")
  #elif __CUDA_ARCH__ >= 720
    #pragma message("Compiling with CUDA Architecture: 7.2")
  #elif __CUDA_ARCH__ >= 710
    #pragma message("Compiling with CUDA Architecture: 7.1")
  #elif __CUDA_ARCH__ >= 700
    #pragma message("Compiling with CUDA Architecture: 7.x")
  #elif __CUDA_ARCH__ >= 650
    #pragma message("Compiling with CUDA Architecture: 6.5")
  #elif __CUDA_ARCH__ >= 600
    #pragma message("Compiling with CUDA Architecture: 6.x")
  #elif __CUDA_ARCH__ >= 530
    #pragma message("Compiling with CUDA Architecture: 5.3")
  #elif __CUDA_ARCH__ >= 520
    #pragma message("Compiling with CUDA Architecture: 5.2")
  #elif __CUDA_ARCH__ >= 510
    #pragma message("Compiling with CUDA Architecture: 5.1")
  #elif __CUDA_ARCH__ >= 500
    #pragma message("Compiling with CUDA Architecture: 5.x")
  #elif __CUDA_ARCH__ >= 400
    #pragma message("Compiling with CUDA Architecture: 4.x")
  #elif __CUDA_ARCH__ >= 300
    #pragma message("Compiling with CUDA Architecture: 3.x")
  #else
    #pragma message("Compiling with CUDA Architecture: < 3.x")
  #endif
#endif

// ******************************************
// CONSTANTS
// ******************************************

// d_NAME declares DEVICE constants (live on GPU)
__device__ __constant__ unsigned int d_n_splines;
__device__ __constant__ short int d_spline_size;
#ifndef Weight_On_SplineBySpline_Basis
__device__ __constant__ int d_n_events;
#endif
/// CW: Constant memory needs to be hard-coded on compile time. Could make this texture memory instead, but don't care enough right now...
__device__ __constant__ float val_gpu[_N_SPLINES_];
__device__ __constant__ short int segment_gpu[_N_SPLINES_];

// h_NAME declares HOST constants (live on CPU)
static short int h_spline_size  = -1;
static int h_n_params     = -1;
#ifndef Weight_On_SplineBySpline_Basis
static int h_n_events     = -1;
#endif

// ******************************************
// TEXTURES
// ******************************************
//KS: Textures are L1 cache variables which are well optimised for fetching. Make texture only for variables you often access but rarely overwrite. There are limits on texture memory so don't use huge arrays
cudaTextureObject_t text_coeff_x = 0;
#ifndef Weight_On_SplineBySpline_Basis
//KS: Map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
cudaTextureObject_t text_nParamPerEvent = 0;
#endif


// *******************************************
//              INITIALISE GPU
// *******************************************

// *******************************************
/// Initialiser when using the x array and combined y,b,c,d array
__host__ void InitGPU_SepMany(
// *******************************************
                          float **gpu_x_array,
                          float **gpu_many_array,
                          float **gpu_weights, 

                          short int** gpu_paramNo_arr,
                          unsigned int** gpu_nKnots_arr,

                 #ifndef Weight_On_SplineBySpline_Basis
                          float **cpu_total_weights, 
                          float **gpu_total_weights, 
                          int n_events,                              
                          unsigned int** gpu_nParamPerEvent,
                  #endif   
                          unsigned int sizeof_array,
                          unsigned int n_splines,
                          int Eve_size) {

  // Allocate chunks of memory to GPU
  cudaMalloc((void **) gpu_paramNo_arr, n_splines*sizeof(short int));
  CudaCheckError();

  cudaMalloc((void **) gpu_nKnots_arr, n_splines*sizeof(unsigned int));
  CudaCheckError();

  cudaMalloc((void **) gpu_x_array, Eve_size*sizeof(float));
  CudaCheckError();

  cudaMalloc((void **) gpu_many_array, _nCoeff_*sizeof_array*sizeof(float));
  CudaCheckError();

  // Allocate memory for the array of weights to be returned to CPU
  cudaMalloc((void **) gpu_weights, n_splines*sizeof(float));
  CudaCheckError();
#ifndef Weight_On_SplineBySpline_Basis
  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) cpu_total_weights, n_events*sizeof(float));
  CudaCheckError();

  //KS: Allocate memory for the array of total weights to be returned to CPU
  cudaMalloc((void **) gpu_total_weights, n_events*sizeof(float));
  CudaCheckError();
  
  //KS: Allocate memory for the map keeping track how many splines each parameter has
  cudaMalloc((void **) gpu_nParamPerEvent, 2*n_events*sizeof(unsigned int));
  CudaCheckError();
  
#endif
  
  // Print allocation info to user
  printf("  Allocated %i entries for paramNo and nKnots arrays, size = %f MB\n", n_splines, double(sizeof(short int)*n_splines+sizeof(unsigned int)*n_splines)/1.E6);
  printf("  Allocated %i entries for x coeff arrays, size = %f MB\n", Eve_size, double(sizeof(float)*Eve_size)/1.E6);
  printf("  Allocated %i entries for {ybcd} coeff arrays, size = %f MB\n", _nCoeff_*sizeof_array, double(sizeof(float)*_nCoeff_*sizeof_array)/1.E6);

  //KS: Ask CUDA about memory usage
  checkGpuMem();
  PrintNdevices();
}

// *******************************************
/// Initialiser when using the x array and combined y,b,c,d array
__host__ void InitGPU_TF1(
// *******************************************
                          float **gpu_coeffs,
                          short int** gpu_paramNo_arr,
                          short int** gpu_nPoints_arr,
                          float **gpu_weights, 
                             
                    #ifndef Weight_On_SplineBySpline_Basis
                          float **cpu_total_weights, 
                          float **gpu_total_weights, 
                          int n_events,
                              
                          unsigned int** gpu_nParamPerEvent,
                    #endif  
                          unsigned int n_splines) {

  // Holds the parameter number
  cudaMalloc((void **) gpu_paramNo_arr, n_splines*sizeof(short int));
  CudaCheckError();

  // Holds the number of points
  cudaMalloc((void **) gpu_nPoints_arr, n_splines*sizeof(short int));
  CudaCheckError();

  // Holds the coefficients (5th order polynomial and constant term == 1) -> 5
  cudaMalloc((void **) gpu_coeffs, 5*n_splines*sizeof(float));
  CudaCheckError();

  // Allocate memory for the array of weights to be returned to CPU
  cudaMalloc((void **) gpu_weights, n_splines*sizeof(float));
  CudaCheckError();

#ifndef Weight_On_SplineBySpline_Basis
  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) cpu_total_weights, n_events*sizeof(float));
  CudaCheckError();
  
  //KS: Allocate memory for the array of total weights to be returned to CPU
  cudaMalloc((void **) gpu_total_weights, n_events*sizeof(float));
  CudaCheckError();
  
  //KS: Allocate memory for the map keeping track how many splines each parameter has
  cudaMalloc((void **) gpu_nParamPerEvent, 2*n_events*sizeof(unsigned int));
  CudaCheckError();
#endif
  
  // Print allocation info to user
  printf("  Allocated %i entries for paramNo and nPoints arrays, size = %f MB\n", n_splines, double(2.0*sizeof(int)*n_splines)/1.E6);
  printf("  Allocated %i entries for coefficient arrays, size = %f MB\n", 5*n_splines, double(sizeof(float)*5*n_splines)/1.E6);

  //KS: Ask CUDA about memory usage
  checkGpuMem();
  PrintNdevices();
}


// *******************************************
// Allocate memory for spline segments
__host__ void InitGPU_Segments(short int **segment) {
// *******************************************

  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) segment, _N_SPLINES_*sizeof(short int));
  CudaCheckError();
}

// *******************************************
// Allocate memory for spline segments
__host__ void InitGPU_Vals(float **vals) {
// *******************************************

  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) vals, _N_SPLINES_*sizeof(float));
  CudaCheckError();
}


// ******************************************************
//                START COPY TO GPU
// ******************************************************

// ******************************************************
/// Copy to GPU for x array and separate ybcd array
__host__ void CopyToGPU_SepMany(
// ******************************************************
                            short int *gpu_paramNo_arr,
                            unsigned int *gpu_nKnots_arr,
                            float *gpu_x_array,
                            float *gpu_many_array,

                            std::vector<short int> paramNo_arr,
                            std::vector<unsigned int> nKnots_arr,
                            std::vector<float> cpu_x_array,
                            std::vector<float> cpu_many_array,

                    #ifndef Weight_On_SplineBySpline_Basis
                            int n_events,
                            std::vector<unsigned int> cpu_nParamPerEvent,
                            unsigned int *gpu_nParamPerEvent,
                    #endif
                            int n_params, 
                            unsigned int n_splines,
                            short int spline_size,
                            unsigned int sizeof_array) {
  if (n_params != _N_SPLINES_) {
    printf("Number of splines not equal to %i, GPU code for event-by-event splines will fail\n", _N_SPLINES_);
    printf("n_params = %i\n", n_params);
    printf("%s : %i\n", __FILE__, __LINE__);
    exit(-1);
  }

  // Write to the global statics (h_* denotes host stored variable)
  h_n_params = n_params;
  h_spline_size = spline_size;
#ifndef Weight_On_SplineBySpline_Basis
  h_n_events    = n_events;
#endif
  // Copy the constants
  // Total number of valid splines for all loaded events
  cudaMemcpyToSymbol(d_n_splines,   &n_splines,   sizeof(n_splines));
  CudaCheckError();
  // Total spline size per spline; i.e. just the number of points or knots in the spline
  cudaMemcpyToSymbol(d_spline_size, &h_spline_size, sizeof(h_spline_size));
  CudaCheckError();
#ifndef Weight_On_SplineBySpline_Basis
  // Number of events
  cudaMemcpyToSymbol(d_n_events, &h_n_events, sizeof(h_n_events));
  CudaCheckError();
#endif
  // Copy the coefficient arrays to the GPU; this only happens once per entire Markov Chain so is OK to do multiple extensive memory copies
  cudaMemcpy(gpu_many_array, cpu_many_array.data(), sizeof(float)*sizeof_array*_nCoeff_, cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_x_array, cpu_x_array.data(), sizeof(float)*spline_size*n_params, cudaMemcpyHostToDevice);
  CudaCheckError();

  //KS: Bind our texture with the GPU variable
  //KS: Tried also moving gpu_many_array to texture memory it only worked with restricted number of MC runs, most likely hit texture memory limit :(
  struct cudaResourceDesc resDesc_coeff_x;
  memset(&resDesc_coeff_x, 0, sizeof(resDesc_coeff_x));
  resDesc_coeff_x.resType = cudaResourceTypeLinear;
  resDesc_coeff_x.res.linear.devPtr = gpu_x_array;
  resDesc_coeff_x.res.linear.desc = cudaCreateChannelDesc<float>();
  resDesc_coeff_x.res.linear.sizeInBytes = sizeof(float)*spline_size*n_params;

  // Specify texture object parameters
  struct cudaTextureDesc texDesc_coeff_x;
  memset(&texDesc_coeff_x, 0, sizeof(texDesc_coeff_x));
  texDesc_coeff_x.readMode = cudaReadModeElementType;

  // Create texture object
  cudaCreateTextureObject(&text_coeff_x, &resDesc_coeff_x, &texDesc_coeff_x, NULL);
  CudaCheckError();

  // Also copy the parameter number for each spline onto the GPU; i.e. what spline parameter are we calculating right now
  cudaMemcpy(gpu_paramNo_arr, paramNo_arr.data(), n_splines*sizeof(short int), cudaMemcpyHostToDevice);
  CudaCheckError();

  // Also copy the knot map for each spline onto the GPU;
  cudaMemcpy(gpu_nKnots_arr, nKnots_arr.data(), n_splines*sizeof(unsigned int), cudaMemcpyHostToDevice);
  CudaCheckError();

  #ifndef Weight_On_SplineBySpline_Basis
  //KS: Keep track how much splines each event has  
  cudaMemcpy(gpu_nParamPerEvent, cpu_nParamPerEvent.data(), 2*n_events*sizeof(unsigned int), cudaMemcpyHostToDevice);
  CudaCheckError();
  
  //KS: Bind our texture with the GPU variable
  // create a resource descriptor based on device pointers
  struct cudaResourceDesc resDesc_nParamPerEvent;
  memset(&resDesc_nParamPerEvent, 0, sizeof(resDesc_nParamPerEvent));
  resDesc_nParamPerEvent.resType = cudaResourceTypeLinear;
  resDesc_nParamPerEvent.res.linear.devPtr = gpu_nParamPerEvent;
  resDesc_nParamPerEvent.res.linear.desc = cudaCreateChannelDesc<unsigned int>();
  resDesc_nParamPerEvent.res.linear.sizeInBytes = 2*n_events*sizeof(unsigned int);

  // Specify texture object parameters
  struct cudaTextureDesc texDesc_nParamPerEvent;
  memset(&texDesc_nParamPerEvent, 0, sizeof(texDesc_nParamPerEvent));
  texDesc_nParamPerEvent.readMode = cudaReadModeElementType;

  //Finally create texture object
  cudaCreateTextureObject(&text_nParamPerEvent, &resDesc_nParamPerEvent, &texDesc_nParamPerEvent, NULL);
  CudaCheckError();
  #endif
}

// ******************************************************
/// Copy to GPU for x array and separate ybcd array
__host__ void CopyToGPU_TF1(
// ******************************************************
                            float *gpu_coeffs,
                            short int *gpu_paramNo_arr,
                            short int *gpu_nPoints_arr,

                            std::vector<float> cpu_coeffs,
                            std::vector<short int> paramNo_arr,
                            std::vector<short int> nPoints_arr,

                  #ifndef Weight_On_SplineBySpline_Basis
                            int n_events,
                            std::vector<unsigned int> cpu_nParamPerEvent,
                            unsigned int *gpu_nParamPerEvent,
                  #endif
                            int n_params,
                            unsigned int n_splines,
                            short int _max_knots) {

  if (n_params != _N_SPLINES_) {
    printf("Number of splines not equal to %i, GPU code for event-by-event splines will fail\n", _N_SPLINES_);
    printf("n_params = %i\n", n_params);
    printf("%s : %i\n", __FILE__, __LINE__);
    exit(-1);
  }

  // Write to the global statics (h_* denotes host stored variable)
  h_n_params = n_params;
  h_spline_size = _max_knots;
#ifndef Weight_On_SplineBySpline_Basis
  h_n_events    = n_events;
#endif
  // Copy the constants
  // Total number of valid splines for all loaded events
  cudaMemcpyToSymbol(d_n_splines,   &n_splines,   sizeof(n_splines));
  CudaCheckError();
  // Total spline size per spline; i.e. just the number of points or knots in the spline
  cudaMemcpyToSymbol(d_spline_size, &h_spline_size, sizeof(h_spline_size));
  CudaCheckError();

#ifndef Weight_On_SplineBySpline_Basis
  // Number of events
  cudaMemcpyToSymbol(d_n_events, &h_n_events, sizeof(h_n_events));
  CudaCheckError();
#endif
  // Move the coefficients
  cudaMemcpy(gpu_coeffs, cpu_coeffs.data(), n_splines*5*sizeof(float), cudaMemcpyHostToDevice);
  CudaCheckError();

  // Also copy the parameter number for each spline onto the GPU; i.e. what spline parameter are we calculating right now
  cudaMemcpy(gpu_paramNo_arr, paramNo_arr.data(), n_splines*sizeof(short int), cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_nPoints_arr, nPoints_arr.data(), n_splines*sizeof(short int), cudaMemcpyHostToDevice);
  CudaCheckError();
  
  #ifndef Weight_On_SplineBySpline_Basis
  //KS: Keep track how much splines each event has  
  cudaMemcpy(gpu_nParamPerEvent, cpu_nParamPerEvent.data(), 2*n_events*sizeof(unsigned int), cudaMemcpyHostToDevice);
  CudaCheckError();

  //KS: Bind our texture with the GPU variable
  // create a resource descriptor based on device pointers
  struct cudaResourceDesc resDesc_nParamPerEvent;
  memset(&resDesc_nParamPerEvent, 0, sizeof(resDesc_nParamPerEvent));
  resDesc_nParamPerEvent.resType = cudaResourceTypeLinear;
  resDesc_nParamPerEvent.res.linear.devPtr = gpu_nParamPerEvent;
  resDesc_nParamPerEvent.res.linear.desc = cudaCreateChannelDesc<float>();
  resDesc_nParamPerEvent.res.linear.sizeInBytes = 2*n_events*sizeof(unsigned int);

  // Specify texture object parameters
  struct cudaTextureDesc texDesc_nParamPerEvent;
  memset(&texDesc_nParamPerEvent, 0, sizeof(texDesc_nParamPerEvent));
  texDesc_nParamPerEvent.readMode = cudaReadModeElementType;

  //Lastly create texture object
  cudaCreateTextureObject(&text_nParamPerEvent, &resDesc_nParamPerEvent, &texDesc_nParamPerEvent, NULL);
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
// Should be most efficient at cache hitting and memory coalescence
// But using spline segments rather than the parameter value: avoids doing binary search on GPU
__global__ void EvalOnGPU_SepMany(
    const short int* __restrict__ gpu_paramNo_arr,
    const unsigned int* __restrict__ gpu_nKnots_arr,
    const float* __restrict__ gpu_coeff_many,
    float *gpu_weights,
    const cudaTextureObject_t __restrict__ text_coeff_x) {
//*********************************************************

  // points per spline is the offset to skip in the index to move between splines
  const unsigned int splineNum = (blockIdx.x * blockDim.x + threadIdx.x);

  // Note, the new arrays are arranged as:
  //       gpu_paramNo_arr has length = spln_counter (keeps track of which parameter we're using on this thread)
  //       gpu_nKnots_arrhas length = spln_counter (keeps track where current spline starts)
  //       text_coeff_x has length = n_params * spline_size
  //       gpu_coeff_many has length = nKnots * 4
  //       ...
  //       gpu_weights has length = spln_counter * spline_size

  // this is the stopping condition!
  if (splineNum < d_n_splines) {
    // This is the segment we want for this parameter variation
    // for this particular splineNum; 0 = MACCQE, 1 = pFC, 2 = EBC, etc

    //CW: Which Parameter we are accesing
    const short int Param = gpu_paramNo_arr[splineNum];

    //CW: Avoids doing costly binary search on GPU
    const short int segment = segment_gpu[Param];

    //KS: Segment for coeff_x is simply parameter*max knots + segment as each parmeters has the same spacing
    const short int segment_X = Param*d_spline_size+segment;

    //KS: Find knot position in out monolitical structure
    const unsigned int CurrentKnotPos = gpu_nKnots_arr[splineNum]*_nCoeff_+segment*_nCoeff_;

    // We've read the segment straight from CPU and is saved in segment_gpu
    // polynomial parameters from the monolithic splineMonolith
    const float fY = gpu_coeff_many[CurrentKnotPos];
    const float fB = gpu_coeff_many[CurrentKnotPos+1];
    const float fC = gpu_coeff_many[CurrentKnotPos+2];
    const float fD = gpu_coeff_many[CurrentKnotPos+3];
    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float dx = val_gpu[Param] - tex1Dfetch<float>(text_coeff_x, segment_X);

    //CW: Wooow, let's use some fancy intrinsics and pull down the processing time by <1% from normal multiplication! HURRAY
    gpu_weights[splineNum] = fmaf(dx, fmaf(dx, fmaf(dx, fD, fC), fB), fY);
    // Or for the more "easy to read" version:
    //gpu_weights[splineNum] = (fY+dx*(fB+dx*(fC+dx*fD)));

#ifdef DEBUG
  printf("splineNum = %i/%i, paramNo = %i, variation = %f, segment = %i, fX = %f, fX+1 = %f, dx = %f, d_n_splines = %i, d_spline_size = %i, weight = %f \n", splineNum, d_n_splines, gpu_paramNo_arr[splineNum], val_gpu[Param], segment, tex1Dfetch<float>(text_coeff_x, segment_X), tex1Dfetch<float>(text_coeff_x, segment_X+1), dx, d_n_splines, d_spline_size, gpu_weights[splineNum]);
#endif
  }
}

//*********************************************************
/// Evaluate the TF1 on the GPU Using 5th order polynomial
__global__ void EvalOnGPU_TF1( 
    const float* __restrict__ gpu_coeffs,
    const short int* __restrict__ gpu_paramNo_arr,
    const short int* __restrict__ gpu_nPoints_arr,
    float *gpu_weights) {
//*********************************************************

  // points per spline is the offset to skip in the index to move between splines
  const unsigned int splineNum = (blockIdx.x * blockDim.x + threadIdx.x);

  // Note, the new arrays are arranged as:
  //       gpu_paramNo_arr has length = spln_counter (keeps track of which parameter we're using on this thread)
  //       gpu_coeff_x has length = n_params * spline_size
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
/// KS: Evaluate the total spline event weight on the GPU, as in most cases GPU is faster, even more this significant reduce memory transfer from GPU to CPU
__global__ void EvalOnGPU_TotWeight(
   const float* __restrict__ gpu_weights,
   float *gpu_total_weights,
  const cudaTextureObject_t __restrict__ text_nParamPerEvent) {
//*********************************************************
  const unsigned int EventNum = (blockIdx.x * blockDim.x + threadIdx.x);
  //KS: Accessing shared memory is much much faster than global memory hence we use shared memory for calculation and then write to global memory
  __shared__ float shared_total_weights[_BlockSize_];
  if(EventNum < d_n_events) //stopping condition
  {
    shared_total_weights[threadIdx.x] = 1.f;
    const unsigned int EventOffset = 2*EventNum;
    for (unsigned int id = 0; id < tex1Dfetch<unsigned int>(text_nParamPerEvent, EventOffset); ++id)
    {
      shared_total_weights[threadIdx.x] *= gpu_weights[tex1Dfetch<unsigned int>(text_nParamPerEvent, EventOffset+1) + id];
      #ifdef DEBUG
      printf("Event = %i, Spline_Num = %i, gpu_weights = %f \n",
              EventNum, tex1Dfetch<unsigned int>(text_nParamPerEvent, 2*EventNum+1) + id, gpu_weights[tex1Dfetch<unsigned int>(text_nParamPerEvent, 2*EventNum+1) + id];
      #endif
    }
    gpu_total_weights[EventNum] = shared_total_weights[threadIdx.x];
  }
}
#endif

// *****************************************
// Run the GPU code for the separate many arrays
// As in separate {x}, {y,b,c,d} arrays
// Pass the segment and the parameter values
// (binary search already performed in samplePDFND::FindSplineSegment()
__host__ void RunGPU_SepMany(
    const short int* gpu_paramNo_arr,
    const unsigned int* gpu_nKnots_arr,

    const float *gpu_coeff_many,

    float* gpu_weights, 
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif
    // Holds the changes in parameters
    float *vals,
    // Holds the segments for parameters
    short int *segment,
    const unsigned int h_n_splines) {
// *****************************************

  dim3 block_size;
  dim3 grid_size;

  block_size.x = _BlockSize_;
  grid_size.x = (h_n_splines / block_size.x) + 1;

  // Copy the segment values to the GPU (segment_gpu), which is h_n_params long
  cudaMemcpyToSymbol(segment_gpu, segment, h_n_params*sizeof(short int));
  CudaCheckError();

  // Copy the parameter values values to the GPU (vals_gpu), which is h_n_params long
  cudaMemcpyToSymbol(val_gpu, vals, h_n_params*sizeof(float));
  CudaCheckError();

#ifdef DEBUG
  printf("\n***********************\nGPU DEBUGGING ENABLED\n***********************\n");
  printf("block_size.x = %i, grid_size.x = %i \n", block_size.x, grid_size.x);
  printf("RunGPU_SepMany segments\n");
  for (int i = 0; i < h_n_params; i++) {
    printf("val[%i] = %f in segment %i\n", i, vals[i], segment[i]);
  }
  printf("nParams = %i, n_splines = %i", h_n_params, h_n_splines);
  printf("\n***********************\nAM NOW CALLING KERNEL\n***********************\n");
#endif

  // Set the cache config to prefer L1 for the kernel
  //cudaFuncSetCacheConfig(EvalOnGPU_SepMany, cudaFuncCachePreferL1);
  EvalOnGPU_SepMany<<<grid_size, block_size>>>(
      gpu_paramNo_arr,
      gpu_nKnots_arr,

      gpu_coeff_many,

      gpu_weights,
      text_coeff_x
      );
  CudaCheckError();

#ifdef DEBUG
  printf("Evaluated kernel with SUCCESS (drink beer)\n");
#endif

//KS: We can either copy gpu_weight and calculate total weight in reweighting loop, or not copy and calculate total weight stall at GPU, which means less memory transfer
#ifdef Weight_On_SplineBySpline_Basis
  // Here we have to make a somewhat large GPU->CPU transfer because it's all the splines' response
  cudaMemcpy(cpu_weights, gpu_weights, h_n_splines*sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();
  #ifdef DEBUG
  printf("Copied GPU weights to CPU with SUCCESS (drink more beer)\n");
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
      gpu_weights,
      gpu_total_weights,
      text_nParamPerEvent
      );
  #ifdef DEBUG
  CudaCheckError();
  printf("Evaluated kernel with SUCCESS (drink tea)\n");
  #endif
  //KS: Here we have to make a somewhat large GPU->CPU transfer because it is proportional to number of events
  //KS: Normally code wait for memory transfer to finish before moving further cudaMemcpyAsync means we will continue to execute code and in a meantime keep copying stuff.
  cudaMemcpyAsync(cpu_total_weights, gpu_total_weights, h_n_events * sizeof(float), cudaMemcpyDeviceToHost, 0);

  #ifdef DEBUG
  CudaCheckError();
  printf("Copied GPU total weights to CPU with SUCCESS (drink moar tea)\n");
  printf("Released calculated response from GPU with SUCCESS (drink most tea)\n");
  #endif
#endif
}

// *****************************************
/// Run the GPU code for the TF1
__host__ void RunGPU_TF1(
    const float *gpu_coeffs,
    const short int* gpu_paramNo_arr,
    const short int* gpu_nPoints_arr,

    float* gpu_weights, 
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif

  // Holds the changes in parameters
    float *vals,
    const unsigned int h_n_splines) {
// *****************************************

  dim3 block_size;
  dim3 grid_size;

  block_size.x = _BlockSize_;
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
  //cudaFuncSetCacheConfig(EvalOnGPU_TF1, cudaFuncCachePreferL1);
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
  
//KS: We can either copy gpu_weight and calculate total weight in reweighting loop, or not copy and calculate total weight stall at GPU, which means less memory transfer
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
      gpu_weights,
      gpu_total_weights,
      text_nParamPerEvent
      );
  #ifdef DEBUG
  CudaCheckError();
  printf("Evaluated kernel with SUCCESS (drink tea)\n");
  #endif
  //KS: Here we have to make a somewhat large GPU->CPU transfer because it is proportional to number of events
  //KS: In the future it might be worth to calculate only weight for events which have splines, this should reduce memory transfer
  cudaMemcpy(cpu_total_weights, gpu_total_weights, h_n_events*sizeof(float), cudaMemcpyDeviceToHost);
  #ifdef DEBUG
  CudaCheckError();
  printf("Copied GPU total weights to CPU with SUCCESS (drink moar tea)\n");
  printf("Released calculated response from GPU with SUCCESS (drink most tea)\n");
  #endif
#endif
}

// *****************************************
/// Make sure all Cuda threads finished execution
__host__ void SynchroniseSplines() {
  cudaDeviceSynchronize();
}

// *********************************
// CLEANING
// *********************************

// *********************************
/// Clean up the {x},{ybcd} arrays
__host__ void CleanupGPU_SepMany( 
    short int *gpu_paramNo_arr,
    unsigned int *gpu_nKnots_arr,

    float *gpu_x_array, 
    float *gpu_many_array, 
#ifndef Weight_On_SplineBySpline_Basis
    float *gpu_total_weights,
    unsigned int *gpu_nParamPerEvent,
    float *cpu_total_weights,
#endif
    float *gpu_weights) {
// *********************************
  cudaFree(gpu_paramNo_arr);
  cudaFree(gpu_nKnots_arr);

  // free the coefficient arrays
  cudaDestroyTextureObject(text_coeff_x);
  cudaFree(gpu_x_array);
  cudaFree(gpu_many_array);

  // free weights on the gpu
  cudaFree(gpu_weights);
#ifndef Weight_On_SplineBySpline_Basis
  cudaFree(gpu_total_weights);
  //KS: Before removing variable let's destroy texture
  cudaDestroyTextureObject(text_nParamPerEvent);
  cudaFree(gpu_nParamPerEvent);
  cudaFreeHost(cpu_total_weights);
  cpu_total_weights = nullptr;
#endif
  return;
}

// *******************************************
/// Clean up pinned variables at CPU
__host__ void CleanupGPU_Segments(short int *segment, float *vals) {
// *******************************************
    cudaFreeHost(segment);
    cudaFreeHost(vals);

    segment = nullptr;
    vals = nullptr;

    return;
}

// *********************************
/// Clean up the TF1 arrays
__host__ void CleanupGPU_TF1(
    float *gpu_coeffs,
    short int *gpu_paramNo_arr,
    short int *gpu_nPoints_arr,
    
#ifndef Weight_On_SplineBySpline_Basis
    float *gpu_total_weights,
    float *cpu_total_weights,
#endif
    float *gpu_weights) {
// *********************************
  cudaFree(gpu_coeffs);
  cudaFree(gpu_paramNo_arr);
  cudaFree(gpu_nPoints_arr);
  cudaFree(gpu_weights);
#ifndef Weight_On_SplineBySpline_Basis
  cudaFree(gpu_total_weights);
  cudaFreeHost(cpu_total_weights);
#endif
  
  return;
}
