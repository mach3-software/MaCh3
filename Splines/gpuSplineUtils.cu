//MaCh3 included
#include "Splines/gpuSplineUtils.cuh"

/// Hard code the number of splines
/// Not entirely necessary: only used for val_gpu and segment_gpu being device constants. Could move them to not being device constants
#define _N_SPLINES_ NSplines_GPU

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
/// Number of splines living on GPU
__device__ __constant__ unsigned int d_n_splines;
/// Number of tf1 living on GPU
__device__ __constant__ unsigned int d_n_TF1;
/// Size of splines living on GPU
__device__ __constant__ short int d_spline_size;
#ifndef Weight_On_SplineBySpline_Basis
/// Number of events living on GPU
__device__ __constant__ int d_n_events;
#endif
/// CW: Constant memory needs to be hard-coded on compile time. Could make this texture memory instead, but don't care enough right now...
__device__ __constant__ float val_gpu[_N_SPLINES_];
__device__ __constant__ short int segment_gpu[_N_SPLINES_];


// *****************************************
// Make sure all Cuda threads finished execution
__host__ void SynchroniseSplines() {
  cudaDeviceSynchronize();
}

// *******************************************
//              INITIALISE GPU
// *******************************************

SMonolithGPU::SMonolithGPU(){
  h_n_params     = -1;
  /// Number of events living on CPU
  h_n_events = -1;

  gpu_weights = nullptr;
  gpu_total_weights = nullptr;
  gpu_nParamPerEvent = nullptr;
  gpu_nPoints_arr = nullptr;
  gpu_paramNo_arr = nullptr;
  gpu_nKnots_arr = nullptr;
  gpu_coeff_x = nullptr;
  gpu_coeff_many = nullptr;
  gpu_coeff_TF1_many = nullptr;
  gpu_paramNo_TF1_arr = nullptr;
  gpu_nParamPerEvent_TF1 = nullptr;
  gpu_weights_tf1 = nullptr;
}

SMonolithGPU::~SMonolithGPU(){

}

// *******************************************
// Initialiser when using the x array and combined y,b,c,d array
__host__ void SMonolithGPU::InitGPU_SplineMonolith(
                 #ifndef Weight_On_SplineBySpline_Basis
                          float **cpu_total_weights, 
                          int n_events,                              
                  #endif   
                          unsigned int total_nknots,
                          unsigned int n_splines,
                          unsigned int n_tf1,
                          int Eve_size) {
// *******************************************
  // Allocate chunks of memory to GPU
  cudaMalloc((void **) &gpu_paramNo_arr, n_splines*sizeof(short int));
  CudaCheckError();

  cudaMalloc((void **) &gpu_nKnots_arr, n_splines*sizeof(unsigned int));
  CudaCheckError();

  cudaMalloc((void **) &gpu_coeff_x, Eve_size*sizeof(float));
  CudaCheckError();

  cudaMalloc((void **) &gpu_coeff_many, _nCoeff_*total_nknots*sizeof(float));
  CudaCheckError();

  // Allocate memory for the array of weights to be returned to CPU
  cudaMalloc((void **) &gpu_weights, n_splines*sizeof(float));
  CudaCheckError();

  // Now TF1 specific
  cudaMalloc((void **) &gpu_coeff_TF1_many, _nTF1Coeff_*n_tf1*sizeof(float));
  CudaCheckError();

  cudaMalloc((void **) &gpu_weights_tf1, n_tf1*sizeof(float));
  CudaCheckError();

  cudaMalloc((void **) &gpu_paramNo_TF1_arr, n_tf1*sizeof(short int));
  CudaCheckError();

#ifndef Weight_On_SplineBySpline_Basis
  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) cpu_total_weights, n_events*sizeof(float));
  CudaCheckError();

  //KS: Allocate memory for the array of total weights to be returned to CPU
  cudaMalloc((void **) &gpu_total_weights, n_events*sizeof(float));
  CudaCheckError();
  
  //KS: Allocate memory for the map keeping track how many splines each parameter has
  cudaMalloc((void **) &gpu_nParamPerEvent, 2*n_events*sizeof(unsigned int));
  CudaCheckError();
  
  //KS: Allocate memory for the map keeping track how many TF1 each parameter has
  cudaMalloc((void **) &gpu_nParamPerEvent_TF1, 2*n_events*sizeof(unsigned int));
  CudaCheckError();
#endif
  
  // Print allocation info to user
  printf("Allocated %i entries for paramNo and nKnots arrays, size = %f MB\n", n_splines, double(sizeof(short int)*n_splines+sizeof(unsigned int)*n_splines)/1.E6);
  printf("Allocated %i entries for x coeff arrays, size = %f MB\n", Eve_size, double(sizeof(float)*Eve_size)/1.E6);
  printf("Allocated %i entries for {ybcd} coeff arrays, size = %f MB\n", _nCoeff_*total_nknots, double(sizeof(float)*_nCoeff_*total_nknots)/1.E6);
  printf("Allocated %i entries for TF1 coefficient arrays, size = %f MB\n", _nTF1Coeff_*n_tf1, double(sizeof(float)*_nTF1Coeff_*n_tf1)/1.E6);

  //KS: Ask CUDA about memory usage
  checkGpuMem();
  PrintNdevices();
}

// *******************************************
// Allocate memory for spline segments
__host__ void SMonolithGPU::InitGPU_Segments(short int **segment) {
// *******************************************
  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) segment, _N_SPLINES_*sizeof(short int));
  CudaCheckError();
}

// *******************************************
// Allocate memory for spline segments
__host__ void SMonolithGPU::InitGPU_Vals(float **vals) {
// *******************************************

  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) vals, _N_SPLINES_*sizeof(float));
  CudaCheckError();
}


// ******************************************************
//                START COPY TO GPU
// ******************************************************

// ******************************************************
// Copy to GPU for x array and separate ybcd array
__host__ void SMonolithGPU::CopyToGPU_SplineMonolith(
                            SplineMonoStruct* cpu_spline_handler,

                            // TFI related now
                            std::vector<float> cpu_many_array_TF1,
                            std::vector<short int> cpu_paramNo_arr_TF1,
                    #ifndef Weight_On_SplineBySpline_Basis
                            int n_events,
                            std::vector<unsigned int> cpu_nParamPerEvent,
                            // TFI related now
                            std::vector<unsigned int> cpu_nParamPerEvent_TF1,
                    #endif
                            int n_params, 
                            unsigned int n_splines,
                            short int spline_size,
                            unsigned int total_nknots,
                            unsigned int n_tf1) {
// ******************************************************
  if (n_params != _N_SPLINES_) {
    printf("Number of splines not equal to %i, GPU code for event-by-event splines will fail\n", _N_SPLINES_);
    printf("n_params = %i\n", n_params);
    printf("%s : %i\n", __FILE__, __LINE__);
    exit(-1);
  }

  // Write to the global statics (h_* denotes host stored variable)
  h_n_params = n_params;
#ifndef Weight_On_SplineBySpline_Basis
  h_n_events    = n_events;
#endif
  // Copy the constants
  // Total number of valid splines for all loaded events
  cudaMemcpyToSymbol(d_n_splines, &n_splines, sizeof(n_splines));
  CudaCheckError();

  // Total number of valid TF1 for all loaded events
  cudaMemcpyToSymbol(d_n_TF1,   &n_tf1, sizeof(n_tf1));
  CudaCheckError();

  // Total spline size per spline; i.e. just the number of points or knots in the spline
  cudaMemcpyToSymbol(d_spline_size, &spline_size, sizeof(spline_size));
  CudaCheckError();
#ifndef Weight_On_SplineBySpline_Basis
  // Number of events
  cudaMemcpyToSymbol(d_n_events, &h_n_events, sizeof(h_n_events));
  CudaCheckError();
#endif
  // Copy the coefficient arrays to the GPU; this only happens once per entire Markov Chain so is OK to do multiple extensive memory copies
  cudaMemcpy(gpu_coeff_many, cpu_spline_handler->coeff_many.data(), sizeof(float)*total_nknots*_nCoeff_, cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_coeff_x, cpu_spline_handler->coeff_x.data(), sizeof(float)*spline_size*n_params, cudaMemcpyHostToDevice);
  CudaCheckError();

  //KS: Bind our texture with the GPU variable
  //KS: Tried also moving gpu_many_array to texture memory it only worked with restricted number of MC runs, most likely hit texture memory limit :(
  struct cudaResourceDesc resDesc_coeff_x;
  memset(&resDesc_coeff_x, 0, sizeof(resDesc_coeff_x));
  resDesc_coeff_x.resType = cudaResourceTypeLinear;
  resDesc_coeff_x.res.linear.devPtr = gpu_coeff_x;
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
  cudaMemcpy(gpu_paramNo_arr, cpu_spline_handler->paramNo_arr.data(), n_splines*sizeof(short int), cudaMemcpyHostToDevice);
  CudaCheckError();

  // Also copy the knot map for each spline onto the GPU;
  cudaMemcpy(gpu_nKnots_arr, cpu_spline_handler->nKnots_arr.data(), n_splines*sizeof(unsigned int), cudaMemcpyHostToDevice);
  CudaCheckError();

  //Now TF1
  // Copy the coefficient arrays to the GPU; this only happens once per entire Markov Chain so is OK to do multiple extensive memory copies
  cudaMemcpy(gpu_coeff_TF1_many, cpu_many_array_TF1.data(), sizeof(float)*n_tf1*_nTF1Coeff_, cudaMemcpyHostToDevice);
  CudaCheckError();

  // Also copy the parameter number for each TF1 onto the GPU; i.e. what TF1 parameter are we calculating right now
  cudaMemcpy(gpu_paramNo_TF1_arr, cpu_paramNo_arr_TF1.data(), n_tf1*sizeof(short int), cudaMemcpyHostToDevice);
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

  // Now TF1
  cudaMemcpy(gpu_nParamPerEvent_TF1, cpu_nParamPerEvent_TF1.data(), 2*n_events*sizeof(unsigned int), cudaMemcpyHostToDevice);
  CudaCheckError();

  //KS: Bind our texture with the GPU variable
  // create a resource descriptor based on device pointers
  struct cudaResourceDesc resDesc_nParamPerEvent_tf1;
  memset(&resDesc_nParamPerEvent_tf1, 0, sizeof(resDesc_nParamPerEvent_tf1));
  resDesc_nParamPerEvent_tf1.resType = cudaResourceTypeLinear;
  resDesc_nParamPerEvent_tf1.res.linear.devPtr = gpu_nParamPerEvent_TF1;
  resDesc_nParamPerEvent_tf1.res.linear.desc = cudaCreateChannelDesc<unsigned int>();
  resDesc_nParamPerEvent_tf1.res.linear.sizeInBytes = 2*n_events*sizeof(unsigned int);

  // Specify texture object parameters
  struct cudaTextureDesc texDesc_nParamPerEvent_tf1;
  memset(&texDesc_nParamPerEvent_tf1, 0, sizeof(texDesc_nParamPerEvent_tf1));
  texDesc_nParamPerEvent_tf1.readMode = cudaReadModeElementType;

  //Finally create texture object
  cudaCreateTextureObject(&text_nParamPerEvent_TF1, &resDesc_nParamPerEvent_tf1, &texDesc_nParamPerEvent_tf1, NULL);
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
// Evaluate the spline on the GPU Using one {y,b,c,d} array and one {x} array
// Should be most efficient at cache hitting and memory coalescence
// But using spline segments rather than the parameter value: avoids doing binary search on GPU
__global__ void EvalOnGPU_Splines(
  const short int* __restrict__ gpu_paramNo_arr,
  const unsigned int* __restrict__ gpu_nKnots_arr,
  const float* __restrict__ gpu_coeff_many,
  float* __restrict__ gpu_weights,
  const cudaTextureObject_t __restrict__ text_coeff_x) {
//*********************************************************
  // points per spline is the offset to skip in the index to move between splines
  const unsigned int splineNum = (blockIdx.x * blockDim.x + threadIdx.x);

  // this is the stopping condition!
  if (splineNum < d_n_splines) {
    // This is the segment we want for this parameter variation
    // for this particular splineNum; 0 = MACCQE, 1 = pFC, 2 = EBC, etc

    //CW: Which Parameter we are accessing
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
    const float fB = gpu_coeff_many[CurrentKnotPos + 1];
    const float fC = gpu_coeff_many[CurrentKnotPos + 2];
    const float fD = gpu_coeff_many[CurrentKnotPos + 3];
    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float dx = val_gpu[Param] - tex1Dfetch<float>(text_coeff_x, segment_X);

    //CW: Wooow, let's use some fancy intrinsics and pull down the processing time by <1% from normal multiplication! HURRAY
    gpu_weights[splineNum] = fmaf(dx, fmaf(dx, fmaf(dx, fD, fC), fB), fY);
    // Or for the more "easy to read" version:
    //gpu_weights[splineNum] = (fY+dx*(fB+dx*(fC+dx*fD)));

    //#ifdef DEBUG
    //printf("splineNum = %i/%i, paramNo = %i, variation = %f, segment = %i, fX = %f, fX+1 = %f, dx = %f, d_n_splines = %i, d_spline_size = %i, weight = %f \n", splineNum, d_n_splines, gpu_paramNo_arr[splineNum], val_gpu[Param], segment, tex1Dfetch<float>(text_coeff_x, segment_X), tex1Dfetch<float>(text_coeff_x, segment_X+1), dx, d_n_splines, d_spline_size, gpu_weights[splineNum]);
    //#endif
  }
}

//*********************************************************
// Evaluate the TF1 on the GPU Using 5th order polynomial
__global__ void EvalOnGPU_TF1(
    const float* __restrict__ gpu_coeffs_tf1,
    const short int* __restrict__ gpu_paramNo_arr_tf1,
    float* __restrict__ gpu_weights_tf1) {
//*********************************************************
  // points per spline is the offset to skip in the index to move between splines
  const unsigned int tf1Num = (blockIdx.x * blockDim.x + threadIdx.x);

  if (tf1Num < d_n_TF1) {
    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float x = val_gpu[gpu_paramNo_arr_tf1[tf1Num]];

    // Read the coefficients
    const unsigned int TF1_Index = tf1Num * _nTF1Coeff_;
    const float a = gpu_coeffs_tf1[TF1_Index];
    const float b = gpu_coeffs_tf1[TF1_Index+1];

    gpu_weights_tf1[tf1Num] = fmaf(a, x, b);

    // gpu_weights_tf1[tf1Num] = a*x + b;
    //gpu_weights_tf1[tf1Num] = 1 + a*x + b*x*x + c*x*x*x + d*x*x*x*x + e*x*x*x*x*x;
  }
}

#ifndef Weight_On_SplineBySpline_Basis
//*********************************************************
// KS: Evaluate the total spline event weight on the GPU, as in most cases GPU is faster, even more this significant reduce memory transfer from GPU to CPU
__global__ void EvalOnGPU_TotWeight(
  const float* __restrict__ gpu_weights,
  const float* __restrict__ gpu_weights_tf1,

  float* __restrict__ gpu_total_weights,

  const cudaTextureObject_t __restrict__ text_nParamPerEvent,
  const cudaTextureObject_t __restrict__ text_nParamPerEvent_TF1) {
//*********************************************************
  const unsigned int EventNum = (blockIdx.x * blockDim.x + threadIdx.x);

  //KS: Accessing shared memory is much much faster than global memory hence we use shared memory for calculation and then write to global memory
  __shared__ float shared_total_weights[_BlockSize_];
  if(EventNum < d_n_events) //stopping condition
  {
    shared_total_weights[threadIdx.x] = 1.f;

    const unsigned int EventOffset = 2 * EventNum;

    for (unsigned int id = 0; id < tex1Dfetch<unsigned int>(text_nParamPerEvent, EventOffset); ++id) {
      shared_total_weights[threadIdx.x] *= gpu_weights[tex1Dfetch<unsigned int>(text_nParamPerEvent, EventOffset+1) + id];
    }

    for (unsigned int id = 0; id < tex1Dfetch<unsigned int>(text_nParamPerEvent_TF1, EventOffset); ++id) {
      shared_total_weights[threadIdx.x] *= gpu_weights_tf1[tex1Dfetch<unsigned int>(text_nParamPerEvent_TF1, EventOffset+1) + id];
    }
    gpu_total_weights[EventNum] = shared_total_weights[threadIdx.x];
  }
}
#endif

// *****************************************
// Run the GPU code for the separate many arrays. As in separate {x}, {y,b,c,d} arrays
// Pass the segment and the parameter values
// (binary search already performed in samplePDFND::FindSplineSegment()
__host__ void SMonolithGPU::RunGPU_SplineMonolith(
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
    float* cpu_weights_tf1,
#else
    float* cpu_total_weights,
#endif
    // Holds the changes in parameters
    float *vals,
    // Holds the segments for parameters
    short int *segment,
    const unsigned int h_n_splines,
    const unsigned int h_n_tf1) {
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

  // KS: Consider asynchronous kernel call, this might help EvalOnGPU_Splines and EvalOnGPU_TF1 are independent
  // Set the cache config to prefer L1 for the kernel
  //cudaFuncSetCacheConfig(EvalOnGPU_Splines, cudaFuncCachePreferL1);
  EvalOnGPU_Splines<<<grid_size, block_size>>>(
    gpu_paramNo_arr,
    gpu_nKnots_arr,

    gpu_coeff_many,

    gpu_weights,
    text_coeff_x
  );
  CudaCheckError();

  grid_size.x = (h_n_tf1 / block_size.x) + 1;
  EvalOnGPU_TF1<<<grid_size, block_size>>>(
    gpu_coeff_TF1_many,
    gpu_paramNo_TF1_arr,

    gpu_weights_tf1
  );
  CudaCheckError();

//KS: We can either copy gpu_weight and calculate total weight in reweighting loop, or not copy and calculate total weight stall at GPU, which means less memory transfer
#ifdef Weight_On_SplineBySpline_Basis
  // Here we have to make a somewhat large GPU->CPU transfer because it's all the splines' response
  cudaMemcpy(cpu_weights, gpu_weights, h_n_splines*sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();

  cudaMemcpy(cpu_weights_tf1, gpu_weights_tf1, h_n_tf1*sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();

//KS: Else calculate Total Weight
#else
  grid_size.x = (h_n_events / block_size.x) + 1;

  EvalOnGPU_TotWeight<<<grid_size, block_size>>>(
      gpu_weights,
      gpu_weights_tf1,

      gpu_total_weights,

      text_nParamPerEvent,
      text_nParamPerEvent_TF1
      );
  CudaCheckError();

  //KS: Here we have to make a somewhat large GPU->CPU transfer because it is proportional to number of events
  //KS: Normally code wait for memory transfer to finish before moving further cudaMemcpyAsync means we will continue to execute code and in a meantime keep copying stuff.
  cudaMemcpyAsync(cpu_total_weights, gpu_total_weights, h_n_events * sizeof(float), cudaMemcpyDeviceToHost, 0);
  CudaCheckError();
#endif

  #ifdef DEBUG
    printf("Copied GPU total weights to CPU with SUCCESS (drink more tea)\n");
    printf("Released calculated response from GPU with SUCCESS (drink most tea)\n");
  #endif
}

// *********************************
// CLEANING
// *********************************

// *********************************
// Clean up the {x},{ybcd} arrays
__host__ void SMonolithGPU::CleanupGPU_SplineMonolith(
    #ifndef Weight_On_SplineBySpline_Basis
    float *cpu_total_weights
    #endif
){
// *********************************
  cudaFree(gpu_paramNo_arr);
  cudaFree(gpu_nKnots_arr);

  // free the coefficient arrays
  cudaDestroyTextureObject(text_coeff_x);

  cudaFree(gpu_coeff_x);
  cudaFree(gpu_coeff_many);

  cudaFree(gpu_coeff_TF1_many);
  cudaFree(gpu_paramNo_TF1_arr);
  // free weights on the gpu
  cudaFree(gpu_weights);
  cudaFree(gpu_weights_tf1);
#ifndef Weight_On_SplineBySpline_Basis
  cudaFree(gpu_total_weights);
  //KS: Before removing variable let's destroy texture
  cudaDestroyTextureObject(text_nParamPerEvent);
  cudaDestroyTextureObject(text_nParamPerEvent_TF1);

  cudaFree(gpu_nParamPerEvent);
  cudaFree(gpu_nParamPerEvent_TF1);
  cudaFreeHost(cpu_total_weights);
  cpu_total_weights = nullptr;
#endif
  return;
}

// *******************************************
/// Clean up pinned variables at CPU
__host__ void SMonolithGPU::CleanupGPU_Segments(short int *segment, float *vals) {
// *******************************************
  cudaFreeHost(segment);
  cudaFreeHost(vals);

  segment = nullptr;
  vals = nullptr;

  return;
}
