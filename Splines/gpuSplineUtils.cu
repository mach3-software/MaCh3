//MaCh3 included
#include "Splines/gpuSplineUtils.cuh"

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


// *****************************************
// Make sure all Cuda threads finished execution
__host__ void SynchroniseSplines() {
  cudaDeviceSynchronize();
}

// *******************************************
//              INITIALISE GPU
// *******************************************

SplineMonolithGPU::SplineMonolithGPU() {
  /// Number of events living on CPU
  cpu_n_params = -1;
  cpu_n_events = -1;
  cpu_n_TF1 = 0;
  cpu_n_splines = 0;
  cpu_spline_size = 0;

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

// *******************************************
SplineMonolithGPU::~SplineMonolithGPU() {
// *******************************************
  cudaFree(gpu_paramNo_arr);
  cudaFree(gpu_nKnots_arr);

  // free the coefficient arrays
  cudaDestroyTextureObject(text_coeff_x);
  cudaFree(gpu_coeff_x);
  cudaFree(gpu_coeff_many);

  cudaFree(gpu_par_val);
  cudaFree(gpu_spline_segment);

  cudaFree(gpu_coeff_TF1_many);
  cudaFree(gpu_paramNo_TF1_arr);
  // free weights on the gpu
  cudaFree(gpu_weights);
  cudaFree(gpu_weights_tf1);
  cudaFree(gpu_total_weights);
  //KS: Before removing variable let's destroy texture
  cudaDestroyTextureObject(text_nParamPerEvent);
  cudaDestroyTextureObject(text_nParamPerEvent_TF1);

  cudaFree(gpu_nParamPerEvent);
  cudaFree(gpu_nParamPerEvent_TF1);
}

// *******************************************
// Initialiser when using the x array and combined y,b,c,d array
__host__ void SplineMonolithGPU::InitGPU_SplineMonolith(
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

  // Print allocation info to user
  printf("Allocated %i entries for paramNo and nKnots arrays, size = %f MB\n",
         n_splines, static_cast<double>(sizeof(short int) * n_splines + sizeof(unsigned int) * n_splines) / 1.0e6);
  printf("Allocated %i entries for x coeff arrays, size = %f MB\n",
         Eve_size, static_cast<double>(sizeof(float) * Eve_size) / 1.0e6);
  printf("Allocated %i entries for {ybcd} coeff arrays, size = %f MB\n",
         _nCoeff_ * total_nknots, static_cast<double>(sizeof(float) * _nCoeff_ * total_nknots) / 1.0e6);
  printf("Allocated %i entries for TF1 coefficient arrays, size = %f MB\n",
         _nTF1Coeff_ * n_tf1, static_cast<double>(sizeof(float) * _nTF1Coeff_ * n_tf1) / 1.0e6);

  //KS: Ask CUDA about memory usage
  checkGpuMem();
  PrintNdevices();
}

// *******************************************
__host__ void SplineMonolithGPU::InitGPU_Unbinned_SplineMonolith(
    M3::float_t **cpu_total_weights,
    int n_events) {
// *******************************************
  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) cpu_total_weights, n_events*sizeof(M3::float_t));
  CudaCheckError();

  //KS: Allocate memory for the array of total weights to be returned to CPU
  cudaMalloc((void **) &gpu_total_weights, n_events*sizeof(M3::float_t));
  CudaCheckError();

  //KS: Allocate memory for the map keeping track how many splines each parameter has
  cudaMalloc((void **) &gpu_nParamPerEvent, 2*n_events*sizeof(unsigned int));
  CudaCheckError();

  //KS: Allocate memory for the map keeping track how many TF1 each parameter has
  cudaMalloc((void **) &gpu_nParamPerEvent_TF1, 2*n_events*sizeof(unsigned int));
  CudaCheckError();
}

// *******************************************
// Allocate memory for spline segments
__host__ void SplineMonolithGPU::InitGPU_Segments(short int **segment) {
// *******************************************
  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) segment, cpu_n_params*sizeof(short int));
  CudaCheckError();
}

// *******************************************
// Allocate memory for spline segments
__host__ void SplineMonolithGPU::InitGPU_Vals(float **vals) {
// *******************************************
  //KS: Rather than allocate memory in standard way this fancy cuda tool allows to pin host memory which make memory transfer faster
  cudaMallocHost((void **) vals, cpu_n_params*sizeof(float));
  CudaCheckError();
}

// ******************************************************
//                START COPY TO GPU
// ******************************************************

// ******************************************************
// Copy to GPU for x array and separate ybcd array
__host__ void SplineMonolithGPU::CopyToGPU_SplineMonolith_Unbinned(
                            const SplineMonoStruct* cpu_spline_handler,

                            // TFI related now
                            const std::vector<float>& cpu_many_array_TF1,
                            const std::vector<short int>& cpu_paramNo_arr_TF1,
                            const int n_events,
                            const std::vector<unsigned int>& cpu_nParamPerEvent,
                            // TFI related now
                            const std::vector<unsigned int>& cpu_nParamPerEvent_TF1,

                            const int n_params,
                            const unsigned int n_splines,
                            const short int spline_size,
                            const unsigned int total_nknots,
                            const unsigned int n_tf1) {
// ******************************************************
  // Write to the global statics (h_* denotes host stored variable)
  cpu_n_params  = n_params;
  // Number of events
  cpu_n_events  = n_events;
  // Total number of valid TF1 for all loaded events
  cpu_n_TF1     = n_tf1;
  // Total number of valid splines for all loaded events
  cpu_n_splines = n_splines;
  /// Size of splines living
  cpu_spline_size = spline_size;

  //CW: Allocate memory for the frequently copied objects
  cudaMalloc(&gpu_par_val, n_params * sizeof(float));
  cudaMalloc(&gpu_spline_segment, n_params * sizeof(short int));

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
  cudaCreateTextureObject(&text_coeff_x, &resDesc_coeff_x, &texDesc_coeff_x, nullptr);
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
  cudaCreateTextureObject(&text_nParamPerEvent, &resDesc_nParamPerEvent, &texDesc_nParamPerEvent, nullptr);
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
  cudaCreateTextureObject(&text_nParamPerEvent_TF1, &resDesc_nParamPerEvent_tf1, &texDesc_nParamPerEvent_tf1, nullptr);
  CudaCheckError();
}

// ******************************************************
// Copy to GPU for x array and separate ybcd array
__host__ void SplineMonolithGPU::CopyToGPU_SplineMonolith_Binned(
  M3::float_t *manycoeff_arr,
  M3::float_t *xcoeff_arr,
  const std::vector<short int>& uniquesplinevec_Monolith,
  const std::vector<unsigned int>& coeffindexvec,

  const int n_params,
  const unsigned int n_splines,
  const short int spline_size,
  const unsigned int total_nknots) {
// ******************************************************
  // Write to the global statics (h_* denotes host stored variable)
  cpu_n_params  = n_params;
  // Total number of valid splines for all loaded events
  cpu_n_splines = n_splines;
  /// Size of splines living
  cpu_spline_size = spline_size;


  //CW: Allocate memory for the frequently copied objects
  cudaMalloc(&gpu_par_val, n_params * sizeof(float));
  cudaMalloc(&gpu_spline_segment, n_params * sizeof(short int));
  #ifndef _LOW_MEMORY_STRUCTS_
  cpu_tmp_weights.resize(cpu_n_splines);
  // Convert M3::float_t -> float for GPU
  std::vector<float> coeff_many_float(manycoeff_arr, manycoeff_arr + total_nknots * _nCoeff_);
  std::vector<float> coeff_x_float(xcoeff_arr, xcoeff_arr + spline_size * n_params);

  // Copy the coefficient arrays to the GPU; this only happens once per entire Markov Chain so is OK to do multiple extensive memory copies
  cudaMemcpy(gpu_coeff_many, coeff_many_float.data(), sizeof(float)*total_nknots*_nCoeff_, cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_coeff_x, coeff_x_float.data(), sizeof(float)*spline_size*n_params, cudaMemcpyHostToDevice);
  CudaCheckError();
  #else
  // Copy the coefficient arrays to the GPU; this only happens once per entire Markov Chain so is OK to do multiple extensive memory copies
  cudaMemcpy(gpu_coeff_many, manycoeff_arr, sizeof(float)*total_nknots*_nCoeff_, cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_coeff_x, xcoeff_arr, sizeof(float)*spline_size*n_params, cudaMemcpyHostToDevice);
  CudaCheckError();
  #endif
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
  cudaCreateTextureObject(&text_coeff_x, &resDesc_coeff_x, &texDesc_coeff_x, nullptr);
  CudaCheckError();

  // Also copy the parameter number for each spline onto the GPU; i.e. what spline parameter are we calculating right now
  cudaMemcpy(gpu_paramNo_arr, uniquesplinevec_Monolith.data(), n_splines*sizeof(short int), cudaMemcpyHostToDevice);
  CudaCheckError();

  // Also copy the knot map for each spline onto the GPU;
  cudaMemcpy(gpu_nKnots_arr, coeffindexvec.data(), n_splines*sizeof(unsigned int), cudaMemcpyHostToDevice);
  CudaCheckError();
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
  const unsigned int gpu_n_splines,
  const short int gpu_spline_size,
  const short int* __restrict__ gpu_paramNo_arr,
  const unsigned int* __restrict__ gpu_nKnots_arr,
  const float* __restrict__ gpu_coeff_many,
  const float* __restrict__ gpu_par_val,
  const short int* __restrict__ gpu_spline_segment,
  float* __restrict__ gpu_weights,
  const cudaTextureObject_t __restrict__ text_coeff_x) {
//*********************************************************
  // points per spline is the offset to skip in the index to move between splines
  const unsigned int splineNum = (blockIdx.x * blockDim.x + threadIdx.x);

  // this is the stopping condition!
  if (splineNum < gpu_n_splines) {
    // This is the segment we want for this parameter variation
    // for this particular splineNum; 0 = MACCQE, 1 = pFC, 2 = EBC, etc

    //CW: Which Parameter we are accessing
    const short int Param = gpu_paramNo_arr[splineNum];

    //CW: Avoids doing costly binary search on GPU
    const short int segment = gpu_spline_segment[Param];

    //KS: Segment for coeff_x is simply parameter*max knots + segment as each parmeters has the same spacing
    const short int segment_X = Param*gpu_spline_size+segment;

    //KS: Find knot position in out monolitical structure
    const unsigned int CurrentKnotPos = gpu_nKnots_arr[splineNum]*_nCoeff_+segment*_nCoeff_;

    // We've read the segment straight from CPU and is saved in gpu_spline_segment
    // polynomial parameters from the monolithic splineMonolith
    const float fY = gpu_coeff_many[CurrentKnotPos];
    const float fB = gpu_coeff_many[CurrentKnotPos + 1];
    const float fC = gpu_coeff_many[CurrentKnotPos + 2];
    const float fD = gpu_coeff_many[CurrentKnotPos + 3];
    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float dx = gpu_par_val[Param] - tex1Dfetch<float>(text_coeff_x, segment_X);

    //CW: Wooow, let's use some fancy intrinsics and pull down the processing time by <1% from normal multiplication! HURRAY
    gpu_weights[splineNum] = fmaf(dx, fmaf(dx, fmaf(dx, fD, fC), fB), fY);
    // Or for the more "easy to read" version:
    //gpu_weights[splineNum] = (fY+dx*(fB+dx*(fC+dx*fD)));

    //#ifdef MACH3_DEBUG
    //printf("splineNum = %i/%i, paramNo = %i, variation = %f, segment = %i, fX = %f, fX+1 = %f, dx = %f, gpu_n_splines = %i, gpu_spline_size = %i, weight = %f \n", splineNum, gpu_n_splines, gpu_paramNo_arr[splineNum], gpu_par_val[Param], segment, tex1Dfetch<float>(text_coeff_x, segment_X), tex1Dfetch<float>(text_coeff_x, segment_X+1), dx, gpu_n_splines, gpu_spline_size, gpu_weights[splineNum]);
    //#endif
  }
}

//*********************************************************
// Evaluate the TF1 on the GPU Using 5th order polynomial
__global__ void EvalOnGPU_TF1(
    const unsigned int gpu_n_TF1,
    const float* __restrict__ gpu_coeffs_tf1,
    const short int* __restrict__ gpu_paramNo_arr_tf1,
    const float* __restrict__ gpu_par_val,
    float* __restrict__ gpu_weights_tf1) {
//*********************************************************
  // points per spline is the offset to skip in the index to move between splines
  const unsigned int tf1Num = (blockIdx.x * blockDim.x + threadIdx.x);

  if (tf1Num < gpu_n_TF1) {
    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float x = gpu_par_val[gpu_paramNo_arr_tf1[tf1Num]];

    // Read the coefficients
    const unsigned int TF1_Index = tf1Num * _nTF1Coeff_;
    const float a = gpu_coeffs_tf1[TF1_Index];
    const float b = gpu_coeffs_tf1[TF1_Index+1];

    gpu_weights_tf1[tf1Num] = fmaf(a, x, b);
    // gpu_weights_tf1[tf1Num] = a*x + b;
    // gpu_weights_tf1[tf1Num] = 1 + a*x + b*x*x + c*x*x*x + d*x*x*x*x + e*x*x*x*x*x;
  }
}

//*********************************************************
// KS: Evaluate the total spline event weight on the GPU, as in most cases GPU is faster, even more this significant reduce memory transfer from GPU to CPU
__global__ void EvalOnGPU_TotWeight(
  const int gpu_n_events,
  const float* __restrict__ gpu_weights,
  const float* __restrict__ gpu_weights_tf1,

  M3::float_t* __restrict__ gpu_total_weights,

  const cudaTextureObject_t __restrict__ text_nParamPerEvent,
  const cudaTextureObject_t __restrict__ text_nParamPerEvent_TF1) {
//*********************************************************
  const unsigned int EventNum = (blockIdx.x * blockDim.x + threadIdx.x);

  if(EventNum < gpu_n_events) //stopping condition
  {
    float local_total_weight = 1.f;

    const unsigned int EventOffset = 2 * EventNum;

    for (unsigned int id = 0; id < tex1Dfetch<unsigned int>(text_nParamPerEvent, EventOffset); ++id) {
      local_total_weight *= gpu_weights[tex1Dfetch<unsigned int>(text_nParamPerEvent, EventOffset+1) + id];
    }

    for (unsigned int id = 0; id < tex1Dfetch<unsigned int>(text_nParamPerEvent_TF1, EventOffset); ++id) {
      local_total_weight *= gpu_weights_tf1[tex1Dfetch<unsigned int>(text_nParamPerEvent_TF1, EventOffset+1) + id];
    }
    gpu_total_weights[EventNum] = static_cast<M3::float_t>(local_total_weight);
  }
}

// *****************************************
// Run the GPU code for the separate many arrays. As in separate {x}, {y,b,c,d} arrays
// Pass the segment and the parameter values
// (binary search already performed in SplineBase::FindSplineSegment()
__host__ void SplineMonolithGPU::RunGPU_SplineMonolith_Unbinned(
    M3::float_t* cpu_total_weights,
    // Holds the changes in parameters
    float *vals,
    // Holds the segments for parameters
    short int *segment) {
// *****************************************
  dim3 block_size;
  dim3 grid_size;

  block_size.x = _BlockSize_;
  grid_size.x = (cpu_n_splines / block_size.x) + 1;

  // Copy the segment values to the GPU (segment_gpu), which is cpu_n_params long
  cudaMemcpy(gpu_spline_segment, segment, cpu_n_params * sizeof(short int), cudaMemcpyHostToDevice);
  CudaCheckError();

  // Copy the parameter values values to the GPU (vals_gpu), which is cpu_n_params long
  cudaMemcpy(gpu_par_val, vals, cpu_n_params * sizeof(float), cudaMemcpyHostToDevice);
  CudaCheckError();

  // KS: Consider asynchronous kernel call, this might help EvalOnGPU_Splines and EvalOnGPU_TF1 are independent
  // Set the cache config to prefer L1 for the kernel
  //cudaFuncSetCacheConfig(EvalOnGPU_Splines, cudaFuncCachePreferL1);
  EvalOnGPU_Splines<<<grid_size, block_size>>>(
    cpu_n_splines,
    cpu_spline_size,
    gpu_paramNo_arr,
    gpu_nKnots_arr,
    gpu_coeff_many,
    gpu_par_val,
    gpu_spline_segment,

    gpu_weights,
    text_coeff_x
  );
  CudaCheckError();

  grid_size.x = (cpu_n_TF1 / block_size.x) + 1;
  EvalOnGPU_TF1<<<grid_size, block_size>>>(
    cpu_n_TF1,
    gpu_coeff_TF1_many,
    gpu_paramNo_TF1_arr,
    gpu_par_val,
    gpu_weights_tf1
  );
  CudaCheckError();

  grid_size.x = (cpu_n_events / block_size.x) + 1;

  EvalOnGPU_TotWeight<<<grid_size, block_size>>>(
      cpu_n_events,
      gpu_weights,
      gpu_weights_tf1,

      gpu_total_weights,

      text_nParamPerEvent,
      text_nParamPerEvent_TF1
      );
  CudaCheckError();

  //KS: Here we have to make a somewhat large GPU->CPU transfer because it is proportional to number of events
  //KS: Normally code wait for memory transfer to finish before moving further cudaMemcpyAsync means we will continue to execute code and in a meantime keep copying stuff.
  cudaMemcpyAsync(cpu_total_weights, gpu_total_weights, cpu_n_events * sizeof(M3::float_t), cudaMemcpyDeviceToHost, 0);
  CudaCheckError();

  #ifdef MACH3_DEBUG
    printf("Copied GPU total weights to CPU with SUCCESS (drink more tea)\n");
    printf("Released calculated response from GPU with SUCCESS (drink most tea)\n");
  #endif
}

// *****************************************
// Run the GPU code for the separate many arrays. As in separate {x}, {y,b,c,d} arrays
// Pass the segment and the parameter values
// (binary search already performed in SplineBase::FindSplineSegment()
__host__ void SplineMonolithGPU::RunGPU_SplineMonolith_Binned(
  M3::float_t* cpu_spline_weights,
  // Holds the changes in parameters
  float *vals,
  // Holds the segments for parameters
  short int *segment) {
// *****************************************
  dim3 block_size;
  dim3 grid_size;

  block_size.x = _BlockSize_;
  grid_size.x = (cpu_n_splines / block_size.x) + 1;

  // Copy the segment values to the GPU (segment_gpu), which is cpu_n_params long
  cudaMemcpy(gpu_spline_segment, segment, cpu_n_params * sizeof(short int), cudaMemcpyHostToDevice);
  CudaCheckError();

  // Copy the parameter values values to the GPU (vals_gpu), which is cpu_n_params long
  cudaMemcpy(gpu_par_val, vals, cpu_n_params * sizeof(float), cudaMemcpyHostToDevice);
  CudaCheckError();

  // KS: Consider asynchronous kernel call, this might help EvalOnGPU_Splines and EvalOnGPU_TF1 are independent
  // Set the cache config to prefer L1 for the kernel
  //cudaFuncSetCacheConfig(EvalOnGPU_Splines, cudaFuncCachePreferL1);
  EvalOnGPU_Splines<<<grid_size, block_size>>>(
    cpu_n_splines,
    cpu_spline_size,
    gpu_paramNo_arr,
    gpu_nKnots_arr,
    gpu_coeff_many,
    gpu_par_val,
    gpu_spline_segment,

    gpu_weights,
    text_coeff_x
  );
  CudaCheckError();

  #ifndef _LOW_MEMORY_STRUCTS_
  // Here we have to make a somewhat large GPU->CPU transfer because it's all the splines' response
  cudaMemcpy(cpu_tmp_weights.data(), gpu_weights, cpu_n_splines * sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();
  // KS: we need to perform conversion from float (this is what GPU uses) to M3::float_t
  /// @todo add maybe some ifdef to avoid doing this when we actually have float
  for (unsigned int i = 0; i < cpu_n_splines; ++i) {
    cpu_spline_weights[i] = static_cast<M3::float_t>(cpu_tmp_weights[i]);
  }
  #else
  // Here we have to make a somewhat large GPU->CPU transfer because it's all the splines' response
  /// @todo KS: Make cpu_spline_weight pinned for faster load
  cudaMemcpy(cpu_spline_weights, gpu_weights, cpu_n_splines * sizeof(float), cudaMemcpyDeviceToHost);
  CudaCheckError();
  #endif

  #ifdef MACH3_DEBUG
  printf("Copied GPU total weights to CPU with SUCCESS (drink more tea)\n");
  printf("Released calculated response from GPU with SUCCESS (drink most tea)\n");
  #endif
}

// *******************************************
/// Clean up pinned variables at CPU
__host__ void SplineMonolithGPU::CleanupPinnedMemory(M3::float_t *cpu_total_weights,
                                                     short int *segment, float *vals) {
// *******************************************
  if(cpu_total_weights != nullptr) cudaFreeHost(cpu_total_weights);
  cudaFreeHost(segment);
  cudaFreeHost(vals);

  cpu_total_weights = nullptr;
  segment = nullptr;
  vals = nullptr;
}
