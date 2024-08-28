/// MaCh3 event-by-event cross-section spline code
/// Written by Richard Calland, Asher Kaboth, Clarence Wret, Kamil Skwarczynski
///
/// Contains code to run on CUDA GPUs. Essentially we load up stripped TSpline3 objects to the GPU and do the equivalent of TSpline3->Eval(double) for all events
/// Now also supports TF1 evals
/// Called from samplePDF/samplePDFND.cpp -> splines/SplineMonolith.cpp -> splines/gpuSplineUtils.cu

//MaCh3 included
#include "manager/gpuUtils.cuh"
#include "splines/SplineCommon.h"

/// @brief Allocate memory on gpu for spline monolith
__host__ void InitGPU_SplineMonolith(
                          float **gpu_x_array,
                          float **gpu_many_array,
                          float **gpu_weights,

                          short int** gpu_paramNo_arr,
                          unsigned int** gpu_nKnots_arr,

                          float **gpu_many_TF1_array,
                          float **gpu_weights_tf1,
                          short int** gpu_paramNo_TF1_arr,
                 #ifndef Weight_On_SplineBySpline_Basis
                          float **cpu_total_weights,
                          float **gpu_total_weights,
                          int n_events,
                          unsigned int** gpu_nParamPerEvent,
                          unsigned int** gpu_nParamPerEvent_tf1,
                  #endif
                          unsigned int total_nknots,
                          unsigned int n_splines,
                          unsigned int n_tf1,
                          int Eve_size);


/// @brief Allocate memory for spline segments
/// @param segment Found spline segment for each parameter
__host__ void InitGPU_Segments(short int **segment);

/// @brief Allocate memory for spline segments
/// @param vals Value to which we want reweight for each parameter
__host__ void InitGPU_Vals(float **vals);


/// @brief Copy to GPU for x array and separate ybcd array
__host__ void CopyToGPU_SplineMonolith(
                            short int *gpu_paramNo_arr,
                            unsigned int *gpu_nKnots_arr,
                            float *gpu_x_array,
                            float *gpu_many_array,

                            std::vector<short int> cpu_paramNo_arr,
                            std::vector<unsigned int> cpu_nKnots_arr,
                            std::vector<float> cpu_x_array,
                            std::vector<float> cpu_many_array,
                            // TFI related now
                            float *gpu_many_TF1_array,
                            short int* gpu_paramNo_arr_TF1,

                            std::vector<float> cpu_many_array_TF1,
                            std::vector<short int> cpu_paramNo_arr_TF1,
                            #ifndef Weight_On_SplineBySpline_Basis
                            int n_events,
                            std::vector<unsigned int> cpu_nParamPerEvent,
                            unsigned int *gpu_nParamPerEvent,
                            // TFI related now
                            std::vector<unsigned int> cpu_nParamPerEvent_TF1,
                            unsigned int *gpu_nParamPerEvent_TF1,
                            #endif
                            int n_params,
                            unsigned int n_splines,
                            short int spline_size,
                            unsigned int total_nknots,
                            unsigned int n_tf1);

/// @brief Evaluate the spline on the GPU Using one {y,b,c,d} array and one {x} array
/// Should be most efficient at cache hitting and memory coalescence
/// But using spline segments rather than the parameter value: avoids doing binary search on GPU
/// @param gpu_paramNo_arr has length = spln_counter (keeps track of which parameter we're using on this thread)
/// @param gpu_nKnots_arr has length = spln_counter (keeps track where current spline starts)
/// @param gpu_coeff_many has length = nKnots * 4
/// @param gpu_weights has length = spln_counter * spline_size
/// @param text_coeff_x array storing info about X coeff, uses texture memory. Has length = n_params * spline_size,
__global__ void EvalOnGPU_Splines(
    const short int* __restrict__ gpu_paramNo_arr,
    const unsigned int* __restrict__ gpu_nKnots_arr,
    const float* __restrict__ gpu_coeff_many,
    float *gpu_weights,
    const cudaTextureObject_t __restrict__ text_coeff_x);

/// @brief Evaluate the TF1 on the GPU Using 5th order polynomial
/// @param gpu_coeffs_tf1 coefficents of TF1, has length = tf1 coeef counter
/// @param gpu_paramNo_arr_tf1 has length = spln_counter (keeps track of which parameter we're using on this thread)
/// @param gpu_weights_tf1 has length = spln_counter * spline_size
__global__ void EvalOnGPU_TF1( 
    const float* __restrict__ gpu_coeffs_tf1,
    const short int* __restrict__ gpu_paramNo_arr_tf1,
    float *gpu_weights_tf1);

#ifndef Weight_On_SplineBySpline_Basis
/// @brief KS: Evaluate the total spline event weight on the GPU, as in most cases GPU is faster, even more this significant reduce memory transfer from GPU to CPU
/// @param gpu_weights Weight for each spline object
/// @param gpu_weights_tf1 Weight for each TF1 object
/// @param gpu_total_weights Total weight for each event
/// @param text_nParamPerEvent map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
/// @param text_nParamPerEvent_TF1 map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
__global__ void EvalOnGPU_TotWeight(
  const float* __restrict__ gpu_weights,
  const float* __restrict__ gpu_weights_tf1,

  float *gpu_total_weights,

  const cudaTextureObject_t __restrict__ text_nParamPerEvent,
  const cudaTextureObject_t __restrict__ text_nParamPerEvent_TF1);
#endif

/// @brief Run the GPU code for the separate many arrays. As in separate {x}, {y,b,c,d} arrays
/// Pass the segment and the parameter values
/// (binary search already performed in samplePDFND::FindSplineSegment()
__host__ void RunGPU_SplineMonolith(
  const short int* gpu_paramNo_arr,
  const unsigned int* gpu_nKnots_arr,

  const float *gpu_coeff_many,

  const short int* gpu_paramNo_tf1_arr,
  const float *gpu_coeff_many_tf1,

  float* gpu_weights,
  float* gpu_weights_tf1,
  #ifdef Weight_On_SplineBySpline_Basis
  float* cpu_weights,
  float* cpu_weights_tf1,
  #else
  float* gpu_total_weights,
  float* cpu_total_weights,
  #endif
  // Holds the changes in parameters
  float *vals,
  // Holds the segments for parameters
  short int *segment,
  const unsigned int h_n_splines,
  const unsigned int h_n_tf1);

/// @brief Make sure all Cuda threads finished execution
__host__ void SynchroniseSplines();

/// @brief Clean up the {x},{ybcd} arrays
__host__ void CleanupGPU_SplineMonolith(
  short int *gpu_paramNo_arr,
  unsigned int *gpu_nKnots_arr,

  float *gpu_x_array,
  float *gpu_many_array,

  float *gpu_many_TF1_array,
  short int* gpu_paramNo_arr_TF1,
  #ifndef Weight_On_SplineBySpline_Basis
  float *gpu_total_weights,
  unsigned int *gpu_nParamPerEvent,
  unsigned int *gpu_nParamPerEvent_TF1,
  float *cpu_total_weights,
  #endif
  float *gpu_weights,
  float *gpu_weights_tf1);

/// @brief Clean up pinned variables at CPU
/// @param segment Found spline segment for each parameter
/// @param vals Value to which we want reweight for each parameter
__host__ void CleanupGPU_Segments(short int *segment, float *vals);
