/// @file gpuSplineUtils.cuh
/// @brief MaCh3 event-by-event cross-section spline code
/// @author Richard Calland
/// @author Asher Kaboth
/// @author Clarence Wret
/// @author Kamil Skwarczynski
///
/// Contains code to run on CUDA GPUs. Essentially we load up stripped TSpline3 objects to the GPU and do the equivalent of TSpline3->Eval(double) for all events
/// Now also supports TF1 evals
/// Called from samplePDF/samplePDFND.cpp -> splines/SplineMonolith.cpp -> splines/gpuSplineUtils.cu

//MaCh3 included
#include "manager/gpuUtils.cuh"
#include "splines/SplineCommon.h"


/// @brief Make sure all Cuda threads finished execution
__host__ void SynchroniseSplines();

/// @brief Evaluate the spline on the GPU Using one {y,b,c,d} array and one {x} array
/// Should be most efficient at cache hitting and memory coalescence
/// But using spline segments rather than the parameter value: avoids doing binary search on GPU
/// @param gpu_paramNo_arr has length = spln_counter (keeps track of which parameter we're using on this thread)
/// @param gpu_nKnots_arr has length = spln_counter (keeps track where current spline starts)
/// @param gpu_coeff_many has length = nKnots * 4, stores all coefficients for all splines and knots
/// @param gpu_weights has length = spln_counter * spline_size
/// @param text_coeff_x array storing info about X coeff, uses texture memory. Has length = n_params * spline_size,
__global__ void EvalOnGPU_Splines(
  const short int* __restrict__ gpu_paramNo_arr,
  const unsigned int* __restrict__ gpu_nKnots_arr,
  const float* __restrict__ gpu_coeff_many,
  float *gpu_weights,
  const cudaTextureObject_t __restrict__ text_coeff_x);

/// @brief Evaluate the TF1 on the GPU Using 5th order polynomial
/// @param gpu_coeffs_tf1 coefficients of TF1, has length = tf1 coeef counter
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

/// @brief Class responsible for calculating spline weight on GPU
class SMonolithGPU
{
  public:
    /// @brief constructor
    SMonolithGPU();
    /// @brief destructor
    virtual ~SMonolithGPU();

    /// @brief Allocate memory on gpu for spline monolith
    /// @param gpu_x_array Small array with X coefficients at GPU
    __host__ void InitGPU_SplineMonolith(
      #ifndef Weight_On_SplineBySpline_Basis
      float **cpu_total_weights,
      int n_events,
      #endif
      unsigned int total_nknots,
      unsigned int n_splines,
      unsigned int n_tf1,
      int Eve_size);

    /// @brief Copy to GPU for x array and separate ybcd array
    __host__ void CopyToGPU_SplineMonolith(
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
      unsigned int n_tf1);

    /// @brief Allocate memory for spline segments
    /// @param segment Found spline segment for each parameter
    __host__ void InitGPU_Segments(short int **segment);

    /// @brief Allocate memory for spline segments
    /// @param vals Value to which we want reweight for each parameter
    __host__ void InitGPU_Vals(float **vals);

    /// @brief Run the GPU code for the separate many arrays. As in separate {x}, {y,b,c,d} arrays
    /// Pass the segment and the parameter values
    /// (binary search already performed in SplineMonolith::FindSplineSegment()
    __host__ void RunGPU_SplineMonolith(
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
      const unsigned int h_n_tf1);

    /// @brief Clean up the {x},{ybcd} arrays
    __host__ void CleanupGPU_SplineMonolith(
      #ifndef Weight_On_SplineBySpline_Basis
      float *cpu_total_weights
      #endif
      );

    /// @brief Clean up pinned variables at CPU
    /// @param segment Found spline segment for each parameter
    /// @param vals Value to which we want reweight for each parameter
    __host__ void CleanupGPU_Segments(short int *segment, float *vals);

  private:
    /// KS: GPU map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
    unsigned int *gpu_nParamPerEvent;
    /// KS: GPU map keeping track how many parameters applies to each event, we keep two numbers here {number of TF1 per event, index where TF1 start for a given event}
    unsigned int *gpu_nParamPerEvent_TF1;

    /// KS: GPU arrays to hold X coefficient
    float *gpu_coeff_x;

    /// GPU arrays to hold other coefficients
    float *gpu_coeff_many;

    /// KS: GPU Number of knots per spline
    unsigned int *gpu_nKnots_arr;

    /// CW: GPU array with the number of points per spline (not per spline point!)
    short int *gpu_paramNo_arr;

    /// GPU arrays to hold TF1 coefficients
    float *gpu_coeff_TF1_many;
    /// GPU arrays to hold number of points
    short int *gpu_nPoints_arr;
    /// CW: GPU array with the number of points per TF1 object
    short int *gpu_paramNo_TF1_arr;

    /// GPU arrays to hold weight for event
    float *gpu_total_weights;
    /// GPU arrays to hold weight for each spline
    float *gpu_weights;
    /// GPU arrays to hold weight for each TF1
    float *gpu_weights_tf1;

    // h_NAME declares HOST constants (live on CPU)
    /// Number of params living on CPU
    int h_n_params;
    /// Number of events living on CPU
    int h_n_events;

    // ******************************************
    // TEXTURES
    // ******************************************
    /// KS: Textures are L1 cache variables which are well optimised for fetching. Make texture only for variables you often access but rarely overwrite. There are limits on texture memory so don't use huge arrays
    cudaTextureObject_t text_coeff_x = 0;
    #ifndef Weight_On_SplineBySpline_Basis
    /// KS: Map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
    cudaTextureObject_t text_nParamPerEvent = 0;
    /// KS: Map keeping track how many parameters applies to each event, we keep two numbers here {number of TF1 per event, index where TF1 start for a given event}
    cudaTextureObject_t text_nParamPerEvent_TF1 = 0;
    #endif
};

