#pragma once
#ifdef USE_FPGA
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#endif

#include "splines/SplineBase.h"

//KS: Joy of forward declaration https://gieseanw.wordpress.com/2018/02/25/the-joys-of-forward-declarations-results-from-the-real-world/
class SMonolithGPU;

/// @brief Even-by-event class calculating response for spline parameters. It is possible to use GPU acceleration
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/05.-Splines).
class SMonolith : public SplineBase {
  public:
    /// @brief Constructor
    /// @param MasterSpline Vector of TSpline3 pointers which we strip back
    /// @param SplineType Whether object is TSpline3 or TF1
    /// @param SaveFlatTree Whether we want to save monolith into speedy flat tree
    SMonolith(std::vector<std::vector<TResponseFunction_red*> > &MasterSpline,
              const std::vector<RespFuncType> &SplineType,
              const bool SaveFlatTree = false);
    /// @brief Constructor where you pass path to preprocessed root FileName
    /// @param FileName path to pre-processed root file containing stripped monolith info
    SMonolith(std::string FileName);
    /// @brief Destructor for SMonolith class.
    virtual ~SMonolith();

    /// @brief  CW: This Eval should be used when using two separate x,{y,a,b,c,d} arrays to store the weights; probably the best one here! Same thing but pass parameter spline segments instead of variations
    void Evaluate() override;

    /// @brief Get class name
    inline std::string GetName()const {return "SplineMonolith";};

    /// @brief KS: After calculations are done on GPU we copy memory to CPU. This operation is asynchronous meaning while memory is being copied some operations are being carried. Memory must be copied before actual reweight. This function make sure all has been copied.
    void SynchroniseMemTransfer();

    /// @brief KS: Get pointer to total weight to make fit faster wrooom!
    /// @param event Name event number in used MC
    /// @return Pointer to the total weight
    inline const float* retPointer(const int event) {return &cpu_total_weights[event];}
    
    /// @brief KS: Set pointers to spline params
    /// @param spline_ParsPointers Vector of pointers to spline params
    inline void setSplinePointers(std::vector< const double* > spline_ParsPointers) {
      splineParsPointer = spline_ParsPointers;
      for (_int_ i = 0; i < nParams; ++i) SplineInfoArray[i].splineParsPointer = spline_ParsPointers[i];
    };
    
    /// The returned gpu weights, read by the GPU
    float* cpu_weights;
    /// KS: This holds the total CPU weights that gets read in samplePDFND
    float *cpu_total_weights;

  private:
    /// @brief KS: Set everything to null etc.
    inline void Initialise();
    /// @brief CW: Function to scan through the MasterSpline of TSpline3
    /// @param MasterSpline Vector of TSpline3_red pointers which we strip back
    /// @param NEvents Number of MC events
    /// @param MaxPoints Maximal number of knots per splines
    /// @param numParams Total number of parameters
    /// @param numKnots Total number of knots, which is sum of individual knots per each spline
    /// @param nTF1_coeff Number of TF1 coefficients in all TF1 objects
    /// @param SplineType Whether object is TSpline3 or TF1
    /// @param NSplinesValid Total number of valid (not null) TSpline3
    /// @param nTF1Valid Total number of valid (not null) TF1
    inline void ScanMasterSpline(std::vector<std::vector<TResponseFunction_red*> > & MasterSpline,
                                 unsigned int &nEvents,
                                 int &MaxPoints,
                                 short int &numParams,
                                 int &nSplines,
                                 unsigned int &NSplinesValid,
                                 unsigned int &numKnots,
                                 unsigned int &nTF1Valid,
                                 unsigned int &nTF1_coeff,
                                 const std::vector<RespFuncType> &SplineType);
    /// @brief CW: Prepare the TSpline3_red objects for the GPU
    /// @param MasterSpline Vector of TResponseFunction_red pointers which we strip back
    inline void PrepareForGPU(std::vector<std::vector<TResponseFunction_red*> > &MasterSpline, const std::vector<RespFuncType> &SplineType);
    /// @brief CW: The shared initialiser from constructors of TResponseFunction_red
    inline void MoveToGPU();
        
    /// @brief KS: Print info about how much knots etc has been initialised
    inline void PrintInitialsiation();

    /// @brief CW: This loads up coefficients into two arrays: one x array and one yabcd array
    /// @brief CW: This should maximize our cache hits!
    /// @param spl pointer to TSpline3_red
    /// @param nPoints number of knots
    /// @param xArray array X value for each knot
    /// @param manyArray Array holding coefficients for each knot
    inline void getSplineCoeff_SepMany(TSpline3_red* &spl, int &nPoints, float *&xArray, float *&manyArray);

    /// @brief CW:Code used in step by step reweighting, Find Spline Segment for each param
    inline void FindSplineSegment() override;
    /// @brief CPU based code which eval weight for each spline
    inline void CalcSplineWeights() override;
    /// @brief Calc total event weight
    inline void ModifyWeights() override;
    /// @brief Conversion from valid splines to all
    inline void ModifyWeights_GPU();
    
    /// @brief KS: Prepare spline file that can be used for fast loading
    inline void PrepareSplineFile();
    /// @brief KS: Load preprocessed spline file
    /// @param FileName Path to ROOT file with predefined reduced Spline Monolith
    inline void LoadSplineFile(std::string FileName);

    /// Array of FastSplineInfo structs: keeps information on each xsec spline for fast evaluation
    /// Method identical to TSpline3::Eval(double) but faster because less operations
    FastSplineInfo *SplineInfoArray;
    /// Store currently found segment they are not in FastSplineInfo as in case of GPU we need to copy paste it to GPU
    short int *segments;
    /// Store parameter values they are not in FastSplineInfo as in case of GPU we need to copy paste it to GPU
    float *vals;
    /// This holds pointer to parameter position which we later copy paste it to GPU
    std::vector< const double* > splineParsPointer;

    /// Number of events
    unsigned int NEvents;
    /// Number of NIWG parameters that have splines
    short int nParams;
    /// Max knots for production
    int _max_knots;
    /// holds the index for good splines; don't do unsigned since starts with negative value!
    int *index_cpu;
    /// holds the index for good TF1; don't do unsigned since starts with negative value!
    int *index_TF1_cpu;

    /// Number of valid splines
    unsigned int NSplines_valid;
    /// Number of valid TF1
    unsigned int NTF1_valid;

    /// Number of total splines if each event had every parameter's spline
    unsigned int NSplines_total_large;

    /// Sum of all knots over all splines
    unsigned int nKnots;
    /// Sum of all coefficients over all TF1
    unsigned int nTF1coeff;

    /// CPU arrays to hold weight for each spline
    float *cpu_weights_var;

    /// CPU arrays to hold weight for each TF1
    float *cpu_weights_tf1_var;

    /// KS: CPU map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
    std::vector<unsigned int> cpu_nParamPerEvent;

    /// KS: CPU map keeping track how many parameters applies to each event, we keep two numbers here {number of TF1 per event, index where TF1 start for a given event}
    std::vector<unsigned int> cpu_nParamPerEvent_tf1;

    #ifdef USE_FPGA
    sycl::queue queue;
    SplineMonoUSM* cpu_spline_handler;
    float* cpu_coeff_TF1_many;
    
    short int* cpu_paramNo_TF1_arr;
    #else
    /// KS: Store info about Spline monolith, this allow to obtain better step time. As all necessary information for spline weight calculation are here meaning better cache hits.
    SplineMonoStruct* cpu_spline_handler;
    /// CPU arrays to hold TF1 coefficients
    std::vector<float> cpu_coeff_TF1_many;
    /// CW: CPU array with the number of points per spline (not per spline point!)
    std::vector<short int> cpu_paramNo_TF1_arr;
    #endif

    /// KS: Store info about Spline monolith, this allow to obtain better step time. As all necessary information for spline weight calculation are here meaning better cache hits.
    SMonolithGPU* gpu_spline_handler;



    /// CPU arrays to hold number of points
    std::vector<short int> cpu_nPoints_arr;



    /// Flag telling whether we are saving spline monolith into handy root file
    bool SaveSplineFile;
};
