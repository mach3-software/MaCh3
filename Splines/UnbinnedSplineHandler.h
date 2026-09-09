#pragma once

#include "Splines/SplineBase.h"

/// @brief Even-by-event class calculating response for spline parameters. It is possible to use GPU acceleration
/// @author Clarence Wret
/// @author Kamil Skwarczynski
class UnbinnedSplineHandler : public SplineBase {
  public:
    /// @brief Constructor
    /// @param MasterSpline Vector of TSpline3 pointers which we strip back
    /// @param SplineType Whether object is TSpline3 or TF1
    /// @param SaveFlatTree Whether we want to save monolith into speedy flat tree
    /// @param _FastSplineName Name to which spline file will be saved
    UnbinnedSplineHandler(std::vector<std::vector<TResponseFunction_red*> > &MasterSpline,
              const std::vector<RespFuncType> &SplineType,
              const bool SaveFlatTree = false,
              const std::string& _FastSplineName = "SplineFile.root",
              const bool Use_GPU = true);
    /// @brief Constructor where you pass path to preprocessed root FileName
    /// @param FileName path to pre-processed root file containing stripped monolith info
    UnbinnedSplineHandler(const std::string& FileName);
    /// @brief Destructor for UnbinnedSplineHandler class.
    virtual ~UnbinnedSplineHandler();

    /// @copydoc SplineBase::Evaluate
    void Evaluate() final;

    /// @brief Get class name
    std::string GetName() const override {return "SplineMonolith";};

    /// @brief KS: Get pointer to total weight to make fit faster wrooom!
    /// @param event Name event number in used MC
    /// @return Pointer to the total weight
    const M3::float_t* RetPointer(const int event) const {return &cpu_total_weights[event];}
    
    /// @brief KS: Set pointers to spline params
    /// @param spline_ParsPointers Vector of pointers to spline params
    void SetSplinePointers(std::vector< const M3::float_t* > spline_ParsPointers) {
      for (M3::int_t i = 0; i < nParams; ++i) SplineInfoArray[i].splineParsPointer = spline_ParsPointers[i];
    };
    
    /// @copydoc SplineBase::PrepareSplineFile
    void PrepareSplineFile(std::string FileName) final;
    /// @copydoc SplineBase::LoadSplineFile
    void LoadSplineFile(std::string FileName) final;
  private:
    /// @brief KS: Set everything to null etc.
    void Initialise();
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
    void ScanMasterSpline(std::vector<std::vector<TResponseFunction_red*> > & MasterSpline,
                                 unsigned int &nEvents,
                                 short int &MaxPoints,
                                 short int &numParams,
                                 int &nSplines,
                                 unsigned int &NSplinesValid,
                                 unsigned int &numKnots,
                                 unsigned int &nTF1Valid,
                                 unsigned int &nTF1_coeff,
                                 const std::vector<RespFuncType> &SplineType);
    /// @brief CW: Prepare the TSpline3_red objects for the GPU
    /// @param MasterSpline Vector of TResponseFunction_red pointers which we strip back
    void PrepareForGPU(std::vector<std::vector<TResponseFunction_red*> > &MasterSpline, const std::vector<RespFuncType> &SplineType);
    /// @brief CW: The shared initialiser from constructors of TResponseFunction_red
    void MoveToGPU();

    /// @brief KS: Print info about how much knots etc has been initialised
    void PrintInitialsiation() const;

    /// @brief CW: This loads up coefficients into two arrays: one x array and one yabcd array
    /// @brief CW: This should maximize our cache hits!
    /// @param spl pointer to TSpline3_red
    /// @param nPoints number of knots
    /// @param xArray array X value for each knot
    /// @param manyArray Array holding coefficients for each knot
    void GetSplineCoeff_SepMany(TSpline3_red* &spl, int &nPoints, float *&xArray, float *&manyArray) const;

    /// @copydoc SplineBase::CalcSplineWeights
    void CalcSplineWeights() final;
    /// @brief Calc total event weight
    void CalcTotalEventWeight();

    /// Number of events
    unsigned int NEvents;

    /// Number of valid TF1
    unsigned int NTF1_valid;

    /// Sum of all knots over all splines
    unsigned int nKnots;
    /// Sum of all coefficients over all TF1
    unsigned int nTF1coeff;

    /// CPU arrays to hold weight for each spline
    float *cpu_weights_spline_var;
    /// CPU arrays to hold weight for each TF1
    float *cpu_weights_tf1_var;
    /// KS: This holds the total CPU weights that gets read in SampleHandler
    M3::float_t* cpu_total_weights;

    /// KS: CPU map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
    std::vector<unsigned int> cpu_nParamPerEvent;

    /// KS: CPU map keeping track how many parameters applies to each event, we keep two numbers here {number of TF1 per event, index where TF1 start for a given event}
    std::vector<unsigned int> cpu_nParamPerEvent_tf1;

    /// KS: Store info about Spline monolith, this allow to obtain better step time. As all necessary information for spline weight calculation are here meaning better cache hits.
    SplineMonoStruct* cpu_monolith;

    /// CPU arrays to hold TF1 coefficients
    std::vector<float> cpu_coeff_TF1_many;

    /// CPU arrays to hold number of points
    std::vector<short int> cpu_nPoints_arr;

    /// CW: CPU array with the number of points per spline (not per spline point!)
    std::vector<short int> cpu_paramNo_TF1_arr;

    /// Flag telling whether we are saving spline monolith into handy root file
    bool SaveSplineFile;

    /// Name of Fast Spline to which will be saved
    std::string FastSplineName;
};
