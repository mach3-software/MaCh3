#pragma once

//MaCh3 includes
#include "samplePDF/Structs.h"
#include "splines/SplineBase.h"
#include "manager/MaCh3Modes.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TH3F.h"
_MaCh3_Safe_Include_End_ //}

/// @brief Bin-by-bin class calculating response for spline parameters.
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/05.-Splines).
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Henry Wallace
class splineFDBase : public SplineBase {
  /// @todo ETA - do all of these functions and members actually need to be public?
  public:
    /// @brief Constructor
    splineFDBase(covarianceXsec *xsec_, MaCh3Modes *Modes_);
    /// @brief Destructor
    virtual ~splineFDBase();

    /// @brief CW: This Eval should be used when using two separate x,{y,a,b,c,d} arrays
    /// to store the weights; probably the best one here! Same thing but pass parameter
    /// spline segments instead of variations
    void Evaluate();

    /// @brief add oscillation channel to spline monolith
    void AddSample(const std::string& SampleName,
                   const std::string& DetID,
                   const std::vector<std::string>& OscChanFileNames,
                   const std::vector<std::string>& SplineVarNames);
    /// @brief flatten multidimensional spline array into proper monolith
    void TransferToMonolith();
    /// @brief Remove setup variables not needed for spline evaluations
    void cleanUpMemory();

    /// @brief Loads and processes splines from ROOT files for a given sample.
    /// @note DB Add virtual so it can be overridden in experiment specific (if needed)
    virtual void FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames);
    /// @brief Check if there are any repeated modes. This is used to reduce the number
    /// of modes in case many interaction modes get averaged into one spline
    std::vector< std::vector<int> > StripDuplicatedModes(std::vector< std::vector<int> > InputVector);
    /// @brief Return the splines which affect a given event
    std::vector< std::vector<int> > GetEventSplines(std::string SampleName, int iOscChan, int EventMode, double Var1Val, double Var2Val, double Var3Val);

    /// @brief Grab histograms with spline binning
    std::vector<TAxis*> FindSplineBinning(std::string FileName, std::string SampleName);

    int CountNumberOfLoadedSplines(bool NonFlat=false, int Verbosity=0);
    std::string getDimLabel(const int BinningOpt, const unsigned int Axis) const;
    /// @brief Get index of sample based on name
    int getSampleIndex(const std::string& SampleName) const;
    /// @brief Ensure we have spline for a given bin
    bool isValidSplineIndex(const std::string& SampleName, int iSyst, int iOscChan, int iMode, int iVar1, int iVar2, int iVar3);

    void BuildSampleIndexingArray(const std::string& SampleName);
    void PrepForReweight();
    void getSplineCoeff_SepMany(int splineindex, M3::float_t *& xArray, M3::float_t *&manyArray);
    void PrintBinning(TAxis* Axis) const;
    /// @brief Print info like DetID number of spline params etc.
    void PrintSampleDetails(const std::string& SampleName) const;
    void PrintArrayDetails(const std::string& SampleName) const;

    /// @brief get pointer to spline weight based on bin variables
    const M3::float_t* retPointer(const int sample, const int oscchan, const int syst, const int mode,
                                  const int var1bin, const int var2bin, const int var3bin) const{
      int index = indexvec[sample][oscchan][syst][mode][var1bin][var2bin][var3bin];
      return &weightvec_Monolith[index];
    }

  protected:
    /// @brief CPU based code which eval weight for each spline
    void CalcSplineWeights() override;
    /// @brief Calc total event weight, not used by Bin-by-bin splines
    void ModifyWeights() override {return;};
    /// Pointer to covariance from which we get information about spline params
    covarianceXsec* xsec;

    //And now the actual member variables
    std::vector<std::string> SampleNames;
    std::vector<int> Dimensions;
    std::vector<std::vector<std::string>> DimensionLabels;
    std::vector<std::string> DetIDs;
    std::vector<int> nSplineParams;
    std::vector<int> nOscChans;

    std::vector< std::vector< std::vector<TAxis*> > > SplineBinning;
    std::vector< std::vector<std::string> > SplineFileParPrefixNames;
    /// A vector of vectors of the spline modes that a systematic applies to
    /// This gets compared against the event mode to figure out if a syst should
    /// apply to an event or not
    std::vector< std::vector< std::vector<int> > > SplineModeVecs;
    /// @brief This holds the global spline index and is used to grab the current parameter value
    /// to evaluate splines at. Each internal vector will be of size of the number of spline
    /// systematics which affect that sample.
    std::vector< std::vector<int> > GlobalSystIndex;
    /// @brief spline interpolation types for each sample. These vectors are from
    /// a call to GetSplineInterpolationFromDetID()
    std::vector< std::vector<SplineInterpolation> > SplineInterpolationTypes;

    /// name of each spline parameter
    std::vector<std::string> UniqueSystNames;
    /// Global index of each spline param, it allows us to match spline ordering with global
    std::vector<int> UniqueSystIndices;

    /// @brief Variables related to determined which modes have splines and which piggy-back of other modes
    std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< int > > > > > > > indexvec;
    std::vector<int > coeffindexvec;
    /// Unique coefficient indices
    std::vector<int> uniquecoeffindices;

    /// holds each spline object before stripping into coefficient monolith
    std::vector< TSpline3_red* > splinevec_Monolith;

    int MonolithSize;
    int MonolithIndex;
    int CoeffIndex;

    /// Need to keep track of which splines are flat and which aren't
    bool *isflatarray;
    /// x coefficients for each spline
    M3::float_t *xcoeff_arr;
    /// ybcd coefficients for each spline
    M3::float_t *manycoeff_arr;

    /// Stores weight from spline evaluation for each single spline
    std::vector<M3::float_t> weightvec_Monolith;
    /// Maps single spline object with single parameter
    std::vector<int> uniquesplinevec_Monolith;

    /// pointer to MaCh3 Mode from which we get spline suffix
    MaCh3Modes* Modes;
    enum TokenOrdering{kSystToken,kModeToken,kVar1BinToken,kVar2BinToken,kVar3BinToken,kNTokens};
    virtual std::vector<std::string> GetTokensFromSplineName(std::string FullSplineName) = 0;
};
