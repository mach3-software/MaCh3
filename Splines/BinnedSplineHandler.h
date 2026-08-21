#pragma once

//MaCh3 includes
#include "Splines/SplineBase.h"
#include "Manager/MaCh3Modes.h"

/// @brief Bin-by-bin class calculating response for spline parameters.
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Henry Wallace
class BinnedSplineHandler : public SplineBase {
  public:
    /// @brief Constructor
    BinnedSplineHandler(ParameterHandlerGeneric *ParamHandler, MaCh3Modes *Modes_);
    /// @brief Destructor
    virtual ~BinnedSplineHandler();

    /// @brief CW: This Eval should be used when using two separate x,{y,a,b,c,d} arrays
    /// to store the weights; probably the best one here! Same thing but pass parameter
    /// spline segments instead of variations
    void Evaluate() final;

    /// @brief add oscillation channel to spline monolith
    void AddSample(const std::string& SampleName,
                   const std::string& SampleTitle,
                   const std::vector<std::string>& OscChanFileNames,
                   const std::vector<std::string>& SplineVarNames);
    /// @brief flatten multidimensional spline array into proper monolith
    void TransferToMonolith();
    /// @brief Remove setup variables not needed for spline evaluations
    void CleanUpMemory();

    /// @brief Loads and processes splines from ROOT files for a given sample.
    /// @note DB Add virtual so it can be overridden in experiment specific (if needed)
    virtual void FillSampleArray(const std::string& SampleTitle, const std::vector<std::string>& OscChanFileNames);
    /// @brief Return the splines which affect a given event
    std::vector<SplineIndex> GetEventSplines(const std::string& SampleTitle, int iOscChan, int EventMode, const std::vector<double>& VarVals);
    /// @brief Count how many splines we have
    int CountNumberOfLoadedSplines(bool NonFlat=false, int Verbosity=0) const;

    /// @brief get pointer to spline weight based on bin variables
    const M3::float_t* RetPointer(const SplineIndex& Variables) const;
    /// @brief KS: Prepare spline file that can be used for fast loading
    void PrepareSplineFile(std::string FileName) final;
    /// @brief KS: Load preprocessed spline file
    /// @param FileName Path to ROOT file with predefined reduced Spline Monolith
    void LoadSplineFile(std::string FileName) final;

  protected:
    /// @brief CPU based code which eval weight for each spline
    void CalcSplineWeights() final;
    /// @brief Initialise flat structure
    void PrepForReweight();
    /// Only need 1 indexing array everything else interfaces with this to get binning properties
    void BuildSampleIndexingArray(const std::string& SampleTitle);
    /// @brief Grab histograms with spline binning
    std::vector<TAxis*> FindSplineBinning(const std::string& FileName, const std::string& SampleTitle);

    /// @brief Print spline binning
    void PrintBinning(TAxis* Axis) const;
    /// @brief Print info like Sample ID of spline params etc.
    void PrintSampleDetails(const std::string& SampleTitle) const;
    /// @brief Print info like Sample ID of spline params etc.
    void PrintArrayDetails(const std::string& SampleTitle) const;
    /// @brief Get index of sample based on name
    /// @param SampleTitle The title of the sample to search for.
    int GetSampleIndex(const std::string& SampleTitle) const;
    /// @brief Ensure we have spline for a given bin
    bool isValidSplineIndex(const std::string& SampleTitle, int iSyst, int iOscChan, int iMode, const std::vector<int>& iVar) const;
    /// @brief Creates an array to be filled with monolith indexes for each sample (allows for indexing between 7D binning and 1D Vector).
    /// @brief Check if there are any repeated modes. This is used to reduce the number
    /// of modes in case many interaction modes get averaged into one spline
    std::vector< std::vector<int> > StripDuplicatedModes(const std::vector< std::vector<int> >& InputVector) const;
    /// @brief Rather work with spline coefficients in the splines, let's copy ND and use coefficient arrays
    void GetSplineCoeff_SepMany(int splineindex, M3::float_t *& xArray, M3::float_t *&manyArray);

    /// Pointer to covariance from which we get information about spline params
    ParameterHandlerGeneric* ParHandler;

    //And now the actual member variables
    std::vector<std::string> SampleNames;
    std::vector<std::string> SampleTitles;
    std::vector<int> Dimensions;
    std::vector<std::vector<std::string>> DimensionLabels;
    std::vector<int> nSplineParams;
    std::vector<int> nOscChans;
    /// Holds TAxis for [sample][channel][dimension]
    std::vector< std::vector< std::vector<TAxis*> > > SplineBinning;
    /// [Sample][Syst]
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
    /// a call to GetSplineInterpolationFromSampleID()
    std::vector< std::vector<SplineInterpolation> > SplineInterpolationTypes;

    /// name of each spline parameter
    std::vector<std::string> UniqueSystNames;
    /// Global index of each spline param, it allows us to match spline ordering with global
    std::vector<int> UniqueSystIndices;

    /// @brief Variables related to determined which modes have splines and which piggy-back of other modes
    std::vector<SplineIndex> IndexVect;
    /// @brief Map between spline origin/properties (iSample, iOscChan, iSyst, iMode, iVar1, iVar2, iVar3) and the index of the spline in IndexVect
    std::map<std::tuple<int, int, int, int, std::vector<int>>, int> IndexVectMap;
    /// Number of coefficients for a single flat (after flattening)
    std::vector<unsigned int> coeffindexvec;
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
    std::vector<M3::float_t> cpu_spline_weights;
    /// Maps single spline object with single parameter
    std::vector<short int> uniquesplinevec_Monolith;

    /// pointer to MaCh3 Mode from which we get spline suffix
    MaCh3Modes* Modes;
    enum TokenOrdering{kSystToken,kModeToken,kVarBinToken};
    /// @brief Extract metadata tokens encoded in a spline name.
    /// allows experiment to have different formats of splines
    virtual std::vector<std::string> GetTokensFromSplineName(const std::string& FullSplineName) = 0;

  private:
    /// @brief CW: The shared initialiser from constructors of TResponseFunction_red
    void MoveToGPU();
    /// @brief This function will find missing splines in file
    void InvestigateMissingSplines() const;
    /// @brief KS: Load preprocessed Settings
    /// @param SplineFile File from which we load new tree
    void LoadSettingsDir(std::unique_ptr<TFile>& SplineFile);
    /// @brief KS: Load preprocessed Monolith
    /// @param SplineFile File from which we load new tree
    void LoadMonolithDir(std::unique_ptr<TFile>& SplineFile);
    /// @brief KS: Load preprocessed Index
    /// @param SplineFile File from which we load new tree
    void LoadIndexDir(std::unique_ptr<TFile>& SplineFile);

    /// @brief KS: Prepare Settings Info within SplineFile
    /// @param SplineFile File from which we load new tree
    void PrepareSettingsDir(std::unique_ptr<TFile>& SplineFile) const;
    /// @brief KS: Prepare Monolith Info within SplineFile
    /// @param SplineFile File from which we load new tree
    void PrepareMonolithDir(std::unique_ptr<TFile>& SplineFile) const;
    /// @brief KS: Prepare Index Info within SplineFile
    /// @param SplineFile File from which we load new tree
    void PrepareIndexDir(std::unique_ptr<TFile>& SplineFile) const;
    /// @brief KS: Prepare Other Info within SplineFile
    /// @param SplineFile File from which we load new tree
    void PrepareOtherInfoDir(std::unique_ptr<TFile>& SplineFile) const;
};
