#pragma once

//MaCh3 includes
#include "Splines/BinnedSplineHandler.h"
#include "Parameters/ParameterHandlerGeneric.h"
#include "Samples/SampleHandlerInterface.h"
#include "Samples/OscillationHandler.h"
#include "Samples/FarDetectorCoreInfoStruct.h"
#include "Samples/BinningHandler.h"

_MaCh3_Safe_Include_Start_ //{
#include "THStack.h"
#include "TLegend.h"
_MaCh3_Safe_Include_End_ //}

/// @brief Class responsible for handling implementation of samples used in analysis, reweighting and returning LLH
/// @author Dan Barrow
/// @author Ed Atkin
///
/// @ingroup SamplesAndParameters
class SampleHandlerBase :  public SampleHandlerInterface
{
 public:
  //######################################### Functions #########################################
  /// @brief Constructor
  /// @param ConfigFileName Name of config to initialise the sample object
  SampleHandlerBase(std::string ConfigFileName, ParameterHandlerGeneric* xsec_cov,
                  const std::shared_ptr<OscillationHandler>& OscillatorObj_ = nullptr);
  /// @brief destructor
  virtual ~SampleHandlerBase();

  /// @brief DB Get what dimensionality binning for given sample has
  /// @param Sample Number of sample
  int GetNDim(const int Sample) const final { return SampleDetails[Sample].nDimensions; }
  /// @brief Get name for Sample Handler
  std::string GetName() const final;
  /// @brief Get fancy title for specified samples
  std::string GetSampleTitle(const int Sample) const final {return SampleDetails[Sample].SampleTitle;}

  /// @brief Return Kinematic Variable name for specified sample and dimension for example "Reconstructed_Neutrino_Energy"
  /// @param iSample Sample index
  /// @param Dimension Dimension index
  std::string GetKinVarName(const int iSample, const int Dimension) const final;

  /// @brief Computes and prints the integral breakdown of all modes and oscillation channels for a given sample.
  void PrintIntegral(const int iSample, const TString& OutputName="/dev/null", const int WeightStyle=0, const TString& OutputCSVName="/dev/null");

  //===============================================================================
  // DB Reweighting and Likelihood functions

  //ETA - abstract these to SampleHandlerBase
  //DB Require these four functions to allow conversion from TH1(2)D to array for multi-threaded GetLikelihood
  void AddData(const int Sample, TH1* Data);
  void AddData(const int Sample, const std::vector<double>& Data_Array);

  /// @brief Helper function to print rates for the samples with LLH
  /// @param DataOnly whether to print data only rates
  void PrintRates(const bool DataOnly = false) final;
  /// @brief DB Multi-threaded GetLikelihood
  double GetLikelihood() const override;
  /// @brief Get likelihood for single sample
  double GetSampleLikelihood(const int isample) const override;
  //===============================================================================

  /// @brief Get index of sample based on name
  /// @param SampleTitle The title of the sample to search for.
  int GetSampleIndex(const std::string& SampleTitle) const;

  /// @brief Get Data histogram
  const TH1* GetDataHist(const int Sample) final;
  const TH1* GetDataHist(const std::string& Sample);

  /// @brief Get MC histogram
  const TH1* GetMCHist(const int Sample) final;
  const TH1* GetMCHist(const std::string& Sample);

  /// @brief Get W2 histogram
  const TH1* GetW2Hist(const int Sample) final;
  const TH1* GetW2Hist(const std::string& Sample);
  /// @brief main routine modifying MC prediction based on proposed parameter values
  void Reweight() override;
  /// @brief Computes the total event weight for a given entry.
  M3::float_t GetEventWeight(const int iEntry);

  const M3::float_t* GetNuOscillatorPointers(const int iEvent) const;

  /// @brief Get number of oscillation channels for a single sample
  int GetNOscChannels(const int iSample) const final {return static_cast<int>(SampleDetails[iSample].OscChannels.size());};

  std::string GetFlavourName(const int iSample, const int iChannel) const final {
    if (iChannel < 0 || iChannel > GetNOscChannels(iSample)) {
      MACH3LOG_ERROR("Invalid Channel Requested: {}", iChannel);
      throw MaCh3Exception(__FILE__ , __LINE__);
    }
    return SampleDetails[iSample].OscChannels[iChannel].flavourName;
  }
  std::unique_ptr<TH1> Get1DVarHist(const int iSample, const std::string &ProjectionVar,
                                    const std::vector<KinematicCut> &EventSelectionVec = {}, int WeightStyle = 0,
                                    const std::vector<KinematicCut> &SubEventSelectionVec = {}) final;
  std::unique_ptr<TH2> Get2DVarHist(const int iSample, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                                    const std::vector< KinematicCut >& EventSelectionVec = {},
                                    int WeightStyle = 0, const std::vector< KinematicCut >& SubEventSelectionVec = {}) final;
  std::vector<KinematicCut> BuildModeChannelSelection(const int iSample, const int kModeToFill, const int kChannelToFill) const;

  void Fill1DSubEventHist(const int iSample, TH1D* _h1DVar, const std::string& ProjectionVar,
                          const std::vector< KinematicCut >& SubEventSelectionVec = {},
                          int WeightStyle=0);
  void Fill2DSubEventHist(const int iSample, TH2* _h2DVar, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                          const std::vector< KinematicCut >& SubEventSelectionVec = {}, int WeightStyle = 0);

  std::unique_ptr<TH1> Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str,
                                                    const int kModeToFill = -1, const int kChannelToFill = -1,
                                                    const int WeightStyle = 0) final;
  std::unique_ptr<TH2> Get2DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_StrX,
                                                    const std::string& ProjectionVar_StrY, const int kModeToFill = -1,
                                                    const int kChannelToFill = -1, const int WeightStyle = 0) final;

  std::unique_ptr<TH1> GetModeHist1D(const int iSample, int s, int m, int style = 0) {
    return Get1DVarHistByModeAndChannel(iSample, GetKinVarName(iSample, 0), m, s, style);
  }
  std::unique_ptr<TH2> GetModeHist2D(const int iSample, int s, int m, int style = 0) {
    return Get2DVarHistByModeAndChannel(iSample, GetKinVarName(iSample, 0), GetKinVarName(iSample, 1), m, s, style);
  }

  std::vector<std::unique_ptr<TH1>> ReturnHistsBySelection1D(const int iSample, const std::string& KinematicProjection,
                                                             const int Selection1, const int Selection2 = -1,
                                                             const int WeightStyle = 0);
  std::vector<std::unique_ptr<TH2>> ReturnHistsBySelection2D(const int iSample, const std::string& KinematicProjectionX,
                                                             const std::string& KinematicProjectionY,
                                                             const int Selection1, const int Selection2=-1,
                                                             const int WeightStyle=0);
  std::unique_ptr<THStack> ReturnStackedHistBySelection1D(const int iSample, const std::string& KinematicProjection,
                                          const int Selection1, const int Selection2 = -1, const int WeightStyle = 0);
  /// @brief Return the legend used for stacked histograms with sample info
  const TLegend* ReturnStackHistLegend() const {return THStackLeg;}

  /// @brief ETA function to generically convert a string from xsec cov to a kinematic type
  int ReturnKinematicParameterFromString(const std::string& KinematicStr) const;
  /// @brief ETA function to generically convert a kinematic type from xsec cov to a string
  std::string ReturnStringFromKinematicParameter(const int KinematicVariable) const;

  /// @brief Store additional info in a chan
  void SaveAdditionalInfo(TDirectory* Dir) final;

  /// @brief JM: Convert a kinematic vector name to its corresponding integer ID.
  int ReturnKinematicVectorFromString(const std::string& KinematicStr) const;
  /// @brief JM: Convert a kinematic vector integer ID to its corresponding name as a string.
  std::string ReturnStringFromKinematicVector(const int KinematicVariable) const;
  /// @brief JM: Check if a kinematic parameter string corresponds to a subevent-level variable
  bool IsSubEventVarString(const std::string& VarStr) const;

  /// @brief Return array storing data entries for every bin
  auto GetDataArray() const {
    return SampleHandler_data;
  }
  /// @brief Return array storing MC entries for every bin
  auto GetMCArray() const {
    return SampleHandler_array;
  }
  /// @brief Return array storing W2 entries for every bin
  auto GetW2Array() const {
    return SampleHandler_array_w2;
  }
  /// @brief Return a sub-array for a given sample.
  std::vector<double> GetArrayForSample(const int Sample, std::vector<double> const & array) const;

  /// @brief Return array storing data entries for every bin
  std::vector<double> GetDataArray(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandler_data);
  }
  /// @brief Return array storing MC entries for every bin
  std::vector<double> GetMCArray(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandler_array);
  }
  /// @brief Return array storing W2 entries for single sample
  std::vector<double> GetW2Array(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandler_array_w2);
  }

 protected:
  /// @brief including Dan's magic NuOscillator
  void InitialiseNuOscillatorObjects();
  /// @brief Initialise pointer to oscillation weight to NuOscillator object
  void SetupNuOscillatorPointers();
  /// @brief Load information about sample handler and corresponding samples from config file
  void ReadConfig();
  /// @brief Initialise single sample from config file
  void LoadSingleSample(const int iSample, const YAML::Node& Settings);

  /// @brief DB Function to determine which weights apply to which types of samples
  virtual void AddAdditionalWeightPointers() = 0;

  /// @brief Ensure Kinematic Map is setup and make sure it is initialised correctly
  void SetupKinematicMap();

  /// @todo abstract the spline initialisation completely to core
  /// @brief initialise your splineXX object and then use InitialiseSplineObject to conviently setup everything up
  virtual void SetupSplines() = 0;

  //DB Require all objects to have a function which reads in the MC
  /// @brief Initialise any variables that your experiment specific SampleHandler needs
  virtual void Init() = 0;

  /// @brief Experiment specific setup, returns the number of events which were loaded
  virtual int SetupExperimentMC() = 0;

  /// @brief Function which translates experiment struct into core struct
  virtual void SetupMC() = 0;
  /// @brief Function responsible for loading data from file or loading from file
  virtual void InititialiseData() = 0;

  /// @brief Function which does a lot of the lifting regarding the workflow in creating different MC objects
  void Initialise();

  /// @brief Contains all your splines (binned or unbinned) and handles the setup and the returning of weights from spline evaluations
  std::unique_ptr<SplineBase> SplineHandler;

  /// @brief Contains oscillator handling calculating oscillation probabilities
  std::shared_ptr<OscillationHandler> Oscillator;
  //===============================================================================
  /// @brief Set pointers for each event to appropriate weights, for unbinned based on event number
  /// while for binned based on other kinematical properties
  void SetSplinePointers();
  /// @brief Retrieve the spline bin indices associated with a given event.
  /// @warning ThrowCrititcal argument will be eventually removed
  std::vector< std::vector<int> > GetSplineBins(int Event, BinnedSplineHandler* BinnedSpline, bool& ThrowCrititcal) const;

  //Functions which find the nominal bin and bin edges
  void FindNominalBinAndEdges();

  /// @brief set the binning for 2D sample used for the likelihood calculation
  void SetBinning();

  /// @brief Initialise data, MC and W2 histograms
  void SetupReweightArrays();
  //===============================================================================

  // ----- Functional Parameters -----
  /// @brief ETA - a function to setup and pass values to functional parameters where you need to pass a value to some custom reweight calc or engine
  virtual void SetupFunctionalParameters();
  /// @brief HH - a helper function for RegisterFunctionalParameter
  void RegisterIndividualFunctionalParameter(const std::string& fpName, int fpEnum, FuncParFuncType fpFunc);
  /// @brief HH - a experiment-specific function where the maps to actual functions are set up
  virtual void RegisterFunctionalParameters() = 0;
  /// @brief Update the functional parameter values to the latest proposed values. Needs to be called before every new reweight so is called in fillArray
  virtual void PrepFunctionalParameters(){};
  /// @brief ETA - generic function applying shifts
  virtual void ApplyShifts(const int iEvent);

  /// @brief DB Function which determines if an event is selected based on @ref KinematicCut
  bool IsEventSelected(const int iSample, const int iEvent) _noexcept_;
  /// @brief JM Function which determines if a subevent is selected
  bool IsSubEventSelected(const std::vector<KinematicCut> &SubEventCuts, const int iEvent, unsigned const int iSubEvent, size_t nsubevents);
  /// @brief HH - reset the shifted values to the original values
  virtual void ResetShifts(const int iEvent) {(void)iEvent;};
  /// @brief LP - Optionally calculate derived observables after all shifts have been applied
  /// @details LP - For example, have shifts that varied lepton energy and hadron energy separately
  ///               in a subclass implementation of this method you may add the shifted quantities
  ///               together to build a shifted neutrino energy estimator
  virtual void FinaliseShifts(const int iEvent) {(void)iEvent;};
  /// @brief HH - a grid of vectors of enums for each sample and event
  std::vector<std::vector<FunctionalShifter*>> funcParsGrid;
  /// @brief HH - a map that relates the funcpar enum to pointer of FuncPars
  /// struct
  /// HH - Changed to a vector of pointers since it's faster than unordered_map
  /// and we are using ints as keys
  std::vector<FunctionalShifter> funcParsMap;

  /// @todo KS: Below functional variables are used only on setup, thus we should refactor them in such a way
  /// that they are removed as class members but this would be breaking change thus keep it for the time being.

  /// @brief HH - a vector that stores all the FuncPars struct
  std::vector<FunctionalParameter> funcParsVec;
  /// @brief HH - a map that relates the name of the functional parameter to
  /// funcpar enum
  std::unordered_map<std::string, int> funcParsNamesMap;
  /// @brief HH - a map that relates the funcpar enum to pointer of the actual
  /// function
  std::unordered_map<int, FuncParFuncType> funcParsFuncMap;
  /// @brief HH - a vector of string names for each functional parameter
  std::vector<std::string> funcParsNamesVec = {};

  /// @brief Check whether a normalisation systematic affects an event or not
  void CalcNormsBins(std::vector<NormParameter>& norm_parameters, std::vector< std::vector< int > >& norms_bins);
  template <typename ParT> bool PassesSelection(const ParT& Par, std::size_t iEvent);
  /// @brief Calculate the total weight weight for a given event
  M3::float_t CalcWeightTotal(const EventInfo* _restrict_ MCEvent) const _noexcept_;

  /// @brief Calculate weights for function parameters
  ///
  /// First you need to setup additional pointers in you experiment code in SetupWeightPointers
  /// Then in this function you can calculate whatever fancy function you want by filling weight to which you have pointer
  /// This way func weight shall be used in GetEventWeight
  virtual void CalcWeightFunc(const int iEvent) {return; (void)iEvent;};

  /// @brief Return the value of an associated kinematic parameter for an event
  double ReturnKinematicParameter(const std::string& KinematicParameter, int iEvent) const {
    return ReturnKinematicParameter(ReturnKinematicParameterFromString(KinematicParameter), iEvent);
  }
  virtual double ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const = 0;

  // === JM declare the same functions for kinematic vectors ===
  std::vector<double> ReturnKinematicVector(const std::string& KinematicParameter, const int iEvent) const {
    return ReturnKinematicVector(ReturnKinematicVectorFromString(KinematicParameter), iEvent);
  }
  virtual std::vector<double> ReturnKinematicVector(const int KinematicVariable, const int iEvent) const {
    return {}; (void)KinematicVariable; (void)iEvent;};
  // ===========================================================

  /// @brief Return the binning used to draw a kinematic parameter
  std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string &KinematicParameter) const final;

  const double* GetPointerToKinematicParameter(const std::string& KinematicParameter, int iEvent) const {
    return GetPointerToKinematicParameter(ReturnKinematicParameterFromString(KinematicParameter), iEvent);
  }
  virtual const double* GetPointerToKinematicParameter(const int KinematicVariable, const int iEvent) const = 0;

  /// @brief Get pointer to oscillation channel associated with given event. Osc channel is const
  const double* GetPointerToOscChannel(const int iEvent) const;
  /// @brief Setup the norm parameters by assigning each event with bin
  void SetupNormParameters();
  /// @brief Setup the osc parameters
  void SetupOscParameters();
  //===============================================================================
  /// @brief Fill a histogram with the event-level information used in the fit
  /// @details
  /// DB Functions required for reweighting functions
  /// DB Replace previous implementation with reading bin contents from SampleHandler_array
  void FillHist(const int Sample, TH1* Hist, std::vector<double> &Array);

  /// @brief DB Nice new multi-threaded function which calculates the event weights and fills the relevant bins of an array
#ifdef MULTITHREAD
  /// @brief Function which does the core reweighting, fills the @ref SampleHandlerBase::SampleHandler_array
  /// vector with the weight calculated from reweighting but multithreaded
  void FillArray_MP();
#endif
  /// @brief Function which does the core reweighting, fills the @ref SampleHandlerBase::SampleHandler_array
  /// vector with the weight calculated from reweighting
  void FillArray();

  /// @brief Helper function to reset histograms
  void ResetHistograms();

  //===============================================================================
  //DB Variables required for GetLikelihood
  /// KS: This stores binning information, in future could be come vector to store binning for every used sample
  std::unique_ptr<BinningHandler> Binning;
  /// DB Array to be filled after reweighting
  std::vector<double> SampleHandler_array;
  /// KS Array used for MC stat
  std::vector<double> SampleHandler_array_w2;
  /// DB Array to be filled in AddData
  std::vector<double> SampleHandler_data;
  //===============================================================================

  //===============================================================================
  /// Stores information about every MC event
  std::vector<EventInfo> MCEvents;
  /// Stores info about currently initialised sample
  std::vector<SampleInfo> SampleDetails;
  //===============================================================================

  //===============================================================================
  //DB Covariance Objects
  /// ETA - All experiments will need an xsec, det and osc cov
  ParameterHandlerGeneric *ParHandler = nullptr;

  //===============================================================================
  /// @brief A unique ID for each sample based on which we can define what systematic should be applied
  std::string SampleHandlerName;

  //===========================================================================
  //DB Vectors to store which kinematic cuts we apply
  //like in XsecNorms but for events in sample. Read in from sample yaml file
  //What gets used in IsEventSelected, which gets set equal to user input plus
  //all the vectors in StoreSelection

  /// @brief What gets pulled from config options, these are constant after loading in
  /// this is of length 3: 0th index is the value, 1st is lower bound, 2nd is upper bound
  std::vector< std::vector< KinematicCut > > StoredSelection;
  /// @brief a way to store selection cuts which you may push back in the get1DVar functions
  /// most of the time this is just the same as StoredSelection
  std::vector< std::vector< KinematicCut > > Selection;
   //===========================================================================

  /// Mapping between string and kinematic enum
  const std::unordered_map<std::string, int>* KinematicParameters;
  /// Mapping between kinematic enum and string
  const std::unordered_map<int, std::string>* ReversedKinematicParameters;

  // === JM mapping between string and kinematic vector enum ===
  const std::unordered_map<std::string, int>* KinematicVectors;
  const std::unordered_map<int, std::string>* ReversedKinematicVectors;
  // ===========================================================

  /// The manager object used to read the sample yaml file
  std::unique_ptr<Manager> SampleManager;
  void InitialiseSplineObject();

  std::unordered_map<std::string, double> _modeNomWeightMap;

  //===============================================================================
  /// DB Miscellaneous Variables
  TLegend* THStackLeg = nullptr;
  //===============================================================================

  /// KS:Super hacky to update W2 or not
  bool FirstTimeW2;
  /// KS:Super hacky to update W2 or not
  bool UpdateW2;

  /// @brief Retrieve the initial neutrino PDG code associated with a given input file name.
  NuPDG GetInitPDGFromFileName(const std::string& FileName) const {return FileToInitPDGMap.at(FileName);}
  /// @brief Retrieve the final neutrino PDG code associated with a given input file name.
  NuPDG GetFinalPDGFromFileName(const std::string& FileName) const {return FileToFinalPDGMap.at(FileName);}

 private:
   /// @brief Temporarily extend Selection for a given sample with additional cuts.
   /// Returns the original Selection so the caller can restore it later.
   std::vector<std::vector<KinematicCut>> ApplyTemporarySelection(const int iSample,
                                                                  const std::vector<KinematicCut>& ExtraCuts);

  std::unordered_map<std::string, NuPDG> FileToInitPDGMap;
  std::unordered_map<std::string, NuPDG> FileToFinalPDGMap;

  enum FDPlotType {
    kModePlot = 0,
    kOscChannelPlot = 1
  };
};
