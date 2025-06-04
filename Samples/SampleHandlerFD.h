#pragma once

//MaCh3 includes
#include "Splines/BinnedSplineHandler.h"

#include "Parameters/ParameterHandlerGeneric.h"
#include "Parameters/ParameterHandlerOsc.h"

#include "Samples/SampleHandlerBase.h"
#include "Samples/OscillationHandler.h"
#include "Samples/FarDetectorCoreInfoStruct.h"

#include "THStack.h"
#include "TLegend.h"

/// @brief Class responsible for handling implementation of samples used in analysis, reweighting and returning LLH
/// @author Dan Barrow
/// @author Ed Atkin
class SampleHandlerFD :  public SampleHandlerBase
{
public:
  //######################################### Functions #########################################
  /// @brief Constructor
  /// @param ConfigFileName Name of config to initialise the sample object
  SampleHandlerFD(std::string ConfigFileName, ParameterHandlerGeneric* xsec_cov,
                  ParameterHandlerOsc* osc_cov = nullptr, const std::shared_ptr<OscillationHandler>& OscillatorObj_ = nullptr);
  /// @brief destructor
  virtual ~SampleHandlerFD();

  /// @ingroup SampleHandlerGetters
  int GetNDim(){return nDimensions;} //DB Function to differentiate 1D or 2D binning
  /// @ingroup SampleHandlerGetters
  std::string GetSampleName(int iSample = 0) const override;
  /// @ingroup SampleHandlerGetters
  std::string GetTitle() const override {return SampleTitle;}

  /// @ingroup SampleHandlerGetters
  std::string GetXBinVarName() {return XVarStr;}
  /// @ingroup SampleHandlerGetters
  std::string GetYBinVarName() {return YVarStr;}

  void PrintIntegral(TString OutputName="/dev/null", int WeightStyle=0, TString OutputCSVName="/dev/null");
  
  //===============================================================================
  // DB Reweighting and Likelihood functions

  //ETA - abstract these to SampleHandlerFDBase
  //DB Require these four functions to allow conversion from TH1(2)D to array for multi-threaded GetLikelihood
  void AddData(TH1D* Data) override;
  void AddData(TH2D* Data) override;
  void AddData(std::vector<double> &data) override;
  void AddData(std::vector< std::vector <double> > &data) override;
  /// @brief DB Multi-threaded GetLikelihood
  double GetLikelihood() override;
  //===============================================================================

  void Reweight() override;
  M3::float_t GetEventWeight(const int iEntry) const;

  ///  @brief including Dan's magic NuOscillator
  void InitialiseNuOscillatorObjects();
  void SetupNuOscillatorPointers();

  void ReadSampleConfig();

  /// @ingroup SampleHandlerGetters
  int GetNMCSamples() override {return nSamples;}
  /// @ingroup SampleHandlerGetters
  int GetNOscChannels() {return static_cast<int>(OscChannels.size());}

  /// @ingroup SampleHandlerGetters
  std::string GetFlavourName(const int iChannel) {
    if (iChannel < 0 || iChannel > static_cast<int>(OscChannels.size())) {
      MACH3LOG_ERROR("Invalid Channel Requested: {}", iChannel);
      throw MaCh3Exception(__FILE__ , __LINE__);      
    }
    return OscChannels[iChannel].flavourName;
  }

  /// @ingroup SampleHandlerGetters
  TH1* Get1DVarHist(const std::string& ProjectionVar, const std::vector< KinematicCut >& SelectionVec = std::vector< KinematicCut >(),
                    int WeightStyle=0, TAxis* Axis=nullptr);
  TH1* Get1DSubEventHist(const std::string& ProjectionVar, const std::vector< KinematicCut >& SelectionVec = std::vector< KinematicCut >(),
                    int WeightStyle=0, TAxis* Axis=nullptr);
  TH2* Get2DVarHist(const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                    const std::vector< KinematicCut >& SelectionVec = std::vector< KinematicCut >(),
                    int WeightStyle=0, TAxis* AxisX=nullptr, TAxis* AxisY=nullptr);
  TH2* Get2DSubEventHist(const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                    const std::vector< KinematicCut >& SelectionVec = std::vector< KinematicCut >(),
                    int WeightStyle=0, TAxis* AxisX=nullptr, TAxis* AxisY=nullptr);
  /// @ingroup SampleHandlerGetters
  TH1* Get1DVarHistByModeAndChannel(const std::string& ProjectionVar_Str, int kModeToFill=-1, int kChannelToFill=-1, int WeightStyle=0, TAxis* Axis=nullptr);
  /// @ingroup SampleHandlerGetters
  TH2* Get2DVarHistByModeAndChannel(const std::string& ProjectionVar_StrX, const std::string& ProjectionVar_StrY, int kModeToFill=-1, int kChannelToFill=-1, int WeightStyle=0, TAxis* AxisX=nullptr, TAxis* AxisY=nullptr);

  /// @ingroup SampleHandlerGetters
  TH1 *GetModeHist1D(int s, int m, int style = 0) {
    return Get1DVarHistByModeAndChannel(XVarStr,m,s,style);
  }
  /// @ingroup SampleHandlerGetters
  TH2 *GetModeHist2D(int s, int m, int style = 0) {
    return Get2DVarHistByModeAndChannel(XVarStr,YVarStr,m,s,style);
  }

  /// @ingroup SampleHandlerGetters
  std::vector<TH1*> ReturnHistsBySelection1D(std::string KinematicProjection, int Selection1,int Selection2=-1, int WeightStyle=0, TAxis* Axis=0);
  /// @ingroup SampleHandlerGetters
  std::vector<TH2*> ReturnHistsBySelection2D(std::string KinematicProjectionX, std::string KinematicProjectionY, int Selection1, int Selection2=-1, int WeightStyle=0, TAxis* XAxis=0, TAxis* YAxis=0);
  /// @ingroup SampleHandlerGetters
  THStack* ReturnStackedHistBySelection1D(std::string KinematicProjection, int Selection1, int Selection2=-1, int WeightStyle=0, TAxis* Axis=0);
  /// @ingroup SampleHandlerGetters
  TLegend* ReturnStackHistLegend() {return THStackLeg;}
  
  /// @brief ETA function to generically convert a string from xsec cov to a kinematic type
  /// @ingroup SampleHandlerGetters
  int ReturnKinematicParameterFromString(const std::string& KinematicStr) const;
  /// @brief ETA function to generically convert a kinematic type from xsec cov to a string
  /// @ingroup SampleHandlerGetters
  std::string ReturnStringFromKinematicParameter(const int KinematicVariable) const;

  // === JM declare the same functions for kinematic vectors ===
  int ReturnKinematicVectorFromString(const std::string& KinematicStr) const;
  std::string ReturnStringFromKinematicVector(const int KinematicVariable) const;
  // ===========================================================

 protected:
  /// @brief DB Function to determine which weights apply to which types of samples pure virtual!!
  virtual void SetupWeightPointers() = 0;

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
  virtual void SetupFDMC() = 0;

  /// @brief Function which does a lot of the lifting regarding the workflow in creating different MC objects
  void Initialise();
  
  /// @brief Contains all your binned splines and handles the setup and the returning of weights from spline evaluations
  std::unique_ptr<BinnedSplineHandler> SplineHandler;

  /// @brief Contains all your binned splines and handles the setup and the returning of weights from spline evaluations
  std::shared_ptr<OscillationHandler> Oscillator;
  //===============================================================================
  void FillSplineBins();

  //Functions which find the nominal bin and bin edges
  void FindNominalBinAndEdges1D();
  void FindNominalBinAndEdges2D();

  //DB Overrided functions of base class which calculate erec bin and boundaries for reweighting speedup in beam samples
  //ETA - this can be done using the core info stored in the new fdmc_struct
  /// @brief sets the binning used for the likelihood calculation, used for both data and MC
  /// @param nbins number of total bins
  /// @param boundaries the bin edges e.g. 0, 0.1, 0.2, 0.3 
  void Set1DBinning(size_t nbins, double* boundaries);
  /// @brief set the binning used for likelihood calculation using uniform binning
  /// @param nbins number of total bins
  /// @param low lower bound of the binning
  /// @param high upper bound of the binning
  void Set1DBinning(size_t nbins, double low, double high);
  /// @brief set the binning for 2D sample used for the likelihood calculation
  /// @param nbins1 number of bins in axis 1 (the x-axis)
  /// @param nbins2 number of bins in axis 2 (the y-axis)
  /// @param boundaries1 the bin boundaries used in axis 1 (the x-axis)
  /// @param boundaries2 the bin boundaries used in axis 2 (the y-axis)
  void Set2DBinning(size_t nbins1, double* boundaries1, size_t nbins2, double* boundaries2);
  /// @brief set the binning for 2D sample used for the likelihood calculation
  /// @param nbins1 number of bins in axis 1 (the x-axis)
  /// @param low1 lower bound of binning in axis 1 (the x-axis)
  /// @param high1 upper bound of binning in axis 1 (the x-axis)
  /// @param nbins2 number of bins in axis 2 (the y-axis)
  /// @param low2 lower bound of binning in axis 2 (the y-axis)
  void Set2DBinning(size_t nbins1, double low1, double high1, size_t nbins2, double low2, double high2);
  /// @brief set the binning for 1D sample used for the likelihood calculation
  /// @param XVec vector containing the binning in axis 1 (the x-axis)
  void Set1DBinning(std::vector<double> &XVec) {Set1DBinning(XVec.size()-1, XVec.data());};
  /// @brief set the binning for 2D sample used for the likelihood calculation
  /// @param XVec vector containing the binning in axis 1 (the x-axis)
  /// @param YVec vector containing the binning in axis 2 (the y-axis)
  void Set2DBinning(std::vector<double> &XVec, std::vector<double> &YVec) {Set2DBinning(XVec.size()-1, XVec.data(), YVec.size()-1, YVec.data());};
  /// @brief wrapper to call set binning functions based on sample config info
  void SetupSampleBinning();
  /// @brief Initialise data, MC and W2 histograms
  void SetupReweightArrays(const size_t numberXBins, const size_t numberYBins);

  /// @brief the strings associated with the variables used for the binning e.g. "RecoNeutrinoEnergy"
  std::string XVarStr, YVarStr;
  std::vector<std::string> SplineVarNames;
  /// @brief vector of the binning used in axis 1 (the x-axis)
  std::vector<double> SampleXBins;
  /// @brief vector of the binning used in axis 2 (the y-axis)
  std::vector<double> SampleYBins;
  //===============================================================================

  // ----- Functional Parameters -----
  /// @brief ETA - a function to setup and pass values to functional parameters where you need to pass a value to some custom reweight calc or engine
  virtual void SetupFunctionalParameters();
  /// @brief HH - a helper function for RegisterFunctionalParameter
  void RegisterIndividualFunctionalParameter(const std::string& fpName, int fpEnum, FuncParFuncType fpFunc);
  /// @brief HH - a experiment-specific function where the maps to actual functions are set up
  virtual void RegisterFunctionalParameters() = 0;
  /// @brief Update the functional parameter values to the latest propsed values. Needs to be called before every new reweight so is called in fillArray 
  virtual void PrepFunctionalParameters(){};
  /// @brief ETA - generic function applying shifts
  virtual void ApplyShifts(int iEvent);

  /// @brief DB Function which determines if an event is selected, where Selection double looks like {{ND280KinematicTypes Var1, douuble LowBound}
  bool IsEventSelected(const int iEvent);
  /// @brief JM Function which determines if a subevent is selected
  bool IsSubEventSelected(const std::vector<KinematicCut> &SubEventCuts, const int iEvent, unsigned const int iSubEvent, size_t nsubevents);
  /// @brief JM Check if a kinematic parameter string corresponds to a subevent-level variable
  bool IsSubEventVarString(const std::string& VarStr);
  /// @brief HH - reset the shifted values to the original values
  virtual void resetShifts(int iEvent) {(void)iEvent;};
  /// @brief HH - a vector that stores all the FuncPars struct
  std::vector<FunctionalParameter> funcParsVec;
  /// @brief HH - a map that relates the name of the functional parameter to
  /// funcpar enum
  std::unordered_map<std::string, int> funcParsNamesMap;
  /// @brief HH - a map that relates the funcpar enum to pointer of FuncPars
  /// struct
  // HH - Changed to a vector of pointers since it's faster than unordered_map
  // and we are using ints as keys
  std::vector<FunctionalParameter *> funcParsMap;
  /// @brief HH - a map that relates the funcpar enum to pointer of the actual
  /// function
  std::unordered_map<int, FuncParFuncType> funcParsFuncMap;
  /// @brief HH - a grid of vectors of enums for each sample and event
  std::vector<std::vector<int>> funcParsGrid;
  /// @brief HH - a vector of string names for each functional parameter
  std::vector<std::string> funcParsNamesVec = {};

  /// @brief Check whether a normalisation systematic affects an event or not
  void CalcNormsBins();
  /// @brief Calculate the spline weight for a given event
  M3::float_t CalcWeightSpline(const int iEvent) const;
  /// @brief Calculate the norm weight for a given event
  M3::float_t CalcWeightNorm(const int iEvent) const;

  /// @brief Calculate weights for function parameters
  ///
  /// First you need to setup additional pointers in you experiment code in SetupWeightPointers
  /// Then in this function you can calculate whatever fancy function you want by filling weight to which you have pointer
  /// This way func weight shall be used in GetEventWeight
  virtual void CalcWeightFunc(int iEvent){return; (void)iEvent;};

  /// @brief Return the value of an assocaited kinematic parameter for an event
  virtual double ReturnKinematicParameter(std::string KinematicParamter, int iEvent) = 0;
  virtual double ReturnKinematicParameter(int KinematicVariable, int iEvent) = 0;
  
  // === JM declare the same functions for kinematic vectors ===
  virtual std::vector<double> ReturnKinematicVector(std::string KinematicParameter, int iEvent) {return {}; (void)KinematicParameter; (void)iEvent;};
  virtual std::vector<double> ReturnKinematicVector(int KinematicVariable, int iEvent) {return {}; (void)KinematicVariable; (void)iEvent;};
  // ===========================================================

  /// @brief Return the binning used to draw a kinematic parameter
  virtual std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter) = 0;
  virtual const double* GetPointerToKinematicParameter(std::string KinematicParamter, int iEvent) = 0;
  virtual const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent) = 0;

  const double* GetPointerToOscChannel(const int iEvent) const;

  void SetupNormParameters();

  //===============================================================================
  //DB Functions required for reweighting functions
  //DB Replace previous implementation with reading bin contents from SampleHandlerFD_array
  /// @brief Fill a 1D histogram with the event-level information used in the fit
  void Fill1DHist() override;
  /// @brief Fill a 2D histogram with the event-level information used in the fit
  void Fill2DHist() override;

  /// @brief DB Nice new multi-threaded function which calculates the event weights and fills the relevant bins of an array
#ifdef MULTITHREAD
  /// @brief fills the SampleHandlerFD_array vector with the weight calculated from reweighting but multithreaded
  void FillArray_MP();
#endif
  /// @brief fills the SampleHandlerFD_array vector with the weight calculated from reweighting
  void FillArray();

  /// @brief Helper function to reset histograms
  inline void ResetHistograms();
  
  //===============================================================================
  //DB Variables required for GetLikelihood
  /// Vector to hold x-axis bin-edges
  std::vector<double> XBinEdges;
  /// Vector to hold y-axis bin-edges
  std::vector<double> YBinEdges;

  // ETA: also makes sense to store the number of X and Y bins
  /// Number of X axis bins in the histogram used for likelihood calculation
  size_t nXBins;
  /// Number of Y axis bins in the histogram used for likelihood calculation
  size_t nYBins;

  /// DB Array to be filled after reweighting
  double** SampleHandlerFD_array;
  /// KS Array used for MC stat
  double** SampleHandlerFD_array_w2;
  /// DB Array to be filled in AddData
  double** SampleHandlerFD_data;
  //===============================================================================

  //===============================================================================
  //MC variables
  FarDetectorCoreInfo MCSamples;
  std::vector<OscChannelInfo> OscChannels;
  //===============================================================================

  //===============================================================================
  //DB Covariance Objects
  //ETA - All experiments will need an xsec, det and osc cov
  //these should be added to SampleHandlerBase to be honest
  ParameterHandlerGeneric *ParHandler = nullptr;
  ParameterHandlerOsc *OscParHandler = nullptr;

  //=============================================================================== 

  /// @brief Keep track of the dimensions of the sample binning
  int nDimensions = M3::_BAD_INT_;
  /// @brief A unique ID for each sample based on powers of two for quick binary operator comparisons 
  std::string SampleName;

  /// @brief the name of this sample e.g."muon-like"
  std::string SampleTitle;

  /// @brief Information to store for normalisation pars
  std::vector<NormParameter> norm_parameters;
  //===========================================================================
  //DB Vectors to store which kinematic cuts we apply
  //like in XsecNorms but for events in sample. Read in from sample yaml file 
  //What gets used in IsEventSelected, which gets set equal to user input plus 
  //all the vectors in StoreSelection
  
  /// @brief What gets pulled from config options, these are constant after loading in
  /// this is of length 3: 0th index is the value, 1st is lower bound, 2nd is upper bound
  std::vector< KinematicCut > StoredSelection;
  /// @brief the strings grabbed from the sample config specifying the selections
  std::vector< std::string > SelectionStr;
  /// @brief a way to store selection cuts which you may push back in the get1DVar functions
  /// most of the time this is just the same as StoredSelection
  std::vector< KinematicCut > Selection;
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
  std::unique_ptr<manager> SampleManager;
  /// @brief function to create the member of the FarDetectorInfo struct so
  /// they are the appropriate size.
  void InitialiseSingleFDMCObject();
  void InitialiseSplineObject();

  std::vector<std::string> mc_files;
  std::vector<std::string> spline_files;

  std::unordered_map<std::string, double> _modeNomWeightMap;
  
  //===============================================================================
  /// DB Miscellaneous Variables
  TLegend* THStackLeg = nullptr;
  //===============================================================================

  /// KS:Super hacky to update W2 or not
  bool FirstTimeW2;
  /// KS:Super hacky to update W2 or not
  bool UpdateW2;
};
