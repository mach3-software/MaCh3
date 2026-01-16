#pragma once

//MaCh3 includes
#include "Splines/BinnedSplineHandler.h"

#include "Parameters/ParameterHandlerGeneric.h"

#include "Samples/SampleHandlerBase.h"
#include "Samples/OscillationHandler.h"
#include "Samples/FarDetectorCoreInfoStruct.h"
#include "Samples/BinningHandler.h"

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
                  const std::shared_ptr<OscillationHandler>& OscillatorObj_ = nullptr);
  /// @brief destructor
  virtual ~SampleHandlerFD();

  /// @brief DB Function to differentiate 1D or 2D binning
  /// @ingroup SampleHandlerGetters
  int GetNDim(const int Sample) const { return SampleDetails[Sample].nDimensions; }
  /// @ingroup SampleHandlerGetters
  std::string GetName() const override;
  /// @ingroup SampleHandlerGetters
  std::string GetSampleTitle(const int Sample) const override {return SampleDetails[Sample].SampleTitle;}

  /// @ingroup SampleHandlerGetters
  std::string GetXBinVarName(const int Sample) const {return SampleDetails[Sample].VarStr[0];}
  /// @ingroup SampleHandlerGetters
  std::string GetYBinVarName(const int Sample) const {return SampleDetails[Sample].VarStr[1];}
  /// @brief Get pointer to binning handler
  const BinningHandler* GetBinningHandler() const {return Binning.get();}

  void PrintIntegral(const int iSample, const TString& OutputName="/dev/null", const int WeightStyle=0, const TString& OutputCSVName="/dev/null");
  
  //===============================================================================
  // DB Reweighting and Likelihood functions

  //ETA - abstract these to SampleHandlerFDBase
  //DB Require these four functions to allow conversion from TH1(2)D to array for multi-threaded GetLikelihood
  void AddData(const int Sample, TH1D* Data);
  void AddData(const int Sample, TH2D* Data);
  void AddData(const int Sample, std::vector<double> &data);
  void AddData(const int Sample, std::vector< std::vector <double> > &data);


  /// @brief Helper function to print rates for the samples with LLH
  /// @param DataOnly whether to print data only rates
  void PrintRates(const bool DataOnly = false) override;
  /// @brief DB Multi-threaded GetLikelihood
  double GetLikelihood() const override;
  /// @ingroup SampleHandlerGetters
  /// @brief Get likelihood for single sample
  double GetSampleLikelihood(const int isample) const override;
  //===============================================================================

  /// @brief Get index of sample based on name
  /// @param SampleTitle The title of the sample to search for.
  /// @ingroup SampleHandlerGetters
  int GetSampleIndex(const std::string& SampleTitle) const;

  /// @brief Get MC histogram
  /// @ingroup SampleHandlerGetters
  TH1* GetMCHist(const int Sample, const int Dimension);
  TH1* GetMCHist(const std::string& Sample, const int Dimension);

  /// @brief Get W2 histogram
  /// @ingroup SampleHandlerGetters
  TH1* GetW2Hist(const int Sample, const int Dimension);
  TH1* GetW2Hist(const std::string& Sample, const int Dimension);

  /// @brief Get Data histogram
  /// @ingroup SampleHandlerGetters
  TH1* GetDataHist(const int Sample, const int Dimension);
  TH1* GetDataHist(const std::string& Sample, const int Dimension);

  void Reweight() override;
  M3::float_t GetEventWeight(const int iEntry);

  /// @brief including Dan's magic NuOscillator
  void InitialiseNuOscillatorObjects();
  void SetupNuOscillatorPointers();
  const M3::float_t* GetNuOscillatorPointers(const int iEvent) const;

  void ReadConfig();
  void LoadSingleSample(const int iSample, const YAML::Node& Settings);

  /// @ingroup SampleHandlerGetters
  int GetNOscChannels(const int iSample) const override {return static_cast<int>(SampleDetails[iSample].OscChannels.size());};

  /// @ingroup SampleHandlerGetters
  std::string GetFlavourName(const int iSample, const int iChannel) {
    if (iChannel < 0 || iChannel > GetNOscChannels(iSample)) {
      MACH3LOG_ERROR("Invalid Channel Requested: {}", iChannel);
      throw MaCh3Exception(__FILE__ , __LINE__);      
    }
    return SampleDetails[iSample].OscChannels[iChannel].flavourName;
  }

  /// @ingroup SampleHandlerGetters
  TH1* Get1DVarHist(const int iSample, const std::string& ProjectionVar, const std::vector< KinematicCut >& EventSelectionVec = std::vector< KinematicCut >(),
                    int WeightStyle=0, TAxis* Axis=nullptr, const std::vector< KinematicCut >& SubEventSelectionVec = std::vector< KinematicCut >());
  TH2* Get2DVarHist(const int iSample, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                    const std::vector< KinematicCut >& EventSelectionVec = std::vector< KinematicCut >(),
                    int WeightStyle=0, TAxis* AxisX=nullptr, TAxis* AxisY=nullptr,
                    const std::vector< KinematicCut >& SubEventSelectionVec = std::vector< KinematicCut >());
  
  /// @ingroup SampleHandlerGetters
  void Fill1DSubEventHist(const int iSample, TH1D* _h1DVar, const std::string& ProjectionVar, const std::vector< KinematicCut >& SubEventSelectionVec = std::vector< KinematicCut >(),
                    int WeightStyle=0);
  void Fill2DSubEventHist(const int iSample, TH2D* _h2DVar, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                    const std::vector< KinematicCut >& SubEventSelectionVec = std::vector< KinematicCut >(), int WeightStyle=0);
  
  /// @ingroup SampleHandlerGetters
  TH1* Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str, int kModeToFill=-1, int kChannelToFill=-1, int WeightStyle=0, TAxis* Axis=nullptr);
  /// @ingroup SampleHandlerGetters
  TH2* Get2DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_StrX, const std::string& ProjectionVar_StrY, int kModeToFill=-1, int kChannelToFill=-1, int WeightStyle=0, TAxis* AxisX=nullptr, TAxis* AxisY=nullptr);

  /// @ingroup SampleHandlerGetters
  TH1 *GetModeHist1D(const int iSample, int s, int m, int style = 0) {
    return Get1DVarHistByModeAndChannel(iSample, GetXBinVarName(iSample), m, s, style);
  }
  /// @ingroup SampleHandlerGetters
  TH2 *GetModeHist2D(const int iSample, int s, int m, int style = 0) {
    return Get2DVarHistByModeAndChannel(iSample, GetXBinVarName(iSample),GetYBinVarName(iSample), m, s, style);
  }

  /// @ingroup SampleHandlerGetters
  std::vector<TH1*> ReturnHistsBySelection1D(const int iSample, std::string KinematicProjection, int Selection1,int Selection2=-1, int WeightStyle=0, TAxis* Axis=0);
  /// @ingroup SampleHandlerGetters
  std::vector<TH2*> ReturnHistsBySelection2D(const int iSample, std::string KinematicProjectionX, std::string KinematicProjectionY, int Selection1, int Selection2=-1, int WeightStyle=0, TAxis* XAxis=0, TAxis* YAxis=0);
  /// @ingroup SampleHandlerGetters
  THStack* ReturnStackedHistBySelection1D(const int iSample, const std::string& KinematicProjection, int Selection1, int Selection2=-1, int WeightStyle=0, TAxis* Axis=0);
  /// @ingroup SampleHandlerGetters
  TLegend* ReturnStackHistLegend() {return THStackLeg;}
  
  /// @brief ETA function to generically convert a string from xsec cov to a kinematic type
  /// @ingroup SampleHandlerGetters
  int ReturnKinematicParameterFromString(const std::string& KinematicStr) const;
  /// @brief ETA function to generically convert a kinematic type from xsec cov to a string
  /// @ingroup SampleHandlerGetters
  std::string ReturnStringFromKinematicParameter(const int KinematicVariable) const;

  /// @brief Store additional info in a chan
  void SaveAdditionalInfo(TDirectory* Dir) override;

  // === JM declare the same functions for kinematic vectors ===
  int ReturnKinematicVectorFromString(const std::string& KinematicStr) const;
  std::string ReturnStringFromKinematicVector(const int KinematicVariable) const;
  // ===========================================================
  /// @brief JM Check if a kinematic parameter string corresponds to a subevent-level variable
  bool IsSubEventVarString(const std::string& VarStr);

 protected:
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
  virtual void SetupFDMC() = 0;

  /// @brief Function which does a lot of the lifting regarding the workflow in creating different MC objects
  void Initialise();
  
  /// @brief Contains all your binned splines and handles the setup and the returning of weights from spline evaluations
  std::unique_ptr<BinnedSplineHandler> SplineHandler;

  /// @brief Contains all your binned splines and handles the setup and the returning of weights from spline evaluations
  std::shared_ptr<OscillationHandler> Oscillator;
  //===============================================================================
  /// @brief Finds the binned spline that an event should apply to and stored them in a
  /// a vector for easy evaluation in the fillArray() function.
  void FillSplineBins();

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
  virtual void ApplyShifts(int iEvent);

  /// @brief DB Function which determines if an event is selected, where Selection double looks like {{ND280KinematicTypes Var1, douuble LowBound}
  bool IsEventSelected(const int iSample, const int iEvent);
  /// @brief JM Function which determines if a subevent is selected
  bool IsSubEventSelected(const std::vector<KinematicCut> &SubEventCuts, const int iEvent, unsigned const int iSubEvent, size_t nsubevents);
  /// @brief HH - reset the shifted values to the original values
  virtual void resetShifts(int iEvent) {(void)iEvent;};
  /// @brief HH - a vector that stores all the FuncPars struct
  std::vector<FunctionalParameter> funcParsVec;
  /// @brief HH - a map that relates the name of the functional parameter to
  /// funcpar enum
  std::unordered_map<std::string, int> funcParsNamesMap;
  /// @brief HH - a map that relates the funcpar enum to pointer of FuncPars
  /// struct
  /// HH - Changed to a vector of pointers since it's faster than unordered_map
  /// and we are using ints as keys
  std::vector<FunctionalParameter *> funcParsMap;
  /// @brief HH - a map that relates the funcpar enum to pointer of the actual
  /// function
  std::unordered_map<int, FuncParFuncType> funcParsFuncMap;
  /// @brief HH - a grid of vectors of enums for each sample and event
  std::vector<std::vector<int>> funcParsGrid;
  /// @brief HH - a vector of string names for each functional parameter
  std::vector<std::string> funcParsNamesVec = {};

  /// @brief Check whether a normalisation systematic affects an event or not
  void CalcNormsBins(std::vector<NormParameter>& norm_parameters, std::vector< std::vector< int > >& xsec_norms_bins);
  /// @brief Calculate the total weight weight for a given event
  M3::float_t CalcWeightTotal(const FarDetectorCoreInfo* _restrict_ MCEvent) const;

  /// @brief Calculate weights for function parameters
  ///
  /// First you need to setup additional pointers in you experiment code in SetupWeightPointers
  /// Then in this function you can calculate whatever fancy function you want by filling weight to which you have pointer
  /// This way func weight shall be used in GetEventWeight
  virtual void CalcWeightFunc(int iEvent) {return; (void)iEvent;};

  /// @brief Return the value of an associated kinematic parameter for an event
  virtual double ReturnKinematicParameter(std::string KinematicParamter, int iEvent) = 0;
  virtual double ReturnKinematicParameter(int KinematicVariable, int iEvent) = 0;
  
  // === JM declare the same functions for kinematic vectors ===
  virtual std::vector<double> ReturnKinematicVector(std::string KinematicParameter, int iEvent) {return {}; (void)KinematicParameter; (void)iEvent;};
  virtual std::vector<double> ReturnKinematicVector(int KinematicVariable, int iEvent) {return {}; (void)KinematicVariable; (void)iEvent;};
  // ===========================================================

  /// @brief Return the binning used to draw a kinematic parameter
  std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string& KinematicParameter) const;
  virtual const double* GetPointerToKinematicParameter(std::string KinematicParamter, int iEvent) = 0;
  virtual const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent) = 0;

  /// @brief Get pointer to oscillation channel associated with given event. Osc channel is const
  const double* GetPointerToOscChannel(const int iEvent) const;
  /// @brief Setup the norm parameters by assigning each event with bin
  void SetupNormParameters();

  //===============================================================================
  //DB Functions required for reweighting functions
  //DB Replace previous implementation with reading bin contents from SampleHandlerFD_array
  /// @brief Fill a histogram with the event-level information used in the fit
  void FillMCHist(const int Sample, const int Dimension);

  /// @brief DB Nice new multi-threaded function which calculates the event weights and fills the relevant bins of an array
#ifdef MULTITHREAD
  /// @brief fills the SampleHandlerFD_array vector with the weight calculated from reweighting but multithreaded
  void FillArray_MP();
#endif
  /// @brief Function which does the core reweighting. This assumes that oscillation weights have
  /// already been calculated and stored in SampleHandlerFD.osc_w[iEvent]. This
  /// function takes advantage of most of the things called in setupSKMC to reduce reweighting time.
  /// It also follows the ND code reweighting pretty closely. This function fills the SampleHandlerFD
  /// array array which is binned to match the sample binning, such that bin[1][1] is the
  /// equivalent of SampleDetails._hPDF2D->GetBinContent(2,2) {Noticing the offset}
  void FillArray();

  /// @brief Helper function to reset histograms
  inline void ResetHistograms();
  
  //===============================================================================
  //DB Variables required for GetLikelihood
  /// KS: This stores binning information, in future could be come vector to store binning for every used sample
  std::unique_ptr<BinningHandler> Binning;
  /// DB Array to be filled after reweighting
  double* SampleHandlerFD_array;
  /// KS Array used for MC stat
  double* SampleHandlerFD_array_w2;
  /// DB Array to be filled in AddData
  double* SampleHandlerFD_data;
  //===============================================================================

  //===============================================================================
  /// Stores information about every MC event
  std::vector<FarDetectorCoreInfo> MCSamples;
  /// Stores info about currently initialised sample
  std::vector<SampleInfo> SampleDetails;
  //===============================================================================

  //===============================================================================
  //DB Covariance Objects
  //ETA - All experiments will need an xsec, det and osc cov
  //these should be added to SampleHandlerBase to be honest
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
  std::unique_ptr<manager> SampleManager;
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
  std::unordered_map<std::string, NuPDG> FileToInitPDGMap;
  std::unordered_map<std::string, NuPDG> FileToFinalPDGMap;
  
  enum FDPlotType {
    kModePlot = 0,
    kOscChannelPlot = 1
  };
};
