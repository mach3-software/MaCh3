#pragma once

//MaCh3 includes
#include "splines/splineFDBase.h"

#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"

#include "samplePDF/samplePDFBase.h"
#include "samplePDF/FarDetectorCoreInfoStruct.h"

#include "THStack.h"
#include "TLegend.h"

//forward declare so we don't bleed NuOscillator headers
class OscillatorBase;

/// @brief Class responsible for handling implementation of samples used in analysis, reweighting and returning LLH
/// @author Dan Barrow
/// @author Ed Atkin
class samplePDFFDBase :  public samplePDFBase
{
public:
  //######################################### Functions #########################################
  /// @param ConfigFileName Name of config to initialise the sample object
  samplePDFFDBase(std::string ConfigFileName, covarianceXsec* xsec_cov, covarianceOsc* osc_cov = nullptr);
  /// @brief destructor
  virtual ~samplePDFFDBase();

  int GetNDim(){return nDimensions;} //DB Function to differentiate 1D or 2D binning
  std::string GetName() const {return samplename;}

  std::string GetXBinVarName() {return XVarStr;}
  std::string GetYBinVarName() {return YVarStr;}

  void PrintIntegral(TString OutputName="/dev/null", int WeightStyle=0, TString OutputCSVName="/dev/null");
  
  //===============================================================================
  // DB Reweighting and Likelihood functions

  //ETA - abstract these to samplePDFFDBase
  //DB Require these four functions to allow conversion from TH1(2)D to array for multi-threaded GetLikelihood
  void addData(TH1D* Data) override;
  void addData(TH2D* Data) override;
  void addData(std::vector<double> &data) override;
  void addData(std::vector< std::vector <double> > &data) override;
  /// @brief DB Multi-threaded GetLikelihood
  double GetLikelihood() override;
  //===============================================================================

  void reweight() override;
  M3::float_t GetEventWeight(const int iSample, const int iEntry) const;

  ///  @brief including Dan's magic NuOscillator
  void SetupNuOscillator();

  virtual void setupSplines(FarDetectorCoreInfo *, const char *, int , int ){};

  void ReadSampleConfig();

  int getNMCSamples() override {return int(MCSamples.size());}

  int getNEventsInSample(int iSample) {
    if (iSample < 0 || iSample > getNMCSamples()) {
      MACH3LOG_ERROR("Invalid Sample Requested: {}",iSample);
      throw MaCh3Exception(__FILE__ , __LINE__);
    }
    return MCSamples[iSample].nEvents;
  }
  
  std::string getFlavourName(int iSample) {
    if (iSample < 0 || iSample > getNMCSamples()) {
      MACH3LOG_ERROR("Invalid Sample Requested: {}",iSample);
      throw MaCh3Exception(__FILE__ , __LINE__);      
    }
    return MCSamples[iSample].flavourName;
  }

  TH1* get1DVarHist(std::string ProjectionVar, std::vector< std::vector<double> > SelectionVec = std::vector< std::vector<double> >(), int WeightStyle=0, TAxis* Axis=nullptr);
  TH2* get2DVarHist(std::string ProjectionVarX, std::string ProjectionVarY, std::vector< std::vector<double> > SelectionVec = std::vector< std::vector<double> >(), int WeightStyle=0, TAxis* AxisX=nullptr, TAxis* AxisY=nullptr);

  TH1* get1DVarHistByModeAndChannel(std::string ProjectionVar_Str, int kModeToFill=-1, int kChannelToFill=-1, int WeightStyle=0, TAxis* Axis=nullptr);
  TH2* get2DVarHistByModeAndChannel(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, int kModeToFill=-1, int kChannelToFill=-1, int WeightStyle=0, TAxis* AxisX=nullptr, TAxis* AxisY=nullptr);

  TH1 *getModeHist1D(int s, int m, int style = 0) {
    return get1DVarHistByModeAndChannel(XVarStr,m,s,style);
  }
  TH2 *getModeHist2D(int s, int m, int style = 0) {
    return get2DVarHistByModeAndChannel(XVarStr,YVarStr,m,s,style);
  }

  std::vector<TH1*> ReturnHistsBySelection1D(std::string KinematicProjection, int Selection1, int Selection2=-1, int WeightStyle=0, TAxis* Axis=0);
  std::vector<TH2*> ReturnHistsBySelection2D(std::string KinematicProjectionX, std::string KinematicProjectionY, int Selection1, int Selection2=-1, int WeightStyle=0, TAxis* XAxis=0, TAxis* YAxis=0);
  THStack* ReturnStackedHistBySelection1D(std::string KinematicProjection, int Selection1, int Selection2=-1, int WeightStyle=0, TAxis* Axis=0);
  TLegend* ReturnStackHistLegend() {return THStackLeg;}
  
  /// ETA function to generically convert a string from xsec cov to a kinematic type
  int ReturnKinematicParameterFromString(const std::string& KinematicStr) const;
  /// ETA function to generically convert a kinematic type from xsec cov to a string
  std::string ReturnStringFromKinematicParameter(const int KinematicVariable) const;

 protected:
  /// @brief DB Function to determine which weights apply to which types of samples pure virtual!!
  virtual void SetupWeightPointers() = 0;

  /// @brief Ensure Kinematic Map is setup and make sure it is initialised correctly
  void SetupKinematicMap();

  /// @todo abstract the spline initialisation completely to core
  /// @brief initialise your splineXX object and then use InitialiseSplineObject to conviently setup everything up
  virtual void SetupSplines() = 0;

  //DB Require all objects to have a function which reads in the MC
  /// @brief Initialise any variables that your experiment specific samplePDF needs
  virtual void Init() = 0;

  /// @brief Experiment specific setup, returns the number of events which were loaded
  virtual int setupExperimentMC(int iSample) = 0;

  /// @brief Function which translates experiment struct into core struct
  virtual void setupFDMC(int iSample) = 0;

  /// @brief Function which does a lot of the lifting regarding the workflow in creating different MC objects
  void Initialise();
  
  /// @brief Contains all your binned splines and handles the setup and the returning of weights from spline evaluations
  std::unique_ptr<splineFDBase> SplineHandler;
  //===============================================================================
  void fillSplineBins();

  //Functions which find the nominal bin and bin edges
  void FindNominalBinAndEdges1D();
  void FindNominalBinAndEdges2D();

  //DB Overrided functions of base class which calculate erec bin and boundaries for reweighting speedup in beam samples
  //ETA - this can be done using the core info stored in the new fdmc_struct
  void set1DBinning(int nbins, double* boundaries);
  void set1DBinning(int nbins, double low, double high);
  void set2DBinning(int nbins1, double* boundaries1, int nbins2, double* boundaries2);
  void set2DBinning(int nbins1, double low1, double high1, int nbins2, double low2, double high2);
  void set1DBinning(std::vector<double> &XVec);
  void set2DBinning(std::vector<double> &XVec, std::vector<double> &YVec);
  void SetupSampleBinning();
  std::string XVarStr, YVarStr;
  std::vector<std::string> SplineVarNames;
  std::vector<double> SampleXBins;
  std::vector<double> SampleYBins;
  //===============================================================================

  /// @brief ETA - a function to setup and pass values to functional parameters where you need to pass a value to some custom reweight calc or engine
  virtual void SetupFunctionalParameters(){};
  /// @brief Update the functional parameter values to the latest propsed values. Needs to be called before every new reweight so is called in fillArray 
  virtual void PrepFunctionalParameters(){};
  /// @brief ETA - generic function applying shifts
  virtual void applyShifts(int iSample, int iEvent){(void) iSample; (void) iEvent;};
  /// @brief DB Function which determines if an event is selected, where Selection double looks like {{ND280KinematicTypes Var1, douuble LowBound}
  bool IsEventSelected(const int iSample, const int iEvent);

  /// @brief Check whether a normalisation systematic affects an event or not
  void CalcXsecNormsBins(int iSample);
  /// @brief Calculate the spline weight for a given event
  M3::float_t CalcWeightSpline(const int iSample, const int iEvent) const;
  /// @brief Calculate the norm weight for a given event
  M3::float_t CalcWeightNorm(const int iSample, const int iEvent) const;

  /// @brief Calculate weights for function parameters
  ///
  /// First you need to setup additional pointers in you experiment code in SetupWeightPointers
  /// Then in this function you can calculate whatever fancy function you want by filling weight to which you have pointer
  /// This way func weight shall be used in GetEventWeight
  virtual void CalcWeightFunc(int iSample, int iEvent){return; (void)iSample; (void)iEvent;};

  virtual double ReturnKinematicParameter(std::string KinematicParamter, int iSample, int iEvent) = 0;
  virtual double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) = 0;

  virtual std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter) = 0; //Returns binning for parameter Var
  virtual const double* GetPointerToKinematicParameter(std::string KinematicParamter, int iSample, int iEvent) = 0; 
  virtual const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) = 0;

  void SetupNormParameters();

  //===============================================================================
  //DB Functions required for reweighting functions
  //DB Replace previous implementation with reading bin contents from samplePDF_array
  void fill1DHist();
  void fill2DHist();

  /// @brief DB Nice new multi-threaded function which calculates the event weights and fills the relevant bins of an array
#ifdef MULTITHREAD
  /// @brief fills the samplePDFFD_array vector with the weight calculated from reweighting but multithreaded
  void fillArray_MP();
#endif
  /// @brief fills the samplePDFFD_array vector with the weight calculated from reweighting
  void fillArray();

  /// @brief Helper function to reset histograms
  inline void ResetHistograms();
  
  //===============================================================================
  //DB Variables required for GetLikelihood
  //
  /// DB Vectors to hold bin edges
  std::vector<double> XBinEdges;
  std::vector<double> YBinEdges;

  /// DB Array to be filled after reweighting
  double** samplePDFFD_array;
  /// KS Array used for MC stat
  double** samplePDFFD_array_w2;
  /// DB Array to be filled in AddData
  double** samplePDFFD_data;
  //===============================================================================

  //===============================================================================
  //MC variables
  std::vector<FarDetectorCoreInfo> MCSamples;
  //===============================================================================

  //===============================================================================
  /// DB Variables required for oscillation
  std::vector<OscillatorBase*> NuOscProbCalcers;
  std::string NuOscillatorConfigFile; 
  //=============================================================================== 

  //===============================================================================
  //DB Covariance Objects
  //ETA - All experiments will need an xsec, det and osc cov
  //these should be added to samplePDFBase to be honest
  covarianceXsec *XsecCov = nullptr;
  covarianceOsc *OscCov = nullptr;

  /// @brief flag used to define whether all oscillation channels have a probability calculated using the same binning
  bool EqualBinningPerOscChannel = false;
  
  //=============================================================================== 

  /// @brief Keep track of the dimensions of the sample binning
  int nDimensions = M3::_BAD_INT_;
  /// @brief A unique ID for each sample based on powers of two for quick binary operator comparisons 
  std::string SampleDetID;
  /// holds "TrueNeutrinoEnergy" and the strings used for the sample binning.
  std::vector<std::string> SplineBinnedVars;

  /// @brief the name of this sample e.g."muon-like"
  std::string samplename;

  /// @brief Information to store for normalisation pars
  std::vector<XsecNorms4> xsec_norms;
  /// pointer to osc params, since not all params affect every sample, we perform some operations before hand for speed
  std::vector<const double*> OscParams;
  //===========================================================================
  //DB Vectors to store which kinematic cuts we apply
  //like in XsecNorms but for events in sample. Read in from sample yaml file 
  //What gets used in IsEventSelected, which gets set equal to user input plus 
  //all the vectors in StoreSelection
  /// @brief the Number of selections in the 
  int NSelections = M3::_BAD_INT_;
  
  /// @brief What gets pulled from config options, these are constant after loading in
  /// this is of length 3: 0th index is the value, 1st is lower bound, 2nd is upper bound
  std::vector< std::vector<double> > StoredSelection;
  /// @brief the strings grabbed from the sample config specifying the selections
  std::vector< std::string > SelectionStr; 
  /// @brief the bounds for each selection lower and upper
  std::vector< std::vector<double> > SelectionBounds;
  /// @brief a way to store selection cuts which you may push back in the get1DVar functions
  /// most of the time this is just the same as StoredSelection
  std::vector< std::vector<double> > Selection;
   //===========================================================================

  /// Mapping between string and kinematic enum
  const std::unordered_map<std::string, int>* KinematicParameters;
  /// Mapping between kinematic enum and string
  const std::unordered_map<int, std::string>* ReversedKinematicParameters;

  std::unique_ptr<manager> SampleManager;
  void InitialiseSingleFDMCObject(int iSample, int nEvents);
  void InitialiseSplineObject();

  std::vector<std::string> oscchan_flavnames;
  std::vector<std::string> oscchan_flavnames_Latex;
  std::vector<std::string> mc_files;
  std::vector<std::string> spline_files;
  std::vector<int> sample_nupdg;
  std::vector<int> sample_nupdgunosc;

  std::unordered_map<std::string, double> _modeNomWeightMap;
  
  //===============================================================================
  /// DB Miscellaneous Variables
  TLegend* THStackLeg = nullptr;
  //===============================================================================
};
