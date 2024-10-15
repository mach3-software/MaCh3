#pragma once

//C++ includes
#include <list>

//ROOT includes
#include "THStack.h"
#include "TLegend.h"

//MaCh3 includes
#include "OscProbCalcer/OscProbCalcerBase.h"
#include "Oscillator/OscillatorBase.h"

#include "splines/splineFDBase.h"

#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"

#include "samplePDF/samplePDFBase.h"
#include "samplePDF/FDMCStruct.h"
#include "samplePDF/ShiftFunctors.h"

/// @brief Class responsible for handling implementation of samples used in analysis, reweighting and returning LLH
class samplePDFFDBase :  public samplePDFBase
{
public:
  //######################################### Functions #########################################

  samplePDFFDBase(){};
  samplePDFFDBase(std::string mc_version, covarianceXsec* xsec_cov);
  virtual ~samplePDFFDBase();

  int GetNDim(){return nDimensions;} //DB Function to differentiate 1D or 2D binning
  std::string GetName(){return samplename;}

  //===============================================================================
  // DB Reweighting and Likelihood functions

  //ETA - abstract these to samplePDFFDBase
  //DB Require these four functions to allow conversion from TH1(2)D to array for multi-threaded GetLikelihood
  void addData(TH1D* Data);
  void addData(TH2D* Data);
  void addData(std::vector<double> &data);
  void addData(std::vector< std::vector <double> > &data);
  /// @brief DB Multi-threaded GetLikelihood
  double GetLikelihood();
  //===============================================================================

  void reweight();
  double GetEventWeight(int iSample, int iEntry);

  // Setup and config functions
  void UseNonDoubledAngles(bool ans) {doubled_angle = ans;};
  
  const double **oscpars;
  void SetXsecCov(covarianceXsec* xsec_cov);
  void SetOscCov(covarianceOsc* osc_cov);

  ///  @brief including Dan's magic NuOscillator
  void SetupNuOscillator();

  /// @deprecated The `DumpWeights` function is deprecated and should not be used.
  /// It was kept for backwards compatibility in compiling but has no effect.
  ///
  /// @note This function was marked for deprecation as of 14/01/2015 by KD.
  ///       - DB (27/08/2020): The function is incredibly hardcoded.
  ///       - DB Consider using 'LetsPrintSomeWeights' to achieve the same functionality.
  ///
  /// @param outname The name of the output file.
  virtual void DumpWeights(std::string outname) {(void)outname; return; };
  //================================================================================

  virtual void setupSplines(fdmc_base *skobj, const char *SplineFileName, int nutype, int signal){};
  void ReadSampleConfig();

 protected:
  /// @brief DB Function to determine which weights apply to which types of samples pure virtual!!
  virtual void SetupWeightPointers() = 0;

  /// @todo abstract the spline initialisation completely to core
  /// @brief initialise your splineXX object and then use InitialiseSplineObject to conviently setup everything up
  virtual void SetupSplines() = 0;

  //DB Require all objects to have a function which reads in the MC
  // @brief Initialise any variables that your experiment specific samplePDF needs
  virtual void Init() = 0;

  //DB Experiment specific setup, returns the number of events which were loaded
  virtual int setupExperimentMC(int iSample) = 0;

  //DB Function which translates experiment struct into core struct
  virtual void setupFDMC(int iSample) = 0;

  //DB Function which does a lot of the lifting regarding the workflow in creating different MC objects
  void Initialise();
  
  splineFDBase *splineFile;
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
  virtual void PrepFunctionalParameters(){};
  /// @brief ETA - generic function applying shifts
  virtual void applyShifts(int iSample, int iEvent){(void) iSample; (void) iEvent;};
  /// @brief DB Function which determines if an event is selected, where Selection double looks like {{ND280KinematicTypes Var1, douuble LowBound}
  bool IsEventSelected(const int iSample, const int iEvent);
  bool IsEventSelected(const std::vector<std::string>& ParameterStr, const int iSample, const int iEvent);
  bool IsEventSelected(const std::vector<std::string>& ParameterStr, const std::vector<std::vector<double>> &SelectionCuts, const int iSample, const int iEvent);

  /// @brief Check whether a normalisation systematic affects an event or not
  void CalcXsecNormsBins(int iSample);
  /// @brief Is the sample for when operating in Reverse Horn Current, read in from sample config
  bool GetIsRHC() {return IsRHC;}
  /// @brief Calculate the spline weight for a given event
  double CalcXsecWeightSpline(const int iSample, const int iEvent);
  /// @brief Calculate the norm weight for a given event
  double CalcXsecWeightNorm(const int iSample, const int iEvent);
  /// @brief Virtual so this can be over-riden in an experiment derived class
  virtual double CalcXsecWeightFunc(int iSample, int iEvent){(void)iSample; (void)iEvent; return 1.0;};

  virtual double ReturnKinematicParameter(std::string KinematicParamter, int iSample, int iEvent) = 0;
  virtual double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) = 0;
  virtual std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter) = 0; //Returns binning for parameter Var
  virtual const double* ReturnKinematicParameterByReference(std::string KinematicParamter, int iSample, int iEvent) = 0; 
  virtual const double* ReturnKinematicParameterByReference(double KinematicVariable, int iSample, int iEvent) = 0;

  //ETA - new function to generically convert a string from xsec cov to a kinematic type
  virtual inline int ReturnKinematicParameterFromString(std::string KinematicStr) = 0;
  virtual inline std::string ReturnStringFromKinematicParameter(int KinematicVariable) = 0;

  // Function to setup Functional and shift parameters. This isn't idea but
  // do this in your experiment specific code for now as we don't have a 
  // generic treatment for func and shift pars yet
  virtual void SetupFuncParameters(){return;};
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
  std::vector<struct fdmc_base> MCSamples;
  TFile *_sampleFile;
  TTree *_data;
  //===============================================================================

  //===============================================================================
  /// DB Variables required for oscillation
  std::vector<OscillatorBase*> NuOscProbCalcers;
  std::string NuOscillatorConfigFile;
  
  //===============================================================================
  
  //Variables controlling oscillation parameters
  bool doubled_angle;

  //===============================================================================
  //DB Covariance Objects
  //ETA - All experiments will need an xsec, det and osc cov
  //these should be added to samplePDFBase to be honest
  covarianceXsec *XsecCov;
  covarianceOsc *OscCov;

  //=============================================================================== 

  /// @brief Keep track of the dimensions of the sample binning
  int nDimensions;
  /// @brief A unique ID for each sample based on powers of two for quick binary operator comparisons 
  int SampleDetID;
  /// @breif Is the sample for events collected in Reverse Horn Current. Important for flux systematics
  bool IsRHC;
  /// holds "TrueNeutrinoEnergy" and the strings used for the sample binning.
  std::vector<std::string> SplineBinnedVars;

  std::string samplename;

  /// Information to store for normalisation pars
  std::vector<XsecNorms4> xsec_norms;
  int nFuncParams;
  std::vector<std::string> funcParsNames;
  std::vector<int> funcParsIndex;

  //===========================================================================
  //DB Vectors to store which kinematic cuts we apply
  //like in XsecNorms but for events in sample. Read in from sample yaml file 
  //What gets used in IsEventSelected, which gets set equal to user input plus 
  //all the vectors in StoreSelection
  /// @brief the Number of selections in the 
  int NSelections;
  
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

  manager* SampleManager;
  void InitialiseSingleFDMCObject(int iSample, int nEvents);
  void InitialiseSplineObject();

  double Unity = 1.;
  double Zero = 0.;
  
  std::vector<std::string> mtuple_files;
  std::vector<std::string> spline_files;
  std::vector<int> sample_vecno;
  std::vector<int> sample_oscnutype;
  std::vector<int> sample_nutype;
  std::vector<bool> sample_signal;

  std::string mtupleprefix;
  std::string mtuplesuffix;
  std::string splineprefix;
  std::string splinesuffix;
};
