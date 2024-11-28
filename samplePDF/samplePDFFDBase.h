#ifndef _samplePDFFDBase_h_
#define _samplePDFFDBase_h_
//C++ includes
#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <list>
#include <vector>

//ROOT includes
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraph2DErrors.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

//Propagator includes
#ifdef USE_PROB3
  #include "BargerPropagator.h"
#else
  #include "beamcudapropagator.cuh"
  #include "atmoscudapropagator.cuh"
#endif

//MaCh3 includes
#include "interfacePDFEbE.h"
#include "samplePDFBase.h"

#include "splines/splineBase.h"
#include "splines/splineFDBase.h"

#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"

#include "FDMCStruct.h"
#include "ShiftFunctors.h"

#include "manager/manager.h"



#define USEBETA 0

class samplePDFFDBase : virtual public samplePDFBase , virtual public interfacePDFEbE
{
public:
  //######################################### Functions #########################################

  samplePDFFDBase(){};
  samplePDFFDBase(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFFDBase();

  int getNDim(); //DB Function to differentiate 1D or 2D binning

  //===============================================================================
  // DB Reweighting and Likelihood functions

  //ETA - abstract these to samplePDFFDBase
  //DB Require these four functions to allow conversion from TH1(2)D to array for multi-threaded GetLikelihood
  void addData(TH1D* Data);
  void addData(TH2D* Data);
  void addData(std::vector<double> &data);
  void addData(std::vector< vector <double> > &data);
  //DB Multi-threaded GetLikelihood
  double GetLikelihood();
  //===============================================================================

  void reweight(double *oscpar);
  double getEventWeight(int iSample, int iEntry);

  // Setup and config functions
  void useNonDoubledAngles(bool ans) {doubled_angle = ans;};
  void useBinnedOscReweighting(bool ans);
  void useBinnedOscReweighting(bool ans, int nbins, double *osc_bins);
  
  void setUseBeta(bool ub) { // Apply beta to nuebar
    useBeta = ub; 
    if (useBeta){std::cout << "useBeta implementation hidden behind '#ifdef USEBETA' in samplePDFSKBase - Check USEBETA is true, recompile and remove this comment" << std::endl; throw;}
  };
  void setApplyBetaNue(bool nueb) { // Apply beta to nue instead of nuebar
    applyBetaNue = nueb;
    if (applyBetaNue){std::cout << "useBeta implementation hidden behind '#ifdef USEBETA' in samplePDFSKBase - Check USEBETA is true, recompile and remove this comment" << std::endl; throw;}
  };
  void setApplyBetaDiag(bool diagb){ // Apply 1/beta to nue and beta to nuebar
    applyBetaDiag = diagb;
    if (applyBetaDiag){std::cout << "useBeta implementation hidden behind '#ifdef USEBETA' in samplePDFSKBase - Check USEBETA is true, recompile and remove this comment" << std::endl; throw;}
  }; 

#if defined (USE_PROB3) && defined (CPU_ONLY)
  inline double calcOscWeights(int sample, int nutype, int oscnutype, double en, double *oscpar);
#endif

#if defined (USE_PROB3) && not defined (CPU_ONLY)
  void calcOscWeights(int nutype, int oscnutype, double *en, double *w, int num, double *oscpar);
#endif

#if not defined (USE_PROB3)
  void calcOscWeights(int sample, int nutype, double *w, double *oscpar);
#endif

  std::string GetSampleName(){return samplename;}

  int GetBinningOpt(){return BinningOpt;}

  //============================= Should be deprecated =============================
  // Note: the following functions aren't used any more! (From 14/1/2015) - KD. Just kept in for backwards compatibility in compiling, but they have no effect.
  // DB 27/08/2020 The following functions shouldn't be used (Currently included for backwards compatibility)

  //DB Incredibly hardcoded - Could probably make 'LetsPrintSomeWeights' do the same functionality and remove this?
  //void DumpWeights(std::string outname);

  // ETA - in the future it would be nice to have some generic getHIst functions
  // although, this introduces a root dependence into the core code?
  //TH1D *getModeHist1D(int s, int m, int style = 0);
  //TH2D *getModeHist2D(int s, int m, int style = 0);
  // Direct translation of getModeHist1D(s,m,style) = get1DVarHist(kPDFBinning,m,s,style)
  //TH1D* get1DVarHist(ND280KinematicTypes Var1, int fModeToFill=-1, int fSampleToFill=-1, int WeightStyle=0, TAxis* Axis=0);
  //TH1D* get1DVarHist(ND280KinematicTypes Var1, std::vector< std::vector<double> > Selection, int WeightStyle=0, TAxis* Axis=0);
  // Direct translation of getModeHist2D(s,m,style) = get2DVarHist(kPDFBinning,kPDFBinning,m,s,style)
  //TH2D* get2DVarHist(ND280KinematicTypes Var1, ND280KinematicTypes Var2, int fModeToFill=-1, int fSampleToFill=-1, int WeightStyle=0, TAxis* XAxis=0, TAxis* YAxis=0);
  //TH2D* get2DVarHist(ND280KinematicTypes Var1, ND280KinematicTypes Var2, std::vector< std::vector<double> > Selection, int WeightStyle=0, TAxis* XAxis=0, TAxis* YAxis=0);
  //TH3D* get3DVarHist(ND280KinematicTypes Var1, ND280KinematicTypes Var2, ND280KinematicTypes Var3, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* XAxis=0, TAxis* YAxis=0, TAxis* ZAxis=0);
 
  //================================================================================

  virtual void setupSplines(fdmc_base *skobj, const char *splineFile, int nutype, int signal){};


  // LW - Setup Osc 
  void virtual SetupOscCalc(double PathLength, double Density);

  //get1DVarHist and get2DVarHist functions
  virtual TH1D* get1DVarHist(std::string KinematicVar1, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis) = 0; 
  virtual TH1D* get1DVarHist(std::string KinematicVar1,std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis) = 0; 
  virtual TH2D* get2DVarHist(std::string KinematicVar1,std::string KinematicVar2, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis, TAxis* Axis2) = 0; 
  virtual TH2D* get2DVarHist(std::string KinematicVar1,std::string KinematicVar2,std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis, TAxis* Axis2) = 0; 

 protected:


  //TODO - I think this will be tricky to abstract. fdmc_base will have to contain the pointers to the appropriate weights, can probably pass the number of these weights to constructor?
  //DB Function to determine which weights apply to which types of samples
  //pure virtual!!
  virtual void setupWeightPointers() = 0;


  //===============================================================================
  //DB Functions relating to sample and exec setup  
  //ETA - abstracting these core functions
  //init will setup all the specific variables 
  void init(double pot, std::string mc_version, covarianceXsec *xsec_cov){return;};
  //void setupMC(manager* sample_manager, const char *sampleInputFile, const char *splineFile, fdmc_base *fdobj, double pot, int nutype, int oscnutype, bool signal, int iSample, bool hasfloats=false){std::cout << "SAMPLEPDFFDBase::setupMC " << std::endl; return;};
  //virtual void setupSplines(fdmc_base *skobj, const char *splineFile, int nutype, int signal);

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
  //===============================================================================

  //ETA - generic function applying shifts
  virtual void applyShifts(int iSample, int iEvent){};
  //DB Function which determines if an event is selected, where Selection double looks like {{ND280KinematicTypes Var1, douuble LowBound}
  bool IsEventSelected(std::vector< std::string > ParameterStr, int iSample, int iEvent);
  bool IsEventSelected(std::vector< std::string > ParameterStr, std::vector< std::vector<double> > &Selection, int iSample, int iEvent);
  virtual void reconfigureFuncPars(){};

  void CalcXsecNormsBins(int iSample);
  //This just gets read in from a yaml file
  bool GetIsRHC() {return IsRHC;}
  // Calculate the spline weight for a given event
  double CalcXsecWeightSpline(const int iSample, const int iEvent);
  // Calculate the norm weight for a given event
  double CalcXsecWeightNorm(const int iSample, const int iEvent);
  virtual double CalcXsecWeightFunc(int iSample, int iEvent) = 0;


  //virtual double ReturnKinematicParameter(KinematicTypes Var, int i) = 0;       //Returns parameter Var for event j in sample i
  virtual double ReturnKinematicParameter(std::string KinematicParamter, int iSample, int iEvent) = 0;
  virtual std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter) = 0; //Returns binning for parameter Var
  //ETA - new function to generically convert a string from xsec cov to a kinematic type
  //virtual double StringToKinematicVar(std::string kinematic_str) = 0;
  virtual void setXsecCov(covarianceXsec* xsec_cov);

  //===============================================================================
  //DB Functions required for reweighting functions  
  //DB Replace previous implementation with reading bin contents from samplePDF_array
  void fill1DHist();
  void fill2DHist();

  
  //DB Nice new multi-threaded function which calculates the event weights and fills the relevant bins of an array
#ifdef MULTITHREAD
  void fillArray_MP();
#endif
  void fillArray();

  // Helper function to reset histograms
  inline void ResetHistograms();
 
#ifndef USE_PROB3 
  inline cudaprob3::ProbType SwitchToCUDAProbType(CUDAProb_nu CUDAProb_nu);  
#endif
  //===============================================================================
  //DB Variables required for GetLikelihood
  //
  //DB Vectors to hold bin edges
  std::vector<double> XBinEdges;
  std::vector<double> YBinEdges;

  //DB Array to be filled after reweighting
  double** samplePDFFD_array;
  //KS Array used for MC stat
  double** samplePDFFD_array_w2;
  //DB Array to be filled in AddData
  double** samplePDFFD_data;
  //===============================================================================

  //===============================================================================
  //MC variables
  vector<struct fdmc_base> MCSamples;
  TFile *_sampleFile;
  TTree *_data;
  //===============================================================================

  //===============================================================================
  //DB Variables required for oscillation

  // Propagator
  // Will initialise either Prob3++ or CUDAProb3 depending on USE_PROB3 flag
  // CUDA/CPU versions of CUDAProb3 set with CPU_ONLY flag
  // If USEPROB3=true and CPU_ONLY=false, probGpu is setup in source code

  // An axis to set binned oscillation weights
  TAxis *osc_binned_axis ;
  //===============================================================================
  
  //Variables controlling nuebar appearance
  double Beta;
  bool useBeta;
  bool applyBetaNue;
  bool applyBetaDiag;

  //Variables controlling oscillation parameters
  bool doubled_angle;
  bool osc_binned;

  //===============================================================================
  //DB Covariance Objects
  //covarianceSkDet_joint *skdet_joint;
  //ETA -  need to think about this? All experiments will need an xsec, det and osc cov
  covarianceXsec *xsecCov;
  //===============================================================================
  //

  //ETA - binning opt can probably go soon...
  int BinningOpt;
  int SampleDetID;
  bool IsRHC;

  std::string samplename;

  //Information to store for normalisation pars
  std::vector<XsecNorms4> xsec_norms;
  int nFuncParams;
  std::vector<std::string> funcParsNames;
  std::vector<int> funcParsIndex;

  //===============================================================================
  //DB Vectors to store which kinematic cuts we apply
  std::vector< std::string > SelectionStr; //like in XsecNorms but for events in sample. Read in from sample yaml file 
  std::vector< std::vector<double> > SelectionBounds; // like in XsecNorms but for events in sample. Read in from sample yaml file in samplePDFExperimentBase.cpp

  std::vector< std::vector<double> > Selection; //What gets used in IsEventSelected, which gets set equal to user input plus all the vectors in StoreSelection
  int NSelections;
  std::vector< std::vector<double> > StoredSelection; //What gets pulled from config options
  //===============================================================================
  //

  //ETA - trying new way of doing shift parameters
  // leave this for now but coming soon!
  //std::vector< EnergyScale* > ShiftFunctors;
  //std::vector< BaseFuncPar* > ShiftFunctors;

 private:
  int blah;

};
#endif
