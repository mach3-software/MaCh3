#ifndef __MCMCPROCESSOR_H_
#define __MCMCPROCESSOR_H_

#ifndef __UNDEF__
#define __UNDEF__ 1234567890
#endif
// C++ includes
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

// ROOT includes
#include "TObjArray.h"
#include "TObjString.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TColor.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TText.h"
#include "TGaxis.h"

// Class to process MCMC output produced by mcmc::runMCMC
// Useful for when we want to extract values from a previous MCMC 
//  and want to send off another MCMC with those values, or perhaps to LLH scans at central values
// Mostly taken from nd280_utils/DrawComp.cpp
//
// Return: Postfit parameters, Postfit covariance matrix
//
// Make postfit pmu cosmu dists
// Make LLH scans around output

enum ParameterEnum {
  kXSecPar  = 0, //KS: This hold both xsec and flux
  kND280Par = 1,
  kFDDetPar = 2,
  kOSCPar   = 3,
  
  kNParameterEnum = 4 //KS: keep it at the end to keep track of all parameters
};

class MCMCProcessor {
  public:
    MCMCProcessor(const std::string &InputFile, bool MakePostfitCorr);
    ~MCMCProcessor();

    void Initialise();
    // Get the post-fit results (arithmetic and Gaussian)
    void GetPostfit(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Central_Gauss, TVectorD *&Errors_Gauss, TVectorD *&Peaks);
    // Get the post-fit covariances and correlations
    void GetCovariance(TMatrixDSym *&Cov, TMatrixDSym *&Corr);
    // Or the individual post-fits
    void GetPostfit_Ind(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Peaks, ParameterEnum kParam);

    void MakeCovariance();
    void MakeCovariance_MP();
    
    void DiagMCMC();
    
    void MakePostfit();
    
    //KS:By caching each step we use multithreading
    void CacheSteps();
    void ResetHistograms();
    
    // Get the number of parameters
    int GetNParams() { return nDraw; };
    int GetNFlux() { return nFlux; };
    int GetNXSec() { return nParam[kXSecPar]; };
    int GetNND() { return nParam[kND280Par]; };
    int GetNFD() { return nParam[kFDDetPar]; };
    int GetOSC() { return nParam[kOSCPar]; };
        
    inline TH1D* const GetHpost(int i) { return hpost[i]; };

    std::string const & GetXSecCov()  const { return CovPos[kXSecPar]; };
    std::string const & GetND280Cov() const { return CovPos[kND280Par]; };
    std::string const & GetFDCov()    const { return CovPos[kFDDetPar]; };
    std::string const & GetOscCov()   const { return CovPos[kOSCPar]; };
    std::string const & GetNDruns()   const { return NDruns; };
    std::vector<std::string> const & GetNDsel() const {return NDsel;};

    // Draw the post-fit comparisons
    void DrawPostfit();
    // Make and Draw Violin
    void MakeViolin();
    // Make and Draw Credible intervals
    void MakeCredibleIntervals();
    // Draw the post-fit covariances
    void DrawCovariance();
    // Make and Draw Credible Regions
    void MakeCredibleRegions();
    //Make fancy triangle plot for selected parameters
    void MakeTrianglePlot(std::vector<std::string> ParamNames);

    // Get the vector of branch names
    const std::vector<TString>& GetBranchNames() const { return BranchNames;};

    void GetNthParameter(const int param, double &Prior, double &PriorError, TString &Title);
    void GetBayesFactor(std::string ParName, double M1_min, double M1_max, std::string M1Name, double M2_min, double M2_max, std::string M2Name);

    int GetnEntries(){return nEntries;};
    int GetnSteps(){return nSteps;};
    //KS: Many setters which in future will be loaded via config
    // Set the step cutting
    // Either by string
    void SetStepCut(std::string Cuts);
    // Or by int
    void SetStepCut(const int Cuts);

    void SetMakeOnlyXsecCorr(bool PlotOrNot){MakeOnlyXsecCorr = PlotOrNot; };
    void SetMakeOnlyXsecCorrFlux(bool PlotOrNot){MakeOnlyXsecCorrFlux = PlotOrNot; };

    void SetPlotRelativeToPrior(bool PlotOrNot){plotRelativeToPrior = PlotOrNot; };
    void SetPrintToPDF(bool PlotOrNot){printToPDF = PlotOrNot; };
    void SetPlotDet(bool PlotOrNot){PlotDet = PlotOrNot; };
    void SetPlotErrorForFlatPrior(bool PlotOrNot){PlotFlatPrior = PlotOrNot; };
    void SetPlotBinValue(bool PlotOrNot){plotBinValue = PlotOrNot; };
    void SetFancyNames(bool PlotOrNot){FancyPlotNames = PlotOrNot; };
    void SetSmoothing(bool PlotOrNot){ApplySmoothing = PlotOrNot; };
    void SetPost2DPlotThreshold(double Threshold){Post2DPlotThreshold = Threshold; };

    void SetnBatches(int Batches){nBatches = Batches; };
    void SetOutputSuffix(std::string Suffix){OutputSuffix = Suffix; };
    
  private:
    inline TH1D* MakePrefit();
    inline void MakeOutputFile();

    // Read Matrices
    inline void ReadInputCov();
    inline void FindInputFiles();
    inline void ReadXSecFile();
    inline void ReadND280File();
    inline void ReadFDFile();
    inline void ReadOSCFile();
   
    // Scan Input etc.
    inline void ScanInput();
    inline void ScanParameterOrder();
    inline void SetupOutput();

    //Analyse posterior distribution
    inline void GetArithmetic(TH1D * const hpost, const int i);
    inline void GetGaussian(TH1D *& hpost, const int i);
    inline void GetHPD(TH1D * const hpost, const int i, const double coverage = 0.6827);
    inline void GetCredibleInterval(TH1D* const hpost, TH1D* hpost_copy, const double coverage = 0.6827);
    inline void GetCredibleRegion(TH2D* hpost, const double coverage = 0.6827);
    inline std::string GetJeffreysScale(const double BayesFactor);


    // MCMC Diagnsotic
    inline void PrepareDiagMCMC();
    inline void ParamTraces();
    inline void AutoCorrelation();
    inline void CalculateESS(const int nLags);
    inline void BatchedAnalysis();
    inline void BatchedMeans();
    inline void GewekeDiagnostic();
    inline void AcceptanceProbabilities();
    
    std::string MCMCFile;
    std::string OutputSuffix;
    // Covariance matrix name position
    std::vector<std::string> CovPos;
    // ND runs
    std::string NDruns;
    // ND selections
    std::vector<std::string> NDsel;

    TChain *Chain;
    std::string StepCut;
    int BurnInCut;
    int nBranches;
    //KS: For merged chains number of entires will be different fron nSteps
    int nEntries;
    int nSteps;
    int nSamples;
    int nSysts;

    std::vector<TString> BranchNames;
    // Is the ith parameter varied
    std::vector<bool> IamVaried;
    std::vector<std::vector<TString>> ParamNames;
    std::vector<std::vector<double>>  ParamCentral;
    std::vector<std::vector<double>>  ParamNom;
    std::vector<std::vector<double>>  ParamErrors;
    std::vector<std::vector<bool>>    ParamFlat;
    // Make an enum for which class this parameter belongs to so we don't have to keep string comparing
    std::vector<ParameterEnum> ParamType;
    //KS: in MCMC output there is order of parameters so for example first goes xsec then nd det etc.
    //Idea is that this paraemter will keep track of it so code is flexible
    std::vector<int> ParamTypeStartPos;
    
    //In XsecMatrix we have both xsec and flux parameters, this is just for some ploting options
    std::vector<bool>   IsXsec; 
    
    // Vector of each systematic
    std::vector<TString> SampleName_v;
    std::vector<TString> SystName_v;
    
    std::string OutputName;
    TString CanvasName;

    bool PlotXSec;
    bool PlotDet;
    bool PlotFlatPrior;

    bool MakeCorr;
    bool MakeOnlyXsecCorr;
    bool MakeOnlyXsecCorrFlux; //Makes only the cov matrix w/ flux
    bool plotRelativeToPrior;
    bool MadePostfit;
    bool printToPDF;
    bool FancyPlotNames;
    bool plotBinValue; //If true it will print value on each bin of covariance matrix
    //KS: Set Threshold when to plot 2D posterior as by default we get a LOT of plots
    double Post2DPlotThreshold;

    // The output file
    TFile *OutputFile;

    int nDraw;
    std::vector<int> nParam;
    int nFlux; // This keep number of Flux params in xsec matrix

    std::vector< int > NDSamplesBins;
    std::vector< std::string > NDSamplesNames;

    // Gaussian fitter
    TF1 *Gauss;

    TCanvas *Posterior;

    TVectorD *Central_Value;
    TVectorD *Means;
    TVectorD *Errors;
    TVectorD *Means_Gauss;
    TVectorD *Errors_Gauss;
    TVectorD *Means_HPD;
    TVectorD *Errors_HPD; 
    TVectorD *Errors_HPD_Positive; 
    TVectorD *Errors_HPD_Negative; 

    TMatrixDSym *Covariance;
    TMatrixDSym *Correlation;

    bool CacheMCMC;
    bool doDiagMCMC;
    bool ApplySmoothing;
    // Holds Posterior Distributions
    TH1D **hpost;
    TH2D ***hpost2D;
    TH2D *hviolin;
    
    double** ParStep;
    int* StepNumber;

    // Number of bins
    int nBins;
    // Drawrange for SetMaximum
    double DrawRange;
    
    int nBatches;
    
    // Holds all the parameter variations
    double *ParamSums;
    double **BatchedAverages;
    double **LagL;

    // Holds the sample values
    double **SampleValues;
    // Holds the systs values
    double **SystValues;

    // Holds all accProb
    double *AccProbValues;
    double *AccProbBatchedAverages;
    
  //Only if GPU is enabled
  #ifdef CUDA
    void PrepareGPU_AutoCorr(const int nLags);

    float* ParStep_cpu;
    float* NumeratorSum_cpu;
    float* ParamSums_cpu;
    float* DenomSum_cpu;

    float* ParStep_gpu;
    float* NumeratorSum_gpu;
    float* ParamSums_gpu;
    float* DenomSum_gpu;
  #endif
};

#endif
