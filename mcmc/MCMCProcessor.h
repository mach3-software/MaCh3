#pragma once

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
#include "TObjString.h"
#include "TTree.h"
#include "TROOT.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TGraphPolar.h"
#include "TMath.h"

// Class to process MCMC output produced by mcmc::runMCMC
// Useful for when we want to extract values from a previous MCMC 
//  and want to send off another MCMC with those values, or perhaps to LLH scans at central values
// Mostly taken from nd280_utils/DrawComp.cpp
//
// Return: Postfit parameters, Postfit covariance matrix
//
// Make LLH scans around output

//KS: Enum for different covariance classes
enum ParameterEnum {
  kXSecPar  = 0, //KS: This hold both xsec and flux
  kNDPar = 1,
  kFDDetPar = 2,
  kOSCPar   = 3,
  
  kNParameterEnum = 4 //KS: keep it at the end to keep track of all parameters
};

class MCMCProcessor {
  public:
    MCMCProcessor(const std::string &InputFile, bool MakePostfitCorr);
    ~MCMCProcessor();

    // Scan chain, what paramters we have and load information from covariance matrices
    void Initialise();
    //Make 1D projection for each parameter and prepare structure
    void MakePostfit();
    //Calcuate covariance by making 2D projection of each combination of parameters
    void MakeCovariance();
    //KS:By caching each step we use multithreading
    void CacheSteps();
    void MakeCovariance_MP();
    //Reset 2D posteriors, in case we would like to calcualte in again with different BurnInCut
    void ResetHistograms();
        
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
    
    void GetPolarPlot(std::vector<std::string> ParNames);

    //Bayesian statistic hypotheis testing
    void GetBayesFactor(std::vector<std::string> ParName, std::vector<std::vector<double>> Model1Bounds, std::vector<std::vector<double>> Model2Bounds, std::vector<std::vector<std::string>> ModelNames);
    void GetSavageDickey(std::vector<std::string> ParName, std::vector<double> EvaluationPoint, std::vector<std::vector<double>> Bounds);
    void ReweightPrior(std::vector<std::string> Names, std::vector<double> NewCentral, std::vector<double> NewError);
    
    //KS: Perform MCMC diagnostic including AutoCorrelation, Trace etc.
    void DiagMCMC();
    
    // Get the number of parameters
    inline int GetNParams() { return nDraw; };
    inline int GetNFlux() { return nFlux; };
    inline int GetNXSec() { return nParam[kXSecPar]; };
    inline int GetNND() { return nParam[kNDPar]; };
    inline int GetNFD() { return nParam[kFDDetPar]; };
    inline int GetOSC() { return nParam[kOSCPar]; };
        
    //Posterior getters
    inline TH1D* GetHpost(const int i) { return hpost[i]; };
    inline TH2D* GetHpost2D(const int i, const int j) { return hpost2D[i][j]; };
    inline TH2D* GetViolin() { return hviolin; };
    inline TH2D* GetViolinPrior() { return hviolin_prior; };

    //Covariance getters
    inline std::string const & GetXSecCov()  const { return CovPos[kXSecPar]; };
    inline std::string const & GetNDCov() const { return CovPos[kNDPar]; };
    inline std::string const & GetFDCov()    const { return CovPos[kFDDetPar]; };
    inline std::string const & GetOscCov()   const { return CovPos[kOSCPar]; };
    //inline std::string const & GetNDruns()   const { return NDruns; };
    //inline std::vector<std::string> const & GetNDsel() const {return NDsel;};

    // Get the post-fit results (arithmetic and Gaussian)
    void GetPostfit(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Central_Gauss, TVectorD *&Errors_Gauss, TVectorD *&Peaks);
    // Get the post-fit covariances and correlations
    void GetCovariance(TMatrixDSym *&Cov, TMatrixDSym *&Corr);
    // Or the individual post-fits
    void GetPostfit_Ind(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Peaks, ParameterEnum kParam);
    
    // Get the vector of branch names
    const std::vector<TString>& GetBranchNames() const { return BranchNames;};
    //Getters
    void GetNthParameter(const int param, double &Prior, double &PriorError, TString &Title);
    int GetParamIndexFromName(const std::string Name);
    inline int GetnEntries(){return nEntries;};
    inline int GetnSteps(){return nSteps;};
    
    //KS: Many setters which in future will be loaded via config
    // Set the step cutting
    // Either by string
    void SetStepCut(std::string Cuts);
    // Or by int
    void SetStepCut(const int Cuts);

    //Setter related to plotting
    inline void SetPlotRelativeToPrior(const bool PlotOrNot){plotRelativeToPrior = PlotOrNot; };
    inline void SetPrintToPDF(const bool PlotOrNot){printToPDF = PlotOrNot; };
    inline void SetPlotErrorForFlatPrior(const bool PlotOrNot){PlotFlatPrior = PlotOrNot; };
    inline void SetPlotBinValue(const bool PlotOrNot){plotBinValue = PlotOrNot; };
    inline void SetFancyNames(const bool PlotOrNot){FancyPlotNames = PlotOrNot; };
    inline void SetSmoothing(const bool PlotOrNot){ApplySmoothing = PlotOrNot; };
    //Code will only plot 2D posteriors if Correlation are larger than defined threshold
    inline void SetPost2DPlotThreshold(const double Threshold){Post2DPlotThreshold = Threshold; };

    //Setter related what parameters we want to exclude from analysis
    inline void SetExcludedTypes(std::vector<std::string> Name){ExcludedTypes = Name; };
    inline void SetExcludedNames(std::vector<std::string> Name){ExcludedNames = Name; };

    //DiagMCMC-related setter
    inline void SetnBatches(const int Batches){nBatches = Batches; };
    inline void SetnLags(const int nLags){AutoCorrLag = nLags; };
    
    inline void SetOutputSuffix(const std::string Suffix){OutputSuffix = Suffix; };
    inline void SetPosterior1DCut(const std::string Cut){Posterior1DCut = Cut; };

    //Setter related to Credible Intervals or Regions
    inline void SetCredibleIntervals(std::vector<double> Intervals)
    {
      if(Intervals.size() > 1)
      {
        for(unsigned int i = 1; i < Intervals.size(); i++ )
        {
          if(Intervals[i] > Intervals[i-1])
          {
            std::cerr<<" Interval "<<i<<" is smaller than "<<i-1<<std::endl;
            std::cerr<<Intervals[i] <<" "<<Intervals[i-1]<<std::endl;
            std::cerr<<" They should be grouped in decreasing order"<<std::endl;
            std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
            throw;
          }
        }
      }
      Credible_Intervals = Intervals;
    };
    inline void SetCredibleIntervalsColours(std::vector<Color_t> Intervals){Credible_IntervalsColours = Intervals; };

    inline void SetCredibleRegions(std::vector<double> Intervals)
    {
      if(Intervals.size() > 1)
      {
        for(unsigned int i = 1; i < Intervals.size(); i++ )
        {
          if(Intervals[i] > Intervals[i-1])
          {
            std::cerr<<" Interval "<<i<<" is smaller than "<<i-1<<std::endl;
            std::cerr<<Intervals[i] <<" "<<Intervals[i-1]<<std::endl;
            std::cerr<<" They should be grouped in decreasing order"<<std::endl;
            std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
            throw;
          }
        }
      }
      Credible_Regions = Intervals;
    };
    inline void SetCredibleRegionStyle(std::vector<Style_t> Intervals){Credible_RegionStyle = Intervals; };
    inline void SetCredibleRegionColor(std::vector<Color_t> Intervals){Credible_RegionColor = Intervals; };
    inline void SetCredibleInSigmas(const bool Intervals){CredibleInSigmas = Intervals; };

  private:
    inline TH1D* MakePrefit();
    inline void MakeOutputFile();
    inline void DrawCorrelations1D();

    // Read Matrices
    inline void ReadInputCov();
    inline void FindInputFiles();
    inline void ReadXSecFile();
    inline void ReadNDFile();
    inline void ReadFDFile();
    inline void ReadOSCFile();
    inline void RemoveParameters();
   
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
    inline std::string GetDunneKaboth(const double BayesFactor);
    inline double GetSigmaValue(int sigma);

    // MCMC Diagnsotic
    inline void PrepareDiagMCMC();
    inline void ParamTraces();
    inline void AutoCorrelation();
    inline void CalculateESS(const int nLags);
    inline void BatchedAnalysis();
    inline void BatchedMeans();
    inline void GewekeDiagnostic();
    inline void AcceptanceProbabilities();
    
    //Usefull strings teling us about output etc
    std::string MCMCFile;
    std::string OutputSuffix;
    // Covariance matrix name position
    std::vector<std::string> CovPos;
    // ND runs
    std::string NDruns;
    // ND selections
    std::vector<std::string> NDsel;

    //Main chain storing all steps etc
    TChain *Chain;
    //BurnIn Cuts
    std::string StepCut;
    std::string Posterior1DCut;
    int BurnInCut;
    int nBranches;
    //KS: For merged chains number of entires will be different fron nSteps
    int nEntries;
    int nSteps;
    int nSamples;
    int nSysts;

    //Name of all branches as well as branches we don't want to incldue in the analysis
    std::vector<TString> BranchNames;
    std::vector<std::string> ExcludedTypes;
    std::vector<std::string> ExcludedNames;
    
    //Number of all parameters used in the analysis
    int nDraw;
    
    // Is the ith parameter varied
    std::vector<bool> IamVaried;
    std::vector<std::vector<TString>> ParamNames;
    std::vector<std::vector<double>>  ParamCentral;
    std::vector<std::vector<double>>  ParamNom;
    std::vector<std::vector<double>>  ParamErrors;
    std::vector<std::vector<bool>>    ParamFlat;
    //Number of parameters per type
    std::vector<int> nParam;
    // Make an enum for which class this parameter belongs to so we don't have to keep string comparing
    std::vector<ParameterEnum> ParamType;
    //KS: in MCMC output there is order of parameters so for example first goes xsec then nd det etc.
    //Idea is that this paraemter will keep track of it so code is flexible
    std::vector<int> ParamTypeStartPos;
    
    //In XsecMatrix we have both xsec and flux parameters, this is just for some ploting options
    std::vector<bool>   IsXsec; 
    int nFlux; // This keep number of Flux params in xsec matrix

    // Vector of each systematic
    std::vector<TString> SampleName_v;
    std::vector<TString> SystName_v;
    
    std::string OutputName;
    TString CanvasName;

    //Plotting flags
    bool PlotXSec;
    bool PlotDet;
    bool PlotFlatPrior;
    bool PlotJarlskog;
    
    //Even more falgse
    bool MakeCorr;
    bool plotRelativeToPrior;
    bool MadePostfit;
    bool printToPDF;
    bool FancyPlotNames;
    bool plotBinValue; //If true it will print value on each bin of covariance matrix
    bool ApplySmoothing;
    //KS: If true credible stuff is done in sigmas not %
    bool CredibleInSigmas;
    //KS: Set Threshold when to plot 2D posterior as by default we get a LOT of plots
    double Post2DPlotThreshold;

    std::vector< int > NDSamplesBins;
    std::vector< std::string > NDSamplesNames;

    // Gaussian fitter
    TF1 *Gauss;

    // The output file
    TFile *OutputFile;
    
    //Fancy canvas used for our beatufull plots
    TCanvas *Posterior;

    //Vector of best fit points and errors obtained with different methods
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

    // Holds Posterior Distributions
    TH1D **hpost;
    TH2D ***hpost2D;
    TH2D *hviolin;
    TH2D *hviolin_prior;

    double** ParStep;
    int* StepNumber;

    // Number of bins
    int nBins;
    // Drawrange for SetMaximum
    double DrawRange;
    
    //Flags related with MCMC diagnostic
    bool CacheMCMC;
    bool doDiagMCMC;
    
    //Number of batches and LagL used in MCMC diagnostic
    int nBatches;
    int AutoCorrLag;
    
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
    inline void PrepareGPU_AutoCorr(const int nLags);

    float* ParStep_cpu;
    float* NumeratorSum_cpu;
    float* ParamSums_cpu;
    float* DenomSum_cpu;

    float* ParStep_gpu;
    float* NumeratorSum_gpu;
    float* ParamSums_gpu;
    float* DenomSum_gpu;
  #endif

  //KS: Credible Intervals/Regions related
  std::vector<double> Credible_Intervals;
  std::vector<Color_t> Credible_IntervalsColours;

  std::vector<double> Credible_Regions;
  std::vector<Style_t> Credible_RegionStyle;
  std::vector<Color_t> Credible_RegionColor;
};
