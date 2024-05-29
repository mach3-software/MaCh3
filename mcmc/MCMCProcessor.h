#pragma once

#ifndef _UNDEF_
#define _UNDEF_ 1234567890
#endif

// C++ includes
#include <complex>

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
#include "TMatrixDSymEigen.h"

// MaCh3 includes
#include "mcmc/StatisticalUtils.h"

// Class to process MCMC output produced by mcmc::runMCMC
// Useful for when we want to extract values from a previous MCMC 
//  and want to send off another MCMC with those values, or perhaps to LLH scans at central values
// Mostly taken from nd280_utils/DrawComp.cpp
//
// Return: Postfit parameters, Postfit covariance matrix
//
// TODO
// Apply reweighted weight to plotting and Bayes Factor
// Implement Diagnostics/GetPenaltyTerm.cpp here

/// KS: Enum for different covariance classes
enum ParameterEnum {
  kXSecPar  = 0, //KS: This hold both xsec and flux
  kNDPar = 1,
  kFDDetPar = 2,
  kOSCPar   = 3,
  
  kNParameterEnum = 4 //KS: keep it at the end to keep track of all parameters
};

/// @brief Class responsible for processing MCMC chains and performing diagnostic, making plots etc
class MCMCProcessor {
  public:
    /// @brief Constructs an MCMCProcessor object with the specified input file and options.
    /// @param InputFile The path to the input file containing MCMC data.
    /// @param MakePostfitCorr A boolean indicating whether to apply post-fit corrections during processing.
    MCMCProcessor(const std::string &InputFile, bool MakePostfitCorr);
    /// @brief Destroys the MCMCProcessor object.
    virtual ~MCMCProcessor();

    /// @brief Scan chain, what parameters we have and load information from covariance matrices
    void Initialise();
    /// @brief Make 1D projection for each parameter and prepare structure
    void MakePostfit();
    /// @brief Calculate covariance by making 2D projection of each combination of parameters
    void MakeCovariance();
    /// @brief KS:By caching each step we use multithreading
    void CacheSteps();
    /// @brief Calculate covariance by making 2D projection of each combination of parameters using multithreading
    void MakeCovariance_MP(bool Mute = false);
    /// @brief Make and Draw SubOptimality
    void MakeSubOptimality(int NIntervals = 10);

    /// @brief Reset 2D posteriors, in case we would like to calculate in again with different BurnInCut
    void ResetHistograms();
        
    /// @brief Draw the post-fit comparisons
    void DrawPostfit();
    /// @brief Make and Draw Violin
    void MakeViolin();
    /// @brief Make and Draw Credible intervals
    void MakeCredibleIntervals(std::vector<double> CredibleIntervals = {0.99, 0.90, 0.68 },
                               std::vector<Color_t> CredibleIntervalsColours = {kCyan+4, kCyan-2, kCyan-10},
                               bool CredibleInSigmas = false
                               );
    /// @brief Draw the post-fit covariances
    void DrawCovariance();
    /// @brief Make and Draw Credible Regions
    void MakeCredibleRegions(std::vector<double> CredibleRegions = {0.99, 0.90, 0.68},
                             std::vector<Style_t> CredibleRegionStyle = {kDashed, kSolid, kDotted},
                             std::vector<Color_t> CredibleRegionColor = {kGreen-3, kGreen-10, kGreen},
                             bool CredibleInSigmas = false
                             );
    /// @brief Make fancy triangle plot for selected parameters
    void MakeTrianglePlot(std::vector<std::string> ParNames,
                          // 1D
                          std::vector<double> CredibleIntervals = {0.99, 0.90, 0.68 },
                          std::vector<Color_t> CredibleIntervalsColours = {kCyan+4, kCyan-2, kCyan-10},
                          //2D
                          std::vector<double> CredibleRegions = {0.99, 0.90, 0.68},
                          std::vector<Style_t> CredibleRegionStyle = {kDashed, kSolid, kDotted},
                          std::vector<Color_t> CredibleRegionColor = {kGreen-3, kGreen-10, kGreen},
                          // Other
                          bool CredibleInSigmas = false
                          );
    /// @brief Make funny polar plot
    void GetPolarPlot(std::vector<std::string> ParNames);

    /// @brief Calculate Bayes factor for vector of params, and model boundaries
    void GetBayesFactor(std::vector<std::string> ParName, std::vector<std::vector<double>> Model1Bounds, std::vector<std::vector<double>> Model2Bounds, std::vector<std::vector<std::string>> ModelNames);
    /// @brief Calculate Bayes factor for point like hypothesis using SavageDickey
    void GetSavageDickey(std::vector<std::string> ParName, std::vector<double> EvaluationPoint, std::vector<std::vector<double>> Bounds);
    /// @brief Reweight Prior by giving new central value and new error
    void ReweightPrior(std::vector<std::string> Names, std::vector<double> NewCentral, std::vector<double> NewError);
    
    /// @brief KS: Perform MCMC diagnostic including AutoCorrelation, Trace etc.
    void DiagMCMC();
    
    // Get the number of parameters
    inline int GetNParams() { return nDraw; };
    inline int GetNFlux() { return nFlux; };
    inline int GetNXSec() { return nParam[kXSecPar]; };
    inline int GetNND() { return nParam[kNDPar]; };
    inline int GetNFD() { return nParam[kFDDetPar]; };
    inline int GetOSC() { return nParam[kOSCPar]; };
        
    /// @brief Get 1D posterior for a given parameter
    inline TH1D* GetHpost(const int i) { return hpost[i]; };
    /// @brief Get 2D posterior for a given parameter combination
    inline TH2D* GetHpost2D(const int i, const int j) { return hpost2D[i][j]; };
    /// @brief Get Violin plot for all parameters with posterior values
    inline TH2D* GetViolin() { return hviolin; };
    /// @brief Get Violin plot for all parameters with prior values
    inline TH2D* GetViolinPrior() { return hviolin_prior; };

    //Covariance getters
    inline std::vector<std::string> const & GetXSecCov()  const { return CovPos[kXSecPar]; };
    inline std::string const & GetNDCov() const { return CovPos[kNDPar].back(); };
    inline std::string const & GetFDCov()    const { return CovPos[kFDDetPar].back(); };
    inline std::string const & GetOscCov()   const { return CovPos[kOSCPar].back(); };

    /// @brief Get the post-fit results (arithmetic and Gaussian)
    void GetPostfit(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Central_Gauss, TVectorD *&Errors_Gauss, TVectorD *&Peaks);
    /// @brief Get the post-fit covariances and correlations
    void GetCovariance(TMatrixDSym *&Cov, TMatrixDSym *&Corr);
    /// @brief Or the individual post-fits
    void GetPostfit_Ind(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Peaks, ParameterEnum kParam);
    
    /// @brief Get the vector of branch names from root file
    const std::vector<TString>& GetBranchNames() const { return BranchNames;};
    /// @brief Get properties of parameter by passing it number
    void GetNthParameter(const int param, double &Prior, double &PriorError, TString &Title);
    /// @brief Get parameter number based on name
    int GetParamIndexFromName(const std::string Name);
    /// @brief Get Number of entries that Chain has, for merged chains will not be the same Nsteps
    inline int GetnEntries(){return nEntries;};
    /// @brief Get Number of Steps that Chain has, for merged chains will not be the same nEntries
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
    /// @biref Code will only plot 2D posteriors if Correlation are larger than defined threshold
    inline void SetPost2DPlotThreshold(const double Threshold){Post2DPlotThreshold = Threshold; };

    //Setter related what parameters we want to exclude from analysis
    inline void SetExcludedTypes(std::vector<std::string> Name){ExcludedTypes = Name; };
    inline void SetExcludedNames(std::vector<std::string> Name){ExcludedNames = Name; };

    //DiagMCMC-related setter
    inline void SetnBatches(const int Batches){nBatches = Batches; };
    inline void SetnLags(const int nLags){AutoCorrLag = nLags; };
    
    inline void SetOutputSuffix(const std::string Suffix){OutputSuffix = Suffix; };
    inline void SetPosterior1DCut(const std::string Cut){Posterior1DCut = Cut; };
  private:
    /// @brief Prepare prefit histogram for parameter overlay plot
    inline TH1D* MakePrefit();
    /// @brief prepare output root file and canvas to which we will save EVERYTHING
    inline void MakeOutputFile();
    /// @brief Draw 1D correlations which might be more helpful than looking at huge 2D Corr matrix
    inline void DrawCorrelations1D();

    /// @brief Read Matrices
    inline void ReadInputCov();
    inline void FindInputFiles();
    inline void ReadXSecFile();
    inline void ReadNDFile();
    inline void ReadFDFile();
    inline void ReadOSCFile();
    inline void RemoveParameters();
   
    /// @brief Scan Input etc.
    inline void ScanInput();
    /// @brief Scan order of params from a different groups
    inline void ScanParameterOrder();
    /// @brief Prepare all objects used for output
    inline void SetupOutput();

    //Analyse posterior distribution
    /// @brief Get Arithmetic mean from posterior
    inline void GetArithmetic(TH1D * const hist, const int i);
    /// @brief Fit Gaussian to posterior
    inline void GetGaussian(TH1D *& hist, const int i);
    /// @brief Get Highest Posterior Density (HPD)
    inline void GetHPD(TH1D * const hist, const int i, const double coverage = 0.6827);
    /// @brief Get 1D Credible Interval
    inline void GetCredibleInterval(TH1D* const hist, TH1D* hpost_copy, const double coverage = 0.6827);
    /// @brief Get 2D Credible Region
    inline void GetCredibleRegion(TH2D* hist2D, const double coverage = 0.6827);

    // MCMC Diagnostic
    /// @brief CW: Prepare branches etc. for DiagMCMC
    inline void PrepareDiagMCMC();
    /// @brief CW: Draw trace plots of the parameters i.e. parameter vs step
    inline void ParamTraces();
    /// @brief KS: Calculate autocorrelations supports both OpenMP and CUDA :)
    inline void AutoCorrelation();
    /// @brief KS: calc Effective Sample Size Following https://mc-stan.org/docs/2_18/reference-manual/effective-sample-size-section.html
    inline void CalculateESS(const int nLags);
    /// @brief Get the batched means variance estimation and variable indicating if number of batches is sensible
    inline void BatchedAnalysis();
    /// @brief CW: Batched means, literally read from an array and chuck into TH1D
    inline void BatchedMeans();
    /// @brief Geweke Diagnostic based on https://www.math.arizona.edu/~piegorsch/675/GewekeDiagnostics.pdf
    inline void GewekeDiagnostic();
    /// @brief Acceptance Probability
    inline void AcceptanceProbabilities();
    /// @brief RC: Perform spectral analysis of MCMC based on http://arxiv.org/abs/astro-ph/0405462
    inline void PowerSpectrumAnalysis();

    /// Name of MCMC file
    std::string MCMCFile;
    /// Output file suffix useful when running over same file with different settings
    std::string OutputSuffix;
    /// Covariance matrix name position
    std::vector<std::vector<std::string>> CovPos;

    /// Main chain storing all steps etc
    TChain *Chain;
    /// BurnIn Cuts
    std::string StepCut;
    /// Cut used when making 1D Posterior distribution
    std::string Posterior1DCut;
    /// KS: Used only for SubOptimality
    int UpperCut;
    /// Value of burn in cut
    int BurnInCut;
    /// Number of branches in a TTree
    int nBranches;
    /// KS: For merged chains number of entries will be different from nSteps
    int nEntries;
    /// KS: For merged chains number of entries will be different from nSteps
    int nSteps;
    /// Number of sample PDF objects
    int nSamples;
    /// Number of covariance objects
    int nSysts;

    //Name of all branches as well as branches we don't want to include in the analysis
    std::vector<TString> BranchNames;
    std::vector<std::string> ExcludedTypes;
    std::vector<std::string> ExcludedNames;
    
    /// Number of all parameters used in the analysis
    int nDraw;
    
    /// Is the ith parameter varied
    std::vector<bool> IamVaried;
    std::vector<std::vector<TString>> ParamNames;
    std::vector<std::vector<double>>  ParamCentral;
    std::vector<std::vector<double>>  ParamNom;
    std::vector<std::vector<double>>  ParamErrors;
    std::vector<std::vector<bool>>    ParamFlat;
    /// Number of parameters per type
    std::vector<int> nParam;
    /// Make an enum for which class this parameter belongs to so we don't have to keep string comparing
    std::vector<ParameterEnum> ParamType;
    /// KS: in MCMC output there is order of parameters so for example first goes xsec then nd det etc.
    /// Idea is that this parameter will keep track of it so code is flexible
    std::vector<int> ParamTypeStartPos;
    
    //In XsecMatrix we have both xsec and flux parameters, this is just for some plotting options
    std::vector<bool>   IsXsec; 
    /// This keep number of Flux params in xsec matrix
    int nFlux;

    // Vector of each systematic
    std::vector<TString> SampleName_v;
    std::vector<TString> SystName_v;
    
    /// Name of output files
    std::string OutputName;
    /// Name of canvas which help to save to the sample pdf
    TString CanvasName;

    //Plotting flags
    bool PlotXSec;
    bool PlotDet;
    /// whether we plot flat prior or not
    bool PlotFlatPrior;
    /// Will plot Jarlskog Invariant using information in the chain
    bool PlotJarlskog;
    
    //Even more flags
    /// Make correlation matrix or not
    bool MakeCorr;
    /// Whether we plot relative to prior or nominal, in most cases is prior
    bool plotRelativeToPrior;
    /// Sanity check if Postfit is already done to not make several times
    bool MadePostfit;
    /// Will plot all plot to PDF not only to root file
    bool printToPDF;
    /// Whether we want fancy plot names or not
    bool FancyPlotNames;
    /// If true it will print value on each bin of covariance matrix
    bool plotBinValue;
    /// Apply smoothing for 2D histos using root algorithm
    bool ApplySmoothing;
    /// KS: Set Threshold when to plot 2D posterior as by default we get a LOT of plots
    double Post2DPlotThreshold;

    std::vector< int > NDSamplesBins;
    std::vector< std::string > NDSamplesNames;

    /// Gaussian fitter
    TF1 *Gauss;

    /// The output file
    TFile *OutputFile;
    
    /// Fancy canvas used for our beautiful plots
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

    /// Posterior Covariance Matrix
    TMatrixDSym *Covariance;
    /// Posterior Correlation Matrix
    TMatrixDSym *Correlation;

    /// Holds 1D Posterior Distributions
    TH1D **hpost;
    /// Holds 2D Posterior Distributions
    TH2D ***hpost2D;
    /// Holds violin plot for all dials
    TH2D *hviolin;
    /// Holds prior violin plot for all dials,
    TH2D *hviolin_prior;

    /// Array holding values for all parameters
    double** ParStep;
    /// Step number for step, important if chains were merged
    int* StepNumber;

    /// Number of bins
    int nBins;
    /// Drawrange for SetMaximum
    double DrawRange;
    
    /// MCMC Chain has been cached
    bool CacheMCMC;
    /// Doing MCMC Diagnostic
    bool doDiagMCMC;
    
    //Number of batches and LagL used in MCMC diagnostic
    /// Number of batches for Batched Mean
    int nBatches;
    /// LagL used in AutoCorrelation
    int AutoCorrLag;
    
    // Holds all the parameter variations
    double *ParamSums;
    double **BatchedAverages;
    double **LagL;

    /// Holds the sample values
    double **SampleValues;
    /// Holds the systs values
    double **SystValues;

    // Holds all accProb
    double *AccProbValues;
    double *AccProbBatchedAverages;
    
  //Only if GPU is enabled
  #ifdef CUDA
    /// @brief Move stuff to GPU to perform auto correlation calculations there
    inline void PrepareGPU_AutoCorr(const int nLags);

    /// Value of each param that will be copied to GPU
    float* ParStep_cpu;
    float* NumeratorSum_cpu;
    float* ParamSums_cpu;
    float* DenomSum_cpu;

    /// Value of each param at GPU
    float* ParStep_gpu;
    float* NumeratorSum_gpu;
    float* ParamSums_gpu;
    float* DenomSum_gpu;
  #endif
};
