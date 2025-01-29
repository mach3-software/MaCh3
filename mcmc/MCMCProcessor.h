#pragma once

#ifndef _UNDEF_
#define _UNDEF_ 1234567890
#endif

// C++ includes
#include <complex>
#include <cstdio>

// MaCh3 includes
#include "mcmc/StatisticalUtils.h"
#include "samplePDF/HistogramUtils.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TColor.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TText.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TROOT.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TGraphPolar.h"
#include "TCandle.h"
#include "TMath.h"
#include "TMatrixDSymEigen.h"
#include "TVirtualFFT.h"
_MaCh3_Safe_Include_End_ //}


//KS: Joy of forward declaration https://gieseanw.wordpress.com/2018/02/25/the-joys-of-forward-declarations-results-from-the-real-world/
class TChain;
class TF1;

/// @todo KS: Apply reweighted weight to plotting and Bayes Factor.
/// @todo KS: Implement 2D reweighing like DayaBay.
/// @todo KS: Implement Diagnostics/GetPenaltyTerm.cpp here.

/// KS: Enum for different covariance classes
enum ParameterEnum {
  kXSecPar  = 0,
  kNDPar    = 1,
  kFDDetPar = 2,
  kOSCPar   = 3,
  
  kNParameterEnum = 4 //KS: keep it at the end to keep track of all parameters
};

/// @brief Class responsible for processing MCMC chains, performing diagnostics, generating plots, and managing Bayesian analysis.
/// @details This class provides utilities to handle MCMC output generated by mcmc::runMCMC. It is particularly useful for extracting values from previous MCMC runs and initiating new MCMC runs with those values. Inspired by nd280_utils/DrawComp.cpp.
/// @see For more details and examples, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/09.-Bayesian-Analysis,-Plotting-and-MCMC-Processor).
/// @author Clarence Wret
/// @author Kamil Skwarczynski
class MCMCProcessor {
  public:
    /// @brief Constructs an MCMCProcessor object with the specified input file and options.
    /// @param InputFile The path to the input file containing MCMC data.
    MCMCProcessor(const std::string &InputFile);
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
    /// @param Mute Allow silencing many messages, especially important if we calculate matrix many times
    void MakeCovariance_MP(const bool Mute = false);
    /// @brief Make and Draw SubOptimality
    /// @cite roberts2009adaptive
    /// @author Henry Wallace
    void MakeSubOptimality(const int NIntervals = 10);

    /// @brief Reset 2D posteriors, in case we would like to calculate in again with different BurnInCut
    void ResetHistograms();
        
    /// @brief Draw the post-fit comparisons
    void DrawPostfit();
    /// @brief Make and Draw Violin
    void MakeViolin();
    /// @brief Make and Draw Credible intervals
    /// @param CredibleIntervals Vector with values of credible intervals, must be in descending order
    /// @param CredibleIntervalsColours Color_t telling what colour to use for each Interval line
    /// @param CredibleInSigmas Bool telling whether intervals are in percentage or in sigmas, then special conversions is used
    void MakeCredibleIntervals(const std::vector<double>& CredibleIntervals = {0.99, 0.90, 0.68 },
                               const std::vector<Color_t>& CredibleIntervalsColours = {kCyan+4, kCyan-2, kCyan-10},
                               const bool CredibleInSigmas = false
                               );
    /// @brief Draw the post-fit covariances
    void DrawCovariance();
    /// @brief Make and Draw Credible Regions
    /// @param CredibleRegions Vector with values of credible intervals, must be in descending order
    /// @param CredibleRegionStyle Style_t telling what line style to use for each Interval line
    /// @param CredibleRegionColor Color_t telling what colour to use for each Interval line
    /// @param CredibleInSigmas Bool telling whether intervals are in percentage or in sigmas, then special conversions is used
    void MakeCredibleRegions(const std::vector<double>& CredibleRegions = {0.99, 0.90, 0.68},
                             const std::vector<Style_t>& CredibleRegionStyle = {kDashed, kSolid, kDotted},
                             const std::vector<Color_t>& CredibleRegionColor = {kGreen-3, kGreen-10, kGreen},
                             const bool CredibleInSigmas = false
                             );
    /// @brief Make fancy triangle plot for selected parameters
    /// @param CredibleIntervals Vector with values of credible intervals, must be in descending order
    /// @param CredibleIntervalsColours Color_t telling what colour to use for each Interval line
    /// @param CredibleInSigmas Bool telling whether intervals are in percentage or in sigmas, then special conversions is used
    /// @param CredibleRegions Vector with values of credible intervals, must be in descending order
    /// @param CredibleRegionStyle Style_t telling what line style to use for each Interval line
    /// @param CredibleRegionColor Color_t telling what colour to use for each Interval line
    /// @param CredibleInSigmas Bool telling whether intervals are in percentage or in sigmas, then special conversions is used
    void MakeTrianglePlot(const std::vector<std::string>& ParNames,
                          // 1D
                          const std::vector<double>& CredibleIntervals = {0.99, 0.90, 0.68 },
                          const std::vector<Color_t>& CredibleIntervalsColours = {kCyan+4, kCyan-2, kCyan-10},
                          //2D
                          const std::vector<double>& CredibleRegions = {0.99, 0.90, 0.68},
                          const std::vector<Style_t>& CredibleRegionStyle = {kDashed, kSolid, kDotted},
                          const std::vector<Color_t>& CredibleRegionColor = {kGreen-3, kGreen-10, kGreen},
                          // Other
                          const bool CredibleInSigmas = false
                          );

    /// @brief Checks the order and size consistency of the `CredibleIntervals` and `CredibleIntervalsColours` vectors.
    /// @param CredibleIntervals A vector of credible interval values.
    /// @param CredibleIntervalsColours A vector of colors associated with each credible interval.
    /// @throws MaCh3Exception If the sizes are not equal or the intervals are not in decreasing order.
    void CheckCredibleIntervalsOrder(const std::vector<double>& CredibleIntervals, const std::vector<Color_t>& CredibleIntervalsColours);

    /// @brief Checks the order and size consistency of the `CredibleRegions`, `CredibleRegionStyle`, and `CredibleRegionColor` vectors.
    /// @param CredibleRegions A vector of credible region values.
    /// @param CredibleRegionStyle A vector of styles associated with each credible region.
    /// @param CredibleRegionColor A vector of colors associated with each credible region.
    /// @throws MaCh3Exception If the sizes are not equal or the regions are not in decreasing order.
    void CheckCredibleRegionsOrder(const std::vector<double>& CredibleRegions,
                                   const std::vector<Style_t>& CredibleRegionStyle,
                                   const std::vector<Color_t>& CredibleRegionColor);

    /// @brief Make funny polar plot
    /// @param ParNames Vector with parameter names for which Polar Plot will be made
    void GetPolarPlot(const std::vector<std::string>& ParNames);

    /// @brief Calculate Bayes factor for vector of params, and model boundaries
    /// @param ParName Vector with parameter names for which we calculate Bayes factor
    /// @param Model1Bounds Lower and upper bound for hypothesis 1. Within this bound we calculate integral used later for Bayes Factor
    /// @param Model2Bounds Lower and upper bound for hypothesis 2. Within this bound we calculate integral used later for Bayes Factor
    /// @param ModelNames Names for hypothesis 1 and 2
    void GetBayesFactor(const std::vector<std::string>& ParName,
                        const std::vector<std::vector<double>>& Model1Bounds,
                        const std::vector<std::vector<double>>& Model2Bounds,
                        const std::vector<std::vector<std::string>>& ModelNames);
    /// @brief Calculate Bayes factor for point like hypothesis using SavageDickey
    void GetSavageDickey(const std::vector<std::string>& ParName,
                         const std::vector<double>& EvaluationPoint,
                         const std::vector<std::vector<double>>& Bounds);
    /// @brief Reweight Prior by giving new central value and new error
    /// @param ParName Parameter names for which we do reweighting
    /// @param NewCentral New central value for which we reweight
    /// @param NewError New error used for calculating weight
    void ReweightPrior(const std::vector<std::string>& Names,
                       const std::vector<double>& NewCentral,
                       const std::vector<double>& NewError);
    
    /// @brief Make .gif of parameter evolution
    /// @param ParName Parameter names for which we do .gif
    /// @param NIntervals Number of intervals for a gif
    void ParameterEvolution(const std::vector<std::string>& Names,
                            const std::vector<int>& NIntervals);

    /// @brief Thin MCMC Chain, to save space and maintain low autocorrelations.
    /// @param ThinningCut every which entry you want to thin
    /// @param Average If true will perform MCMC averaging instead of thinning
    inline void ThinMCMC(const int ThinningCut) { ThinningMCMC(MCMCFile+".root", ThinningCut); };

    /// @brief KS: Perform MCMC diagnostic including Autocorrelation, Trace etc.
    void DiagMCMC();
    
    // Get the number of parameters
    /// @brief Get total number of used parameters
    inline int GetNParams() { return nDraw; };
    inline int GetNXSec() { return nParam[kXSecPar]; };
    inline int GetNND() { return nParam[kNDPar]; };
    inline int GetNFD() { return nParam[kFDDetPar]; };
    inline int GetOSC() { return nParam[kOSCPar]; };
    /// @brief Number of params from a given group, for example flux
    int GetGroup(const std::string& name) const;

    /// @brief Get 1D posterior for a given parameter
    /// @param i parameter index
    inline TH1D* GetHpost(const int i) { return hpost[i]; };
    /// @brief Get 2D posterior for a given parameter combination
    /// @param i parameter index X
    /// @param j parameter index Y
    inline TH2D* GetHpost2D(const int i, const int j) { return hpost2D[i][j]; };
    /// @brief Get Violin plot for all parameters with posterior values
    inline TH2D* GetViolin() { return hviolin; };
    /// @brief Get Violin plot for all parameters with prior values
    inline TH2D* GetViolinPrior() { return hviolin_prior; };

    //Covariance getters
    inline std::vector<std::string> GetXSecCov()  const { return CovPos[kXSecPar]; };
    inline std::string GetNDCov() const { return CovPos[kNDPar].back(); };
    inline std::string GetFDCov() const { return CovPos[kFDDetPar].back(); };
    inline std::vector<std::string> GetOscCov()   const { return CovPos[kOSCPar]; };

    /// @brief Get the post-fit results (arithmetic and Gaussian)
    void GetPostfit(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Central_Gauss, TVectorD *&Errors_Gauss, TVectorD *&Peaks);
    /// @brief Get the post-fit covariances and correlations
    /// @param Cov Covariance matrix
    /// @param Corr Correlation matrix
    void GetCovariance(TMatrixDSym *&Cov, TMatrixDSym *&Corr);
    /// @brief Or the individual post-fits
    void GetPostfit_Ind(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Peaks, ParameterEnum kParam);
    
    /// @brief Get the vector of branch names from root file
    const std::vector<TString>& GetBranchNames() const { return BranchNames;};
    /// @brief Get properties of parameter by passing it number
    void GetNthParameter(const int param, double &Prior, double &PriorError, TString &Title);
    /// @brief Get parameter number based on name
    int GetParamIndexFromName(const std::string& Name);
    /// @brief Get Number of entries that Chain has, for merged chains will not be the same Nsteps
    inline Long64_t GetnEntries(){return nEntries;};
    /// @brief Get Number of Steps that Chain has, for merged chains will not be the same nEntries
    inline Long64_t GetnSteps(){return nSteps;};
    
    /// @brief Set the step cutting by string
    /// @param Cuts string telling cut value
    void SetStepCut(const std::string& Cuts);
    /// @brief Set the step cutting by int
    /// @param Cuts integer telling cut value
    void SetStepCut(const int Cuts);

    /// @brief You can set relative to prior or relative to generated. It is advised to use relate to prior
    /// @param PlotOrNot bool controlling plotRelativeToPrior argument
    inline void SetPlotRelativeToPrior(const bool PlotOrNot){plotRelativeToPrior = PlotOrNot; };
    inline void SetPrintToPDF(const bool PlotOrNot){printToPDF = PlotOrNot; };
    /// @brief Set whether you want to plot error for parameters which have flat prior
    inline void SetPlotErrorForFlatPrior(const bool PlotOrNot){PlotFlatPrior = PlotOrNot; };
    inline void SetPlotBinValue(const bool PlotOrNot){plotBinValue = PlotOrNot; };
    inline void SetFancyNames(const bool PlotOrNot){FancyPlotNames = PlotOrNot; };
    /// @brief Set whether want to use smoothing for histograms using ROOT algorithm
    inline void SetSmoothing(const bool PlotOrNot){ApplySmoothing = PlotOrNot; };
    /// @brief Code will only plot 2D posteriors if Correlation are larger than defined threshold
    /// @param Threshold This threshold is compared with correlation value
    inline void SetPost2DPlotThreshold(const double Threshold){Post2DPlotThreshold = Threshold; };
    /// @brief Toggle using the FFT-based autocorrelation calculator
    inline void SetUseFFTAutoCorrelation(const bool useFFT){useFFTAutoCorrelation = useFFT; };

    /// @brief Setter related what parameters we want to exclude from analysis, for example if cross-section parameters look like xsec_, then passing "xsec_" will
    /// @param Batches Vector with parameters type names we want to exclude
    inline void SetExcludedTypes(std::vector<std::string> Name){ExcludedTypes = Name; };
    inline void SetExcludedNames(std::vector<std::string> Name){ExcludedNames = Name; };

    /// @brief Set value of Nbatches used for batched mean, this need to be done earlier as batches are made when reading tree
    /// @param Batches Number of batches, default is 20
    inline void SetnBatches(const int Batches){nBatches = Batches; };
    inline void SetnLags(const int nLags){AutoCorrLag = nLags; };
    
    /// @brief Sett output suffix, this way jobs using the same file will have different names
    inline void SetOutputSuffix(const std::string Suffix){OutputSuffix = Suffix; };
    /// @brief Allow to set addtional cuts based on ROOT TBrowser cut, for to only affect one mass ordering
    inline void SetPosterior1DCut(const std::string Cut){Posterior1DCut = Cut; };
  private:
    /// @brief Prepare prefit histogram for parameter overlay plot
    inline std::unique_ptr<TH1D> MakePrefit();
    /// @brief prepare output root file and canvas to which we will save EVERYTHING
    inline void MakeOutputFile();
    /// @brief Draw 1D correlations which might be more helpful than looking at huge 2D Corr matrix
    inline void DrawCorrelations1D();

    /// @brief CW: Read the input Covariance matrix entries. Get stuff like parameter input errors, names, and so on
    inline void ReadInputCov();
    /// @brief Read the output MCMC file and find what inputs were used
    /// @warning There is bit of hardcoding for names so we should revisit it
    inline void FindInputFiles();
    /// @brief Read the xsec file and get the input central values and errors
    inline void ReadXSecFile();
    /// @brief Read the ND cov file and get the input central values and errors
    inline void ReadNDFile();
    /// @brief Read the FD cov file and get the input central values and errors
    inline void ReadFDFile();
    /// @brief Read the Osc cov file and get the input central values and errors
    inline void ReadOSCFile();
    /// @brief Remove parameter specified in config
    inline void RemoveParameters();
    /// @brief Print info like how many params have been loaded etc
    inline void PrintInfo() const;

    /// @brief Scan Input etc.
    inline void ScanInput();
    /// @brief Scan order of params from a different groups
    inline void ScanParameterOrder();
    /// @brief Prepare all objects used for output
    inline void SetupOutput();

    // MCMC Diagnostic
    /// @brief CW: Prepare branches etc. for DiagMCMC
    inline void PrepareDiagMCMC();
    /// @brief CW: Draw trace plots of the parameters i.e. parameter vs step
    inline void ParamTraces();
    /// @brief KS: Calculate autocorrelations supports both OpenMP and CUDA :)
    inline void AutoCorrelation();
    /// @brief MJR: Autocorrelation function using FFT algorithm for extra speed
    /// @author Michael Reh
    inline void AutoCorrelation_FFT();
    /// @brief KS: calc Effective Sample Size
    /// @param nLags Should be the same nLags as used in AutoCorrelation()
    /// @param LagL Value of LagL for each dial and each Lag
    ///
    /// This function computes the Effective Sample Size (ESS) using the autocorrelations
    /// calculated by AutoCorrelation(). Ensure that the parameter nLags here matches
    /// the number of lags used in AutoCorrelation() to obtain accurate results.
    /// @cite StanManual
    /// @cite hanson2008mcmc
    /// @cite gabry2024visual
    inline void CalculateESS(const int nLags, const std::vector<std::vector<double>>& LagL);
    /// @brief Get the batched means variance estimation and variable indicating if number of batches is sensible
    /// @cite chakraborty2019estimating
    /// @cite rossetti2024batch
    inline void BatchedAnalysis();
    /// @brief CW: Batched means, literally read from an array and chuck into TH1D
    inline void BatchedMeans();
    /// @brief Geweke Diagnostic based on the methods described by Fang (2014) and Karlsbakk (2011).
    /// @cite Fang2014GewekeDiagnostics
    /// @cite karlsbakk2011
    inline void GewekeDiagnostic();
    /// @brief Acceptance Probability
    inline void AcceptanceProbabilities();
    /// @brief RC: Perform spectral analysis of MCMC
    /// @cite Dunkley:2004sv
    /// @author Richard Calland
    inline void PowerSpectrumAnalysis();

    /// Name of MCMC file
    std::string MCMCFile;
    /// Output file suffix useful when running over same file with different settings
    std::string OutputSuffix;
    /// Covariance matrix name position
    std::vector<std::vector<std::string>> CovPos;
    /// Covariance matrix config
    std::vector<YAML::Node> CovConfig;

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
    /// Number of all parameters used in the analysis
    int nDraw;

    //Name of all branches as well as branches we don't want to include in the analysis
    std::vector<TString> BranchNames;
    std::vector<std::string> ExcludedTypes;
    std::vector<std::string> ExcludedNames;
    
    /// Is the ith parameter varied
    std::vector<bool> IamVaried;
    /// Name of parameters which we are going to analyse
    std::vector<std::vector<TString>> ParamNames;
    /// Parameters central values which we are going to analyse
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
    // KS: For example flux or detector within matrix
    std::vector<std::string> ParameterGroup;

    /// Vector of each systematic
    std::vector<TString> SampleName_v;
    /// Vector of each sample PDF object
    std::vector<TString> SystName_v;
    
    /// Name of output files
    std::string OutputName;
    /// Name of canvas which help to save to the sample pdf
    TString CanvasName;

    /// Whether we plot flat prior or not, we usually provide error even for flat prior params
    bool PlotFlatPrior;
    /// Will plot Jarlskog Invariant using information in the chain
    bool PlotJarlskog;
    
    //Even more flags
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
    /// MJR: Use FFT-based autocorrelation algorithm (save time & resources)?
    bool useFFTAutoCorrelation;

    std::vector<int> NDSamplesBins;
    std::vector<std::string> NDSamplesNames;

    /// Gaussian fitter
    TF1 *Gauss;

    /// The output file
    TFile *OutputFile;
    
    /// Fancy canvas used for our beautiful plots
    std::unique_ptr<TCanvas> Posterior;

    //Vector of best fit points and errors obtained with different methods
    /// Vector with central value for each parameter
    TVectorD *Central_Value;
    /// Vector with mean values using Arithmetic Mean
    TVectorD *Means;
    /// Vector with errors values using RMS
    TVectorD *Errors;
    /// Vector with mean values using Gaussian fit
    TVectorD *Means_Gauss;
    /// Vector with error values using Gaussian fit
    TVectorD *Errors_Gauss;
    /// Vector with mean values using Highest Posterior Density
    TVectorD *Means_HPD;
    /// Vector with error values using Highest Posterior Density
    TVectorD *Errors_HPD;
    /// Vector with positive error (right hand side) values using Highest Posterior Density
    TVectorD *Errors_HPD_Positive;
    /// Vector with negative error (left hand side) values using Highest Posterior Density
    TVectorD *Errors_HPD_Negative;

    /// Posterior Covariance Matrix
    TMatrixDSym *Covariance;
    /// Posterior Correlation Matrix
    TMatrixDSym *Correlation;

    /// Holds 1D Posterior Distributions
    std::vector<TH1D*> hpost;
    /// Holds 2D Posterior Distributions
    std::vector<std::vector<TH2D*>> hpost2D;
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
    
    /// Total parameter sum for each param
    double *ParamSums;
    /// Values of batched average for every param and batch
    double **BatchedAverages;

    /// Holds the sample values
    double **SampleValues;
    /// Holds the systs values
    double **SystValues;

    /// Holds all accProb
    double *AccProbValues;
    /// Holds all accProb in batches
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
