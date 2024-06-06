#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "samplePDF/samplePDFBase.h"
#include "covariance/covarianceBase.h"
#include "covariance/covarianceOsc.h"

#include "manager/manager.h"
#include "mcmc/MCMCProcessor.h"

//KS: Joy of forward declaration https://gieseanw.wordpress.com/2018/02/25/the-joys-of-forward-declarations-results-from-the-real-world/
class TRandom3;
class TStopwatch;
class TTree;
class TGraphAsymmErrors;
class TDirectory;

/// @brief Base class for implementing fitting algorithms
class FitterBase {

 public:
  /// @brief Constructor
  FitterBase(manager * const fitMan);
  /// @brief Destructor
  virtual ~FitterBase();

  /// @brief This function adds a sample PDF object to the analysis framework. The sample PDF object will be utilized in fitting procedures or likelihood scans.
  /// @param sample A pointer to a sample PDF object derived from samplePDFBase.
  void addSamplePDF(samplePDFBase* sample);

  /// @brief This function adds a Covariance object to the analysis framework. The Covariance object will be utilized in fitting procedures or likelihood scans.
  /// @param cov A pointer to a Covariance object derived from covarianceBase.
  void addSystObj(covarianceBase* cov);

  /// @brief Adds an oscillation handler for covariance objects.
  /// @param oscf A pointer to a covarianceOsc object for forward oscillations.
  void addOscHandler(covarianceOsc* oscf);

  /// @brief Adds two oscillation handlers for covariance objects.
  /// @param osca A pointer to a covarianceOsc object for the first oscillation.
  /// @param oscb A pointer to a covarianceOsc object for the second oscillation.
  void addOscHandler(covarianceOsc* osca, covarianceOsc* oscb);

  /// @brief The specific fitting algorithm implemented in this function depends on the derived class. It could be Markov Chain Monte Carlo (MCMC), MinuitFit, or another algorithm.
  virtual void runMCMC() = 0;

  /// @brief Perform a 1D likelihood scan.
  void RunLLHScan();

  /// @brief LLH scan is good first estimate of step scale
  void GetStepScaleBasedOnLLHScan();

  /// @brief Perform a 2D likelihood scan.
  /// @warning This operation may take a significant amount of time, especially for complex models.
  void Run2DLLHScan();

  /// @brief Perform a 2D and 1D sigma var for all samples.
  /// @warning Code uses TH2Poly
  void RunSigmaVar();

  /// @brief Get name of class
  virtual inline std::string GetName()const {return "FitterBase";};
 protected:
  /// @brief Process MCMC output
  void ProcessMCMC();

  /// @brief Prepare the output file.
  void PrepareOutput();

  /// @brief Save output and close files.
  void SaveOutput();

  /// @brief Save the settings that the MCMC was run with.
  void SaveSettings();

  /// @brief Used by sigma variation, check how 1 sigma changes spectra
  inline TGraphAsymmErrors* MakeAsymGraph(TH1D* sigmaArrayLeft, TH1D* sigmaArrayCentr, TH1D* sigmaArrayRight, std::string title);

  /// The manager
  manager *fitMan;

  /// MaCh3 Modes
  MaCh3Modes* Modes;

  /// current state
  unsigned int step;
  /// current likelihood
  double logLCurr;
  /// proposed likelihood
  double logLProp;
  /// current acceptance prob
  double accProb;
  /// counts accepted steps
  int accCount;

  /// LLH for samples/syst objects
  /// oscillation covariance llh
  double osc_llh;
  /// store the llh breakdowns
  double *sample_llh;
  /// systematic llh breakdowns
  double *syst_llh;

  /// Sample holder
  std::vector<samplePDFBase*> samples;
  /// Total number of samples used
  unsigned int TotalNSamples;

  /// Systematic holder
  std::vector<covarianceBase*> systematics;

  /// handles oscillation parameters
  covarianceOsc *osc;
  /// handles oscillation parameters
  covarianceOsc *osc2;

  /// tells global time how long fit took
  TStopwatch* clock;
  /// tells how long single step/fit iteration took
  TStopwatch* stepClock;
  /// Time of single step
  double stepTime;

  /// Random number
  TRandom3* random;

  /// Output
  TFile *outputFile;
  /// Output cov folder
  TDirectory *CovFolder;
  /// Output tree with posteriors
  TTree *outTree;
  /// auto save every N steps
  int auto_save;

  /// Necessary for some fitting algorithms like PSO
  bool fTestLikelihood;
  /// save nominal matrix info or not
  bool save_nominal;
  /// Save proposal at each step
  bool SaveProposal;

  /// Checks if file saved not repeat some operations
  bool FileSaved;
  /// Checks if setting saved not repeat some operations
  bool SettingsSaved;
  /// Checks if output prepared not repeat some operations
  bool OutputPrepared;

  #ifdef DEBUG
  /// Debugging flag
  bool debug;
  /// Debugging Output file
  std::ofstream debugFile;
  #endif
};

