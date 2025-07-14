#pragma once

// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>

// MaCh3 Includes
#include "Samples/SampleHandlerBase.h"
#include "Parameters/ParameterHandlerBase.h"
#include "Manager/Manager.h"
#include "Fitters/MCMCProcessor.h"

//KS: Joy of forward declaration https://gieseanw.wordpress.com/2018/02/25/the-joys-of-forward-declarations-results-from-the-real-world/
class TRandom3;
class TStopwatch;
class TTree;
class TGraphAsymmErrors;
class TDirectory;

/// @brief Base class for implementing fitting algorithms
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/06.-Fitting-Algorithms).
class FitterBase {
 public:
  /// @brief Constructor
  /// @param fitMan A pointer to a manager object, which will handle all settings.
  FitterBase(manager * const fitMan);
  /// @brief Destructor for the FitterBase class.
  virtual ~FitterBase();

  /// @brief This function adds a sample PDF object to the analysis framework. The sample PDF object will be utilized in fitting procedures or likelihood scans.
  /// @param sample A pointer to a sample PDF object derived from ParameterHandlerBase.
  void AddSampleHandler(SampleHandlerBase* sample);

  /// @brief This function adds a Covariance object to the analysis framework. The Covariance object will be utilized in fitting procedures or likelihood scans.
  /// @param cov A pointer to a Covariance object derived from ParameterHandlerBase.
  void AddSystObj(ParameterHandlerBase* cov);

  /// @brief The specific fitting algorithm implemented in this function depends on the derived class. It could be Markov Chain Monte Carlo (MCMC), MinuitFit, or another algorithm.
  virtual void RunMCMC() = 0;

  /// @brief Calculates the required time for each sample or covariance object in a drag race simulation. Inspired by Dan's feature
  /// @param NLaps number of laps, every part of Fitter will be tested with given number of laps and you will get total and average time
  void DragRace(const int NLaps = 100);

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

  /// @brief Perform a 1D sigma var for all samples.
  /// @warning Code uses SampleHandlerFD
  void RunSigmaVarFD();

  /// @brief Allow to start from previous fit/chain
  /// @param FitName Name of previous chain
  /// @todo implement some check that number of params matches etc
  virtual void StartFromPreviousFit(const std::string& FitName);

  /// @brief Get name of class
  virtual inline std::string GetName()const {return "FitterBase";};
 protected:
  /// @brief Process MCMC output
  void ProcessMCMC();

  /// @brief Prepare the output file.
  void PrepareOutput();

  /// @brief Save output and close files.
  void SaveOutput();

  /// @brief Remove obsolete memory and make other checks before fit starts
  /// @todo consider expanding into ParmaterHandler and add more sanitisers
  void SanitiseInputs();

  /// @brief Save the settings that the MCMC was run with.
  void SaveSettings();

  /// @brief YSP: Set up a mapping to store parameters with user-specified ranges, suggested by D. Barrow
  bool GetScaneRange(std::map<std::string, std::vector<double>>& scanRanges);

  /// @brief KS: Check whether we want to skip parameter using skip vector
  bool CheckSkipParameter(const std::vector<std::string>& SkipVector, const std::string& ParamName) const;


  /// @brief For comparison with P-Theta we usually have to apply different parameter values then usual 1, 3 sigma
  ///
  /// Example YAML format:
  /// @code{.yaml}
  /// SigmaVar:
  ///   Q2_norm_7:
  ///     "3": 2.0
  ///   SRC_Norm_O:
  ///     "-1": 0.5
  ///     "1": 1.5
  ///     "3": 2.0
  /// @endcode
  void CustomRange(const std::string& ParName, const double sigma, double& ParamShiftValue);

  /// The manager
  manager *fitMan;

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
  /// step start if restarting
  int stepStart;

  /// store the llh breakdowns
  std::vector<double> sample_llh;
  /// systematic llh breakdowns
  std::vector<double> syst_llh;

  /// Sample holder
  std::vector<SampleHandlerBase*> samples;
  /// Total number of samples used
  unsigned int TotalNSamples;

  /// Systematic holder
  std::vector<ParameterHandlerBase*> systematics;

  /// tells global time how long fit took
  std::unique_ptr<TStopwatch> clock;
  /// tells how long single step/fit iteration took
  std::unique_ptr<TStopwatch> stepClock;
  /// Time of single step
  double stepTime;

  /// Random number
  std::unique_ptr<TRandom3> random;

  /// Output
  TFile *outputFile;
  /// Output cov folder
  TDirectory *CovFolder;
  /// Output sample folder
  TDirectory *SampleFolder;
  /// Output tree with posteriors
  TTree *outTree;
  /// auto save every N steps
  int auto_save;

  /// Necessary for some fitting algorithms like PSO
  bool fTestLikelihood;

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

