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

  /// @brief Perform a 2D likelihood scan.
  /// @warning This operation may take a significant amount of time, especially for complex models.
  void Run2DLLHScan();

  /// @brief Get name of class
  virtual inline std::string GetName()const {return "FitterBase";};
protected:

  /// @brief Prepare the output file.
  void PrepareOutput();

  /// @brief Save output and close files.
  void SaveOutput();

  /// @brief Save the settings that the MCMC was run with.
  void SaveSettings();

  // The manager
  manager *fitMan;

  // current state
  unsigned int step;
  // current likelihood
  double logLCurr;
  // proposed likelihood
  double logLProp;
  // current acceptance prob
  double accProb;
  // counts accepted steps
  int accCount;

  //LLH for samples/syst objects
  // oscillation covariance llh
  double osc_llh;
  // store the llh breakdowns
  double *sample_llh;
  // systematic llh breakdowns
  double *syst_llh;

  // Sample holder
  std::vector<samplePDFBase*> samples;

  // Systematic holder
  std::vector<covarianceBase*> systematics;

  // handles oscillation parameters
  covarianceOsc *osc;
  covarianceOsc *osc2;

  // benchmarking, file IO, debugging  etc
  TStopwatch* clock;
  TStopwatch* stepClock;
  double stepTime;

  // Random number
  TRandom3* random;

  // Output
  TFile *outputFile;
  TDirectory *CovFolder;
  TTree *outTree;
  int auto_save; // auto save every N steps

  //Necessary for some fitting algorithms like PSO
  bool fTestLikelihood;
  bool save_nominal;

  //Checks to not repeat some operations
  bool FileSaved;
  bool SettingsSaved;
  bool OutputPrepared;

  #ifdef DEBUG
  //For debugging
  bool debug;
  std::ofstream debugFile; // Output file
  #endif
};

