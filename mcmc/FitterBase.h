#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "TTree.h"
#include "TString.h"

#include "TVectorT.h"

#include "samplePDF/samplePDFBase.h"
#include "covariance/covarianceBase.h"
#include "covariance/covarianceOsc.h"
#include "covariance/covarianceXsec.h"

#include "manager/manager.h"
#include "MCMCProcessor.h"

//KS: Joy of forward declaration https://gieseanw.wordpress.com/2018/02/25/the-joys-of-forward-declarations-results-from-the-real-world/
class TRandom3;
class TStopwatch;

class FitterBase {

 public:
  FitterBase(manager * const fitMan);
  virtual ~FitterBase();

  void addSamplePDF(samplePDFBase* sample);
  void addSystObj(covarianceBase* cov);
  void addOscHandler(covarianceOsc* oscf); // change this as covOsc now belongs to covBase
  void addOscHandler(covarianceOsc* osca, covarianceOsc *oscb);

  void PrintInitialState();

  virtual void runMCMC() = 0;
  void RunLLHScan();
  void Run2DLLHScan();
 protected:
  // Prepare the output file
  void PrepareOutput();

  // Save output and close files
  void SaveOutput();

  // Save the output settings and MCMC
  // **********************
  // Save the settings that the MCMC was run with
  void SaveSettings();

  // The manager
  manager *fitMan;

  // current state
  unsigned int step; // current step
  double logLCurr; // current likelihood
  double logLProp; // proposed likelihood
  double accProb; // current acceptance prob
  int accCount; // counts accepted steps

  double osc_llh; // oscillation covariance llh
  double *sample_llh; // store the llh breakdowns
  double *syst_llh; // systematic llh breakdowns

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

  bool save_nominal;
  bool fTestLikelihood;
  //Checks to not repeat some operations
  bool FileSaved;
  bool SettingsSaved;
  bool OutputPrepared;

  #ifdef DEBUG
  bool debug;
  std::ofstream debugFile; // Output file
  #endif
};

