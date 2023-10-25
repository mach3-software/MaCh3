#ifndef __FITTERBASE_H__
#define __FITTERBASE_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "TStopwatch.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom.h"

#include "TVectorT.h"

#include "samplePDF/samplePDFBase.h"
#include "covariance/covarianceBase.h"
#include "covariance/covarianceOsc.h"

#include "manager/manager.h"
#include "MCMCProcessor.h"


class FitterBase {

 public:
  FitterBase(manager * const fitMan);
  ~FitterBase();

  void addSamplePDF(samplePDFBase* sample);
  void addSystObj(covarianceBase* cov);
  void addOscHandler(covarianceOsc* oscf); // change this as covOsc now belongs to covBase
  void addOscHandler(covarianceOsc* osca, covarianceOsc *oscb);

  void setOscOnly(bool yes) { osc_only = yes; };
  // determines whether to force osc1 and osc2 to use same mass hierarchy
  void setEqualMH(bool eq) { equalMH = eq; };
  void setEqualOscPar(int par) { equalOscPar[par] = 1; };

  bool getOscOnly() { return osc_only; };

  void PrintInitialState();

 protected:
  // Prepare the output file
  void PrepareOutput();

  // Save output and close files
  void SaveOutput();

  // Save the output settings and MCMC
  // **********************
  // Save the settings that the MCMC was run with
  inline void SaveSettings() {
    // If we're using the new mcmc constructor which knows about settings
  // **********************
    if (fitMan == NULL) {
      std::cout << "************************" << std::endl;
      std::cout << "************************" << std::endl;
      std::cout << "WARNING WILL NOT SAVE MANAGER OUTPUT TO FILE BECAUSE YOU USED A DEPRECATED CONSTRUCTOR" << std::endl;
      std::cout << "************************" << std::endl;
      std::cout << "************************" << std::endl;
    } else {
    // Save the settings we have in the manager
    //fitMan->SaveSettings(outputFile);
    // Warn if we're running a deprecated constructor (again)
    }
  }

  // The manager
  manager *fitMan;

  // fit types
  bool osc_only; // don't fit nuisance parameters
  // Determines whether oscillation parameters are equal (ie. osc2[i] = osc[i])
  std::vector<int> equalOscPar;
  // determines whether to force osc1 and osc2 to use same mass hierarchy
  bool equalMH;

  // current state
  int step; // current step
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
  bool debug;
  ofstream debugFile; // Output file

  // Random number
  TRandom3* random;

  // Output
  TFile *outputFile;
  TTree *outTree;
  int auto_save; // auto save every N steps

  double stepTime;

  bool save_nominal;
};

#endif
