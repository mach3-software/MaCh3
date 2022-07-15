#ifndef __MCMC_H__
#define __MCMC_H__

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


class mcmc {

 public:
  mcmc(manager * const fitMan);
  ~mcmc();

  void addSamplePDF(samplePDFBase* sample);
  void addSystObj(covarianceBase* cov);
  void addOscHandler(covarianceOsc* oscf); // change this as covOsc now belongs to covBase
  void addOscHandler(covarianceOsc* osca, covarianceOsc *oscb);

  void runMCMC();

  void setChainLength(int L) { chainLength = L; };
  void setOscOnly(bool yes) { osc_only = yes; };
  // determines whether to force osc1 and osc2 to use same mass hierarchy
  void setEqualMH(bool eq) { equalMH = eq; };
  void setEqualOscPar(int par) { equalOscPar[par] = 1; };

  bool getOscOnly() { return osc_only; };

  // initial parameters
  bool StartingParsLoaded(){return init_pos;};

  void PrintInitialState();

  void setInitialStepNumber(int stepNum = 0){stepStart = stepNum;};
  
 private:
  // Helper function called from constructors
  inline void init(std::string name);

  // Prepare the output file
  inline void PrepareOutput();
  
  // Process MCMC output
  inline void ProcessMCMC();

  // Propose a step
  inline void ProposeStep();

  // Do we accept the step
  inline void CheckStep();

  // Print the progress
  inline void PrintProgress();

  // Load starting positions from the end of a previous chain
  inline void ReadParsFromFile(std::string file);
  // Find starting positions from the end of a previous chain
  inline double FindStartingValue(std::string par_file);

  // Save the output settings and MCMC
  inline void SaveSettings();

  inline void SaveChain();

  // The manager
  manager *fitMan;

  bool reject; // Do we reject based on hitting boundaries in systs
  int accCount; // counts accepted steps

  // handles oscillation parameters
  covarianceOsc *osc;
  covarianceOsc *osc2;

  // Determines whether oscillation parameters are equal (ie. osc2[i] = osc[i])
  std::vector<int> equalOscPar;

  // determines whether to force osc1 and osc2 to use same mass hierarchy
  bool equalMH;
  
  // Sample holder
  std::vector<samplePDFBase*> samples;

  // Systematic holder
  std::vector<covarianceBase*> systematics;

  int chainLength; // number of steps in chain

  // current state
  int step; // current step
  double logLCurr; // current likelihood
  double logLProp; // proposed likelihood

  double osc_llh; // oscillation covariance llh
  double *sample_llh; // store the llh breakdowns
  double *syst_llh; // systematic llh breakdowns

  double accProb; // current acceptance prob

  // benchmarking, file IO, debugging  etc
  TStopwatch clock;
  TStopwatch stepClock;
  bool debug;
  ofstream debugFile; // Output file

  // Random number
  TRandom3* random;      

  // Output
  TFile *outputFile;
  TTree *outTree;
  int auto_save; // auto save every N steps
  bool save_nominal;

  // fit types
  bool osc_only; // don't fit nuisance parameters

  // simulated annealing
  bool anneal; 
  double AnnealTemp;

  // starting positions
  bool init_pos; // have the starting parameters been set manually?
  std::map< TString, double > init_pars;

  double stepTime;
  int stepStart;

};

#endif
