#ifndef __MCMC_H__
#define __MCMC_H__

#include "FitterBase.h"

class mcmc : public FitterBase {
 public:
  mcmc(manager * const fitMan);
  ~mcmc();

  void runMCMC();

  void setChainLength(int L) { chainLength = L; };

  // initial parameters
  bool StartingParsLoaded(){return init_pos;};

  void setInitialStepNumber(int stepNum = 0){stepStart = stepNum;};
  
 private:

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

  bool reject; // Do we reject based on hitting boundaries in systs

  int chainLength; // number of steps in chain

  // simulated annealing
  bool anneal; 
  double AnnealTemp;

  // starting positions
  bool init_pos; // have the starting parameters been set manually?
  std::map< TString, double > init_pars;

  int stepStart;
};

#endif
