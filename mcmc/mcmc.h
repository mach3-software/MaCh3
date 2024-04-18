#pragma once

#include "FitterBase.h"

class mcmc : public FitterBase {
 public:
  mcmc(manager * const fitMan);
  virtual ~mcmc();

  void runMCMC() override;

  inline void setChainLength(unsigned int L) { chainLength = L; };

  inline void setInitialStepNumber(const unsigned int stepNum = 0){stepStart = stepNum;};
  
  inline std::string GetName()const {return "MCMC";};
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

  bool reject; // Do we reject based on hitting boundaries in systs

  unsigned int chainLength; // number of steps in chain

  // simulated annealing
  bool anneal; 
  double AnnealTemp;

  int stepStart;
};

