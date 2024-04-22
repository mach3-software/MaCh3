#pragma once

#include "FitterBase.h"

class mcmc : public FitterBase {
 public:
   /// @brief Constructor
  mcmc(manager * const fitMan);
  /// @brief Destructor
  virtual ~mcmc();

  /// @brief Actual implementation of MCMC fitting algorithm
  void runMCMC() override;

  /// @brief Set how long chain should be
  inline void setChainLength(unsigned int L) { chainLength = L; };

  /// @brief Set initial step number, used when starting from another chain
  inline void setInitialStepNumber(const unsigned int stepNum = 0){stepStart = stepNum;};
  
  /// @brief Get name of class
  inline std::string GetName()const {return "MCMC";};
 private:

  /// @brief Process MCMC output
  inline void ProcessMCMC();

  /// @brief Propose a step
  inline void ProposeStep();

  /// @brief Do we accept the step
  inline void CheckStep();

  /// @brief Print the progress
  inline void PrintProgress();

  /// @brief Load starting positions from the end of a previous chain
  inline void ReadParsFromFile(std::string file);

  bool reject; // Do we reject based on hitting boundaries in systs

  unsigned int chainLength; // number of steps in chain

  // simulated annealing
  bool anneal; 
  double AnnealTemp;

  int stepStart;
};

