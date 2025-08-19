#pragma once

#include "Fitters/FitterBase.h"

/// @brief Implementation of MR2T2 algorithm
/// @author Asher Kaboth
class mcmc : public FitterBase {
 public:
   /// @brief Constructor
   /// @param fitMan A pointer to a manager object, which will handle all settings.
  mcmc(manager * const fitMan);
  /// @brief Destructor
  virtual ~mcmc();

  /// @brief Actual implementation of MCMC fitting algorithm
  void RunMCMC() override;

  /// @brief Set how long chain should be
  inline void setChainLength(unsigned int L) { chainLength = L; };
  
  /// @brief Allow to start from previous fit/chain
  /// @param FitName Name of previous chain
  /// @todo implement some check that number of params matches etc
  void StartFromPreviousFit(const std::string& FitName) override;

  /// @brief Get name of class
  inline std::string GetName() const override {return "MCMC";};
 private:
  /// @brief Propose a step
  inline void ProposeStep();

  /// @brief Do we accept the step
  inline void CheckStep();

  /// @brief Print the progress
  inline void PrintProgress();

  /// Do we reject based on hitting boundaries in systs
  bool reject;
  /// number of steps in chain
  unsigned int chainLength;

  /// simulated annealing
  bool anneal;
  /// simulated annealing temperature
  double AnnealTemp;

  /// Penalty term for proposed step
  std::vector<double> syst_llh_prop;

  /// Sample likelihood for proposed step
  std::vector<double> sample_llh_prop;
};

