#pragma once

#include "Fitters/MR2T2.h"

/*
MPI demonstration of Differential Evolution Markov Chain Monte Carlo with snooker updater and Z-matrix crossover
Based on https://arxiv.org/abs/1005.5256


Basic algorithm:
   1. Wait for all fitters to initialise
   2. Propose intial step as MR2T2 step to get going
   ---
   3. Get the first DEMCMC step through  x_0 + g*(y_0+z_0) + e where
      x_0 is the current position of the chain
      y_0 and z_0 are two other randomly selected chains
      g is a scaling factor typically between 0.4 and 1.0
      e is a small random perturbation typically Gaussian with tiny variance
   4. Get likelihood
   5. Broadcast likelihood to all processes
   6. Accept or reject step based on likelihood ratio
   7. Broadcast that this chain has done a step
   8. Repeat
   ---

*/

class DEMCMC : public MR2T2 {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
    DEMCMC(manager *manager);
   
    /// @brief Destructor
    virtual ~DEMCMC() = default;
 protected:
   void RunMCMC() override;

    /// @brief Propose a step
    void ProposeStep() override;
    /// @brief Step acceptance probability


   double gamma;
   double scaling_matrix;
   int n_params;
   std::vector<int> curr_step_idx;
   std::vector<double> transfer_vec;
   std::vector<double> all_params;
};