#pragma once

#include "Fitters/FitterBase.h"

/// @brief Base class for MCMC fitting algorithms
/// @details Inherits from `FitterBase` and defines the interface for MCMC-based fitting, including chain management and step handling.
/// @author Asher Kaboth
class MCMCBase : public FitterBase {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
    MCMCBase(manager * const fitMan);

    /// @brief Destructor
    virtual ~MCMCBase() = default;

    /// @brief Actual implementation of MCMC fitting algorithm
    void RunMCMC() override;

    /// @brief Allow to start from previous fit/chain
    /// @param FitName Name of previous chain
    /// @todo implement some check that number of params matches etc
    void StartFromPreviousFit(const std::string &FitName) override;

    /// @brief Set how long chain should be
    /// @param L new chain length
    inline void setChainLength(unsigned int L) { chainLength = L; };
 protected:
   void Init();   

   /// @brief The full StartStep->DoStep->EndStep chain
    inline void DoMCMCStep();

    /// @brief Propose a step
    virtual void ProposeStep()=0;

    /// @brief Actions before step proposal [start stopwatch]
    inline void PreStepProcess();

    /// @brief Actions after step proposal [end stopwatch, fill tree]
   void PostStepProcess();

    /// @brief The MCMC step proposal and acceptance
    virtual void DoStep()=0;

    /// @brief Step acceptance probability
    virtual double AcceptanceProbability()=0;

    /// @brief Is step accepted?
    /// @param acc_prob used for telling whether step is accepted or not
    bool IsStepAccepted(const double acc_prob);

    /// @brief Accept a step
    virtual void AcceptStep();

    /// @brief Adaptive MCMC step
    inline void AdaptiveStep();

    /// @brief Print the progress
    virtual void PrintProgress();

    /// Do we reject based on hitting boundaries in systs
    bool out_of_bounds;

    /// Accept
    bool accept;

    /// number of steps in chain
    unsigned int chainLength;

    /// simulated annealing
    bool anneal;
    /// simulated annealing temperature
    double AnnealTemp;

};
