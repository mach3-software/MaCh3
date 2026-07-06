#pragma once

#include "Fitters/FitterBase.h"
#include "Fitters/MulticanonicalMCMCHandler.h"
/// @brief Base class for MCMC fitting algorithms
/// @details Inherits from `FitterBase` and defines the interface for MCMC-based fitting, including chain management and step handling.
/// @author Asher Kaboth
class MCMCBase : public FitterBase {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
    MCMCBase(Manager * const fitMan);
    
    /// @brief Destructor
    virtual ~MCMCBase() = default;

    /// @brief Actual implementation of MCMC fitting algorithm
    void RunMCMC() final;

    /// @brief Allow to start from previous fit/chain
    /// @param FitName Name of previous chain
    void StartFromPreviousFit(const std::string &FitName) final;

    /// @brief Set how long chain should be
    /// @param L new chain length
   void setChainLength(unsigned int L) { chainLength = L; };
 protected:
    /// @brief The full StartStep->DoStep->EndStep chain
    void DoMCMCStep();

    /// @brief Propose a step
    virtual void ProposeStep()=0;

    /// @brief Actions before step proposal [start stopwatch]
    void PreStepProcess();

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
    void AcceptStep();

    /// @brief Adaptive MCMC step
    void AdaptiveStep();

    /// @brief Print the progress
    /// @param StepsPrint whether to print info about accepted steps and -LogL
    void PrintProgress(const bool StepsPrint = true);

    /// multicanonical handler for umbrella sampling    
    std::unique_ptr<MulticanonicalMCMCHandler> multicanonicalHandler;

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

    /// multi-canonical method toggle on/off
    bool multicanonical;
};
