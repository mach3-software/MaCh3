#pragma once

#include "Fitters/FitterBase.h"

/// @brief Base class for MCMC fitting algorithms
/// @author Asher Kaboth

class MCMCBase : public FitterBase {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
    MCMCBase(manager * const fitMan);
    
    /// @brief Destructor
    virtual ~MCMCBase() = default;

    /// @brief Actual implementation of MCMC fitting algorithm
    virtual void RunMCMC() override;

    /// @brief Allow to start from previous fit/chain
    /// @param FitName Name of previous chain
    /// @todo implement some check that number of params matches etc
    void StartFromPreviousFit(const std::string &FitName) override;

    /// @brief Set how long chain should be
    inline void setChainLength(unsigned int L) { chainLength = L; };

    /// @brief Get name of class
    inline std::string GetName() const { return "MCMC"; };

 protected:
    /// @brief Do a step
    inline void DoMCMCStep();

    /// @brief Propose a step
    virtual void ProposeStep()=0;

    /// @brief Do we accept the step
    virtual double CheckStep()=0;

    inline bool AcceptStep(double acc_prob);

    /// @brief Adaptive MCMC step
    inline void AdaptiveStep();

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
};