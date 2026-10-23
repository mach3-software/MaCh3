#pragma once

#include "Algorithms/MCMCBase.h"

/// @brief Hamiltonian Monte Carlo.
/// @author Henry Wallace
///
/// @ingroup FittingAlgorithms
class HamiltonianMC : public MCMCBase {
public:
    /// @brief Constructor
    /// @param FitManager A pointer to a manager object, which will handle all settings.
    HamiltonianMC(Manager *const FitManager);

    /// @brief Destructor
    virtual ~HamiltonianMC() = default;

    /// @copydoc MCMCBase::RunMCMC
    void RunMCMC() override;

protected:
    /// @brief The MCMC step proposal and acceptance
    void DoStep() override;

    /// @brief Draw a fresh momentum and integrate a full leapfrog trajectory from the current state
    /// @details Sets logLProp to \f$U\f$ at the end of the trajectory, as the base class expects
    void ProposeStep() override;

    /// @brief Step acceptance probability, from the change in the Hamiltonian over the trajectory
    double AcceptanceProbability() override;

    /// @brief Work out the flat parameter ordering and enable gradients in the samples
    /// @details Called once at the start of the fit. The flat ordering used by #momentum and
    /// #gradient is every systematic's parameters concatenated in the order they were added.
    void SetupGradients();

    /// @brief Refresh #gradient, i.e. \f$\nabla U\f$, from the parameter priors and the samples
    /// @note Expects the samples to have been reweighted for the current proposed values already
    void UpdateGradU();

    /// @brief Half a leapfrog kick, \f$p \rightarrow p - \frac{\Delta t}{2} \nabla U\f$
    void UpdateMomentum();

    /// @brief One leapfrog drift, \f$x \rightarrow x + \Delta t\, p\f$, then reweight the samples
    /// @details Parameters which would land out of bounds are reflected back in, see the class docs
    void UpdatePosition();

    /// How many times one drift may bounce off a parameter's bounds before we give up on the step
    static constexpr int MaxReflections = 100;

    /// Number of leapfrog iterations per step
    int n_leapfrog;
    /// Leapfrog step size
    double delta_t;

    /// Kinetic energy of the momentum the current trajectory started from
    double ke_init;
    /// Potential energy, i.e. -logL, where the current trajectory started from
    double logL_init;

    /// Momentum conjugate to every parameter, in the flat ordering from SetupGradients()
    std::vector<double> momentum;
    /// \f$\nabla U\f$ for every parameter, in the flat ordering from SetupGradients()
    std::vector<double> gradient;

    /// Where each systematic's parameters start in #momentum and #gradient
    std::vector<int> syst_offset;
    /// Where each sample's GetGradLikelihood() entries start in #gradient
    std::vector<int> sample_offset;
};
