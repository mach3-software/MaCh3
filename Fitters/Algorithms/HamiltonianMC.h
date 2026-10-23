#pragma once

#include "Algorithms/MCMCBase.h"

/// @brief Hamiltonian Monte Carlo.
/// @details The mass matrix is each parameter handler's throw matrix, picking up the adapted one
/// when adaption is on, so DeltaT is in units of sigma rather than of the raw parameter values.
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

    /// @brief Draw a fresh velocity and integrate a full leapfrog trajectory from the current state
    /// @details Sets logLProp to \f$U\f$ at the end of the trajectory, as the base class expects
    void ProposeStep() override;

    /// @brief Step acceptance probability, from the change in the Hamiltonian over the trajectory
    double AcceptanceProbability() override;

    /// @brief Work out the flat parameter ordering and enable gradients in the samples
    /// @details Called once at the start of the fit. The flat ordering used by #velocity and
    /// #gradient is every systematic's parameters concatenated in the order they were added.
    void SetupGradients();

    /// @brief Refresh #gradient, i.e. \f$\nabla U\f$, from the parameter priors and the samples
    /// @note Expects the samples to have been reweighted for the current proposed values already
    void UpdateGradU();

    /// @brief Build every mass matrix block from its parameter handler's throw matrix
    void SetupMassMatrix();

    /// @brief Rebuild one handler's mass matrix block from its current throw matrix
    /// @param s Index into #systematics
    void BuildMassBlock(const size_t s);

    /// @brief Rebuild the blocks of any handler whose throw matrix has been re-adapted
    void RefreshMassMatrix();

    /// @brief Draw a fresh velocity \f$w \sim N(0, M^{-1})\f$
    void ThrowVelocity();

    /// @brief Kinetic energy \f$\frac{1}{2} w^T M w\f$ of the current velocity
    double KineticEnergy();

    /// @brief Kinetic energy, and the momentum \f$p = M w\f$ that goes with the current velocity
    /// @param Momentum Filled with \f$p\f$ in the flat ordering, resized if needed
    /// @details Gradient based stopping rules need the momentum rather than the velocity. It costs
    /// one back substitution on top of the forward one KineticEnergy() already does
    double KineticEnergyAndMomentum(std::vector<double>& Momentum);

    /// @brief Reset the proposed parameter values back to the last accepted ones
    void RestoreCurrentPosition();

    /// @brief Half a leapfrog kick, \f$w \rightarrow w - \frac{h}{2} M^{-1} \nabla U\f$
    /// @param Step The integrator step, negative to integrate backwards
    void UpdateVelocity(const double Step);

    /// @brief A leapfrog drift, \f$x \rightarrow x + h\, w\f$, then reweight the samples
    /// @param Step The integrator step, negative to integrate backwards
    /// @details Parameters which would land out of bounds are reflected back in
    void UpdatePosition(const double Step);

    /// @brief Pin fixed parameters' velocities to zero so they never drift
    void ZeroFixedVelocities();

    /// @brief Reweight every sample for the current proposed parameter values
    /// @details All that is needed after a position change, it drives the oscillator, splines and
    /// functional parameters internally
    void ReweightSamples();

    /// @brief Potential energy \f$U = -\log\mathcal{L}\f$ at the current proposed values
    /// @details Also fills syst_llh and sample_llh, which the output tree writes. Returns
    /// M3::_LARGE_LOGL_ or more when the prior penalty alone is already mad, without evaluating the
    /// samples; what that means for the step is the caller's call
    /// @note Expects the samples to have been reweighted already
    double PotentialEnergy();

    /// @brief Copy the proposed parameter values into the flat ordering from SetupGradients()
    /// @param Out Resized if needed
    void SavePosition(std::vector<double>& Out) const;

    /// @brief Write a flat position back into the parameter handlers
    /// @param In Position in the flat ordering from SetupGradients()
    /// @note Writes in place, the samples and splines cache raw pointers into the proposed values
    void RestorePosition(const std::vector<double>& In);

    /// One parameter handler's contribution to the block diagonal mass matrix
    struct MassBlock {
      /// Where this handler's parameters start in the flat vectors
      int Offset = 0;
      /// How many parameters this handler has
      int NPars = 0;
      /// \f$M^{-1}\f$, the handler's throw matrix with fixed parameters decoupled, row major
      std::vector<double> InvMass;
      /// Value of GetThrowMatrixUpdates() this block was built from
      int BuiltAtUpdate = -1;
      /// Lower triangular \f$L\f$ with \f$L L^T = M^{-1}\f$, row major
      std::vector<double> Chol;
    };

    /// How many times one drift may bounce off a parameter's bounds before we give up on the step
    static constexpr int MaxReflections = 100;

    /// Number of leapfrog iterations per step
    int n_leapfrog;
    /// Leapfrog step size
    double delta_t;

    /// Kinetic energy of the velocity the current trajectory started from
    double ke_init;
    /// Potential energy, i.e. -logL, where the current trajectory started from
    double logL_init;

    /// Velocity \f$w = M^{-1} p\f$ of every parameter, in the flat ordering from SetupGradients()
    std::vector<double> velocity;
    /// \f$\nabla U\f$ for every parameter, in the flat ordering from SetupGradients()
    std::vector<double> gradient;

    /// One mass matrix block per systematic, in the same order as #syst_offset
    std::vector<MassBlock> mass;
    /// Whether each parameter is fixed, cached when the mass matrix is built
    std::vector<char> par_fixed;
    /// Standard normal throws ThrowVelocity() correlates with the Cholesky factor
    std::vector<double> rand_throw;
    /// Scratch space for the triangular solve in KineticEnergy()
    std::vector<double> chol_solve;

    /// Where each systematic's parameters start in #velocity and #gradient
    std::vector<int> syst_offset;
    /// Where each sample's GetGradLikelihood() entries start in #gradient
    std::vector<int> sample_offset;
};
