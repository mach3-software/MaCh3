#pragma once

#include "Algorithms/HamiltonianMC.h"
#include "Algorithms/DualAveraging.h"

/// @brief No-U-Turn Sampler, a Hamiltonian MC which picks its own trajectory length.
/// @details Builds the trajectory by doubling and stops when it starts to curve back on itself, so
/// NLeapFrog is no longer needed. Samples from the trajectory multinomially rather than taking its
/// end point.
/// @author Henry Wallace
///
/// @ingroup FittingAlgorithms
class NoUTurnMC : public HamiltonianMC {
public:
    /// @brief Constructor
    /// @param FitManager A pointer to a manager object, which will handle all settings.
    explicit NoUTurnMC(Manager *const FitManager);

    /// @brief Destructor
    virtual ~NoUTurnMC() = default;

    /// @copydoc MCMCBase::RunMCMC
    void RunMCMC() override;

protected:
    /// @brief Build a trajectory and take the sample it selects
    void DoStep() override;

    /// @brief Evaluate the state the next trajectory will start from
    /// @details Only the initial state, not the trajectory. The base class calls this once before
    /// the chain starts and then seeds logLCurr from logLProp without accepting, so leaving the
    /// proposed values anywhere but the current ones would start the chain inconsistent.
    void ProposeStep() override;

    /// @brief The mean Metropolis acceptance over the last trajectory
    /// @note Multinomial sampling is the accept step, so this is a diagnostic rather than a test
    double AcceptanceProbability() override;

    /// One end of the trajectory, enough to resume integrating from it without a reweight
    struct TrajectoryEnd {
      /// Parameter values, in the flat ordering from SetupGradients()
      std::vector<double> Position;
      /// Velocity at this end
      std::vector<double> Velocity;
      /// \f$\nabla U\f$ at this end
      std::vector<double> Gradient;
    };

    /// A state which could become the next sample in the chain
    struct Candidate {
      /// Parameter values, in the flat ordering from SetupGradients()
      std::vector<double> Position;
      /// \f$U\f$ here
      double Potential = 0.;
      /// \f$\frac{1}{2} w^T M w\f$ here
      double Kinetic = 0.;
      /// Per systematic likelihood breakdown, which the output tree writes
      std::vector<double> SystLlh;
      /// Per sample likelihood breakdown, which the output tree writes
      std::vector<double> SampleLlh;
    };

    /// One level of the doubling. Allocated once and reused, so a step needs no allocation
    struct TreeLevel {
      /// Log of the summed trajectory weight over this subtree, relative to the starting state
      double LogSumWeight = 0.;
      /// \f$\rho = \sum p\f$ over this subtree's states
      std::vector<double> Rho;
      /// Velocity of the first state in this subtree
      std::vector<double> VelBeg;
      /// Velocity of the last state in this subtree
      std::vector<double> VelEnd;
      /// The sample this subtree selected
      Candidate Cand;
    };

    /// @brief Size every trajectory buffer, once, before the chain starts
    void SetupTree();

    /// @brief Refuse the fitter settings NUTS can't honour
    void ValidateCompatibility() const;

    /// @brief Add the per-step sampler diagnostics to the output tree
    void RegisterDiagnosticBranches();

    /// @brief Extend the trajectory by \f$2^{Depth}\f$ leapfrog steps, writing into #tree_level
    /// @param Depth Subtree depth, 0 being a single leapfrog step
    /// @param SignedStep The integrator step, negative to extend backwards in time
    /// @param HamiltonianInit \f$H\f$ where the trajectory started, for the weights and divergences
    /// @return false if the subtree U-turned or diverged, in which case it must be discarded whole
    bool BuildTree(const int Depth, const double SignedStep, const double HamiltonianInit);

    /// @brief Is the trajectory still growing, i.e. has it not yet turned back on itself?
    /// @param VelBeg Velocity at one end, @param VelEnd Velocity at the other
    /// @param Rho \f$\sum p\f$ over the span between them
    /// @details The metric correct criterion is \f$\Delta x \cdot p\f$, not \f$\Delta x \cdot w\f$,
    /// which would weight each parameter by its variance and so depend on the parameter units.
    /// Since \f$\rho \approx M \Delta x / \Delta t\f$, testing \f$w \cdot \rho\f$ gives it.
    bool NoUTurn(const std::vector<double>& VelBeg, const std::vector<double>& VelEnd,
                 const std::vector<double>& Rho) const;

    /// @brief Store the live state as a trajectory end
    void SaveEnd(TrajectoryEnd& End);
    /// @brief Put a stored trajectory end back, ready to integrate onwards from
    /// @note Deliberately does not reweight. The next leapfrog step kicks with the stored gradient
    /// and its drift reweights, so the stale sample weights are never read
    void RestoreEnd(const TrajectoryEnd& End);
    /// @brief Store the live state as a candidate sample
    void SaveCandidate(Candidate& Cand, const double Potential, const double Kinetic);

    /// @brief Put the selected sample back and set what the output tree reads
    /// @param HamiltonianInit \f$H\f$ where the trajectory started
    void WriteBackSelected(const double HamiltonianInit);

    /// @brief \f$\log(e^A + e^B)\f$ without overflowing
    static double LogAddExp(const double A, const double B);

    /// @brief Pick a starting step size and arm the dual averaging
    void InitialiseStepSize();

    /// @brief Double or halve the step size until one leapfrog step crosses 0.5 acceptance
    /// @return A step size in the right ballpark, so the warmup doesn't spend itself getting there
    /// @details Leaves the sampler back where it started, reweighted
    double FindReasonableStepSize();

    /// @brief Take one leapfrog step from a saved state and report how well energy was conserved
    /// @return \f$H_0 - H\f$, or minus infinity if the step was unusable
    double TryStepSize(const TrajectoryEnd& Start, const double StepSize, const double HamiltonianInit);

    /// @brief Feed the last trajectory's acceptance to the dual averaging, or freeze it
    void UpdateStepSize();

    /// @brief Summed throw matrix versions, to notice when a mass matrix block was rebuilt
    int ThrowMatrixVersionSum() const;

    /// Tunes #delta_t during warmup
    DualAveragingStepSize dual_averaging;
    /// Whether to adapt #delta_t at all, or take DeltaT as given
    bool adapt_step_size;
    /// Whether to hunt for a starting step size before the chain begins
    bool find_initial_step_size;
    /// Whether a rebuilt mass matrix restarts the adaptation
    bool restart_on_metric_update;
    /// Steps of adaptation before #delta_t is frozen
    unsigned int warmup_steps;
    /// Set once the warmup is over and #delta_t stops moving
    bool step_size_frozen;
    /// Set when the last ProposeStep() saw a mass matrix block rebuilt
    bool metric_changed;

    /// Most doublings per step, so at most \f$2^{MaxTreeDepth} - 1\f$ leapfrog steps
    int MaxTreeDepth;
    /// A state is divergent once \f$H - H_0\f$ exceeds this
    double MaxDeltaH;

    /// The two ends of the current trajectory
    TrajectoryEnd end_minus, end_plus;
    /// The sample the current trajectory has selected
    Candidate selected;
    /// One entry per doubling level, see TreeLevel
    std::vector<TreeLevel> tree_level;
    /// \f$\sum p\f$ over the whole current trajectory
    std::vector<double> rho;
    /// \f$p = M w\f$ of whichever state was last evaluated
    std::vector<double> momentum;

    /// Summed \f$\min(1, e^{H_0 - H})\f$ over every leapfrog step this trajectory took
    double sum_metro_prob;
    /// Largest \f$H - H_0\f$ seen this trajectory
    double max_energy_error;
    /// \f$H - H_0\f$ at the selected sample
    double energy_error_step;
    /// \f$H\f$ where this trajectory started
    double hamiltonian_init_step;
    /// Doublings this trajectory managed
    int tree_depth_step;
    /// Leapfrog steps this trajectory took, discarded subtrees included
    int n_leapfrog_step;
    /// Divergent states this trajectory hit
    int n_divergent_step;
    /// Steps which ran out of doublings rather than U-turning, over the whole chain
    int n_max_depth_hit;
};
