#pragma once

#include "Algorithms/NoUTurnMC.h"
#include "Algorithms/HybridSettings.h"
#include "Algorithms/HybridBlocks.h"
#include "Algorithms/RandomWalkTuner.h"

#include <memory>

/// @brief NUTS on the parameters with a likelihood gradient, a correlated random walk on the rest
///
/// @details Some parameters have no usable likelihood gradient. A functional parameter which
/// shifts reconstructed energy is the clearest case: the shift moves events across bin edges, so
/// a binned likelihood is piecewise constant in it, and the derivative is zero almost everywhere
/// and undefined on the edges. Oscillation parameters and unbinned splines can land in the same
/// place for a more ordinary reason: nothing hands back a derivative.
///
/// Plain HMC over such a parameter is wrong, not just slow: the velocity draw moves it with no
/// restoring force, so what gets explored is its prior, not its posterior. So the parameters are
/// split in two (HybridBlocks) and each block gets its own kernel:
/// - the gradient free block is updated by a correlated random walk with a Metropolis accept;
/// - everything else is updated by NUTS, exactly as NoUTurnMC does.
///
/// Each is an update of its own block with the other held, so each leaves the posterior invariant.
/// HybridSettings::Compose chooses how they are put together:
/// - **composition** (the default) runs HybridSettings::Sweeps random walk sweeps and then one NUTS
///   step, every step. Neither block is ever frozen, dual averaging sees every step, and since a
///   sweep costs one likelihood evaluation against up to \f$2^{\rm MaxTreeDepth}\f$ for NUTS, extra
///   sweeps are nearly free.
/// - **mixture** tosses a coin with state independent weight HybridSettings::RandomWalkProb and runs
///   one kernel. Each block is frozen on the other's steps, which inflates autocorrelation times.
///
/// Optional extra moves, each also posterior invariant:
/// - flip moves (DoFlipStep()): a reflection between separated modes, such as the theta_23 octant;
/// - delayed acceptance of the random walk (DoDelayedRandomWalkStep()), screening each proposal on a
///   cheap surrogate first.
///
/// This class owns the random walk's proposal width outright (EnforceRandomWalkStepScale()): under
/// this sampler only the random walk proposes from the throw matrix, and the right width,
/// \f$2.38/\sqrt{n_{rw}}\f$, depends on how many parameters it moves, not on how many its handler
/// holds. During warmup the width is tuned onto a target acceptance (RandomWalkTuner), and the number
/// of sweeps onto a share of each step's time; both are frozen when warmup ends.
///
/// The code is split by topic: HybridMCMC.cpp (set-up and the schemes), HybridRandomWalk.cpp,
/// HybridDelayedAcceptance.cpp and HybridFlipMoves.cpp.
///
/// @author Henry Wallace
///
/// @ingroup FittingAlgorithms
class HybridMCMC : public NoUTurnMC {
public:
    /// @brief Constructor
    /// @param FitManager A pointer to a manager object, which will handle all settings.
    HybridMCMC(Manager *const FitManager);

    /// @brief Destructor
    virtual ~HybridMCMC() = default;

    /// @copydoc MCMCBase::RunMCMC
    void RunMCMC() override;

protected:
    /// @name Set-up and the two schemes (HybridMCMC.cpp)
    /// @{
    /// @brief Work out the blocks and everything that depends on them
    /// @details Runs through the hook at the top of HamiltonianMC::SetupMassMatrix(), so the metric
    /// already knows to decouple the random walk block by the time it is built
    void SetupExternalBlock() override;
    /// @brief Refuse adaptive handlers with Robbins-Monro on, which would see both kernels' accepts
    void RefuseRobbinsMonro() const;
    /// @brief Under mixture, warn that dual averaging only sees the NUTS steps of the warmup
    void WarnAboutMixtureWarmup() const;

    /// @brief Run a step of whichever scheme was configured, plus any flip moves and tuning
    void DoStep() override;
    /// @brief Run the random walk sweeps and then one NUTS step
    void DoCompositionStep();
    /// @brief Toss the coin and run one of the two kernels
    void DoMixtureStep();
    /// @brief One NUTS step on its block, counted for the end of run summary
    void DoNutsStep();

    /// @brief Point every handler's proposal mask at one of the two blocks
    /// @param RandomWalk True for the random walk's block, false for NUTS's
    void ApplyProposalMask(const bool RandomWalk);
    /// @brief Leave no proposal mask behind, so nothing outside a kernel is affected
    void ClearProposalMasks();

    /// @brief The prior -lnL of every handler at its proposed values, each stored in syst_llh
    double EvaluatePriors();
    /// @brief Reweight every sample at the proposed values, without its gradient
    /// @details The random walk only needs the likelihood. NUTS reweights at its own starting point
    /// before it asks for a gradient (NoUTurnMC::ProposeStep), so leaving it out here is safe
    void ReweightSamplesWithoutGradient();
    /// @brief Reweight every sample (ReweightSamplesWithoutGradient()) and add their -lnL to a total
    /// @param LogL What to add to, the prior part
    /// @return LogL plus every sample's -lnL, each also stored in sample_llh
    double AddSampleLikelihoods(double LogL);
    /// @brief Mark every sample's -lnL as out of bounds, for a proposal that is never reweighted
    void SetSampleLikelihoodsOutOfBounds();

    /// @brief Add NUTS's own diagnostics, plus a branch saying which kernel ran each step
    void RegisterDiagnosticBranches() override;
    /// @brief Report each kernel's acceptance rate separately
    void PrintKernelRates() const;
    /// @}

    /// @name Random walk kernel (HybridRandomWalk.cpp)
    /// @{
    /// @brief Set the starting proposal width for the block and log how it was chosen
    void SetupRandomWalkScale();
    /// @brief The global step scale the block wants: the tuner's factor times \f$2.38/\sqrt{n_{rw}}\f$
    double RandomWalkStepScale() const;
    /// @brief One random walk sweep on the gradient free block, proposal through to acceptance
    void DoRandomWalkStep();
    /// @brief Propose the gradient free block and evaluate the posterior there
    /// @details Sets logLProp to \f$-\log\mathcal{L}\f$ of the whole model, with NUTS's block left at
    /// its current values, which is what makes this a conditional update
    void ProposeRandomWalkStep();
    /// @brief Metropolis acceptance probability for the random walk step
    /// @details The proposal is a Gaussian with the block's own marginal covariance, so it is
    /// symmetric and the ratio is just the posterior one. See ParameterHandlerBase::SetProposalMask()
    double RandomWalkAcceptanceProbability() const;
    /// @brief Scale the random walk proposal to the block it actually moves
    /// @details Sets the global step scale of every handler holding random walk parameters to
    /// #rw_step_scale and their per parameter step scales to 1. Runs every step: the adaptive
    /// handler resets both the moment it starts throwing, and that has to be undone before the next
    /// sweep
    void EnforceRandomWalkStepScale();
    /// @brief Every RandomWalkCheckEvery steps: tune the width onto its target acceptance during
    /// warmup (Robbins-Monro on its log, see RandomWalkTuner), and warn about a stalled block after
    void TuneRandomWalkScale();
    /// @brief During warmup, set the number of sweeps so the random walk takes
    /// HybridSettings::TimeShare of each step's time. Fixed after warmup
    void TuneSweepCount();
    /// @brief Warn about a frozen kernel's batch that accepted almost nothing, or almost everything
    void WarnAboutBatch(const RandomWalkTuner::Batch& Batch) const;
    /// @}

    /// @name Delayed acceptance of the random walk (HybridDelayedAcceptance.cpp)
    /// @{
    /// @brief Turn delayed acceptance on, if asked for and every random walk parameter is a shift
    void SetupSurrogate();
    /// @brief The random walk step by delayed acceptance (Christen and Fox 2005): screened on each
    /// sample's surrogate first, and only a proposal that passes pays for the exact reweight
    /// @details With \f$\tilde\pi_x\f$ the surrogate centred on the current point x, stage one
    /// accepts y with \f$a_1(x,y) = \min(1, \tilde\pi_x(y)/\tilde\pi_x(x))\f$, stage two with
    /// \f$\min(1, \pi(y) a_1(y,x) / (\pi(x) a_1(x,y)))\f$, which leaves \f$\pi\f$ exactly invariant
    /// whatever the surrogate. \f$a_1(y,x)\f$ needs the surrogate centred on y, which the same
    /// sampled migrations reversed give (SampleHandlerInterface::ReverseSurrogateLikelihood). Valid
    /// only while the block moves shifts alone, so the event weights are the same at x and y
    void DoDelayedRandomWalkStep();
    /// @brief Make the current point the surrogate's anchor: an exact reweight there and its prior
    /// and sample -lnL
    void AnchorSurrogateAtCurrent();
    /// @brief Before the first anchor: with every event sampled the surrogate must be the exact
    /// -lnL, at a proposal and back at the anchor. Throws if it isn't
    void CheckSurrogate();
    /// @}

    /// @name Flip moves (HybridFlipMoves.cpp)
    /// @{
    /// @brief Find the handlers with flips and settle how many flip moves to run per step
    void SetupFlipMoves();
    /// @brief One flip move: propose one flip group's reflection with everything else held, and
    /// accept or reject it on the full posterior
    /// @details Some modes are separated by a region neither kernel crosses in reasonable time,
    /// the theta_23 octant for one. The reflection about the FlipParameter point jumps straight from
    /// one to the other. The move is its own inverse, and its acceptance carries the
    /// \f$|\det J|\f$ of any FunctionalFlip partner, see ParameterHandlerBase::ProposeFlips(). It is
    /// kept out of the random walk sweeps, where it would drag their acceptance down
    void DoFlipStep();
    /// @brief Every RandomWalkCheckEvery steps, log the flip acceptance since the last report
    void ReportFlipBatch();
    /// @}

    /// The configuration, as read and checked
    const HybridSettings settings;
    /// Whether FittingAlgorithm is Auto: every choice made automatically, and plain NUTS steps when
    /// nothing lacks a gradient rather than refusing
    const bool auto_mode;
    /// Which parameters each kernel owns
    HybridBlocks blocks;
    /// Whether every parameter had a gradient, so no random walk runs (Auto only)
    bool nuts_only = false;

    /// @name Random walk state
    /// @{
    /// Random walk sweeps per composition step. Starts at HybridSettings::Sweeps, then tuned in
    /// warmup when HybridSettings::AutoSweeps
    int rw_sweeps;
    /// Whether #rw_sweeps is still being tuned, see TuneSweepCount()
    bool rw_sweeps_tuning;
    /// The random walk's statistics and width factor, see RandomWalkTuner
    std::unique_ptr<RandomWalkTuner> rw_tuner;
    /// The global step scale the random walk block uses, see RandomWalkStepScale()
    double rw_step_scale = 0.;
    /// Whether tuning has been frozen at the end of warmup, so it is only announced once
    bool rw_tune_frozen = false;
    /// What TuneSweepCount() balances: time on sweeps and on NUTS steps since its last call
    struct SweepTiming {
      /// Seconds spent sweeping
      double SweepSeconds = 0.;
      /// Seconds spent on NUTS steps
      double NutsSeconds = 0.;
      /// Sweeps timed
      unsigned int Sweeps = 0;
      /// NUTS steps timed
      unsigned int NutsSteps = 0;
    } sweep_timing;
    /// @}

    /// @name NUTS bookkeeping
    /// @{
    /// Which kernel ran this step, for the output tree: 0 NUTS only, 1 random walk only, 2 both.
    /// Under composition it is always 2, so NUTS's own diagnostics are valid on every step
    int kernel_step = 0;
    /// NUTS steps taken, for the end of run summary. The random walk's are in #rw_tuner
    unsigned int hmc_steps = 0;
    /// NUTS steps accepted
    unsigned int hmc_accepted = 0;
    /// @}

    /// @name Delayed acceptance state, see DoDelayedRandomWalkStep()
    /// @{
    struct DelayedAcceptance {
      /// Whether delayed acceptance is on
      bool On = false;
      /// Whether CheckSurrogate() has run
      bool Checked = false;
      /// Whether the anchor is the current point
      bool Anchored = false;
      /// The prior -lnL at the anchor
      double PriorAtAnchor = 0.;
      /// Each sample's -lnL at the anchor
      std::vector<double> SampleLLHAtAnchor;
      /// Proposals screened on the surrogate
      unsigned int Proposed = 0;
      /// Proposals that passed to the exact stage
      unsigned int Passed = 0;
      /// Proposals accepted at the exact stage
      unsigned int Accepted = 0;
    } da;
    /// @}

    /// @name Flip move state, see DoFlipStep()
    /// @{
    struct FlipMoves {
      /// Flip moves per step. Defaults to 1 when some handler has a FlipParameter, 0 otherwise
      int PerStep = 0;
      /// Handlers with at least one flip group
      std::vector<size_t> Handlers;
      /// Proposed over the chain
      unsigned int Proposed = 0;
      /// Accepted over the chain
      unsigned int Accepted = 0;
      /// Rejected on bounds without a likelihood evaluation
      unsigned int OutOfBounds = 0;
      /// Proposed since the last ReportFlipBatch()
      unsigned int BatchProposed = 0;
      /// Accepted since the last ReportFlipBatch()
      unsigned int BatchAccepted = 0;
    } flips;
    /// @}
};
