#pragma once

#include "Algorithms/NoUTurnMC.h"
#include "Parameters/ParameterStructs.h"

/// @brief A mixture of NUTS and a random walk, each sampling its own block of the parameters
///
/// @details Some parameters have no usable likelihood gradient. A functional parameter which
/// shifts reconstructed energy is the clearest case: the shift moves events across bin edges, so
/// a binned likelihood is piecewise constant in it, and the derivative is zero almost everywhere
/// and undefined on the edges. That is not a missing implementation, it is the shape of the
/// problem, and no amount of care in the reweight would give HMC a force to integrate. Oscillation
/// parameters and unbinned splines land in the same place for a more ordinary reason: nothing hands
/// back a derivative, and SampleHandlerBase quietly drops their gradient.
///
/// Running plain HMC over such a parameter is not merely inefficient, it is wrong. The velocity
/// draw gives it a velocity whatever the gradient says, so it drifts ballistically with no
/// restoring force while dragging the likelihood of everything else around with it. What the
/// sampler then explores is the prior, not the posterior.
///
/// So the parameter vector is split in two and each block gets the kernel it deserves:
/// - the gradient free block is updated by a correlated random walk with a Metropolis accept,
///   which needs no derivatives at all
/// - everything else is updated by NUTS, exactly as before
///
/// Both are updates of their own block conditional on the other block sitting at its current
/// value, so each leaves the full posterior invariant on its own. There are two ways to put them
/// together, chosen by Scheme:
///
/// - "mixture" tosses a coin and runs one of the two. A mixture of invariant kernels with state
///   independent weights is invariant, which is why #rw_prob must not look at the state.
/// - "composition" runs #rw_sweeps random walk sweeps and then one NUTS step, every step. A
///   composition of invariant kernels is invariant as well, and needs no reversibility of the
///   product, so this is just as correct.
///
/// Composition is usually the better buy, for three reasons. Neither block is ever frozen, so a
/// trace is a trace of the parameter rather than of the coin, and an autocorrelation time measured
/// in chain steps is not inflated by the roughly 1/#rw_prob of steps on which the parameter could
/// not move. It removes the variance the coin injects into how many updates each block actually
/// received. And dual averaging sees every step, so the step size tunes on the whole warmup window
/// instead of the fraction of it that happened to go to NUTS.
///
/// #rw_sweeps is where the real gain is. A random walk step costs one likelihood evaluation while a
/// NUTS trajectory costs up to 2^MaxTreeDepth gradient evaluations, so extra sweeps are nearly
/// free: against a depth of 5, ten sweeps buy ten times the updates on the gradient free block for
/// about a quarter more work per step. Per step composition costs about what both kernels cost
/// together, so a chain of a given length takes longer than the coin would take; what it does not
/// do is waste that time.
///
/// The random walk block is a natural fit for the adaptive MCMC machinery: adaption accumulates a
/// covariance over the whole chain and the throw matrix it produces becomes this kernel's proposal,
/// so the gradient free parameters get a tuned proposal without anyone choosing step sizes by hand.
///
/// This class owns the proposal width outright, see EnforceRandomWalkStepScale(). It can, because
/// under this sampler the random walk is the only thing which ever proposes from the throw matrix:
/// NUTS reads that matrix as a metric and writes its positions straight out of the leapfrog, and
/// the step scales are applied at proposal time, so nothing NUTS does depends on them. The width
/// a random walk wants depends on how many parameters it moves at once, and both of the numbers
/// the framework would otherwise supply are set by a different count: the global step scale from
/// the config is one number for the whole handler, and the adaptive handler recomputes it as
/// \f$2.38/\sqrt{N}\f$ over every parameter in that handler. With 28 gradient free parameters
/// inside a handler of 531, that is 4.4 times too small, and a per parameter step scale written
/// for some other global scale can be out by far more than that in the other direction. Both are
/// replaced here by \f$2.38/\sqrt{n_{rw}}\f$ over the block this kernel actually moves.
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
    /// @brief Run a step of whichever scheme was configured
    void DoStep() override;

    /// @brief Toss the coin and run one of the two kernels
    void DoMixtureStep();

    /// @brief Run #rw_sweeps random walk sweeps and then one NUTS step
    void DoCompositionStep();

    /// @brief Work out which parameters each kernel owns
    /// @details Fills #rw_mask, #hmc_mask and HamiltonianMC::par_external. Called through the
    /// SetupExternalBlock() hook, which runs at the top of HamiltonianMC::SetupMassMatrix(), so
    /// the metric already knows to decouple the random walk block by the time it is built
    void SetupExternalBlock() override;

    /// @brief Whether a parameter belongs to the gradient free block
    /// @param s Index into #systematics
    /// @param i Parameter index within that handler
    bool InRandomWalkBlock(const size_t s, const int i) const;

    /// @brief One random walk step on the gradient free block, proposal through to acceptance
    void DoRandomWalkStep();

    /// @brief Propose the gradient free block and evaluate the posterior there
    /// @details Sets logLProp to \f$-\log\mathcal{L}\f$ of the whole model, with the block NUTS
    /// owns left at its current values, which is what makes this a conditional update
    void ProposeRandomWalkStep();

    /// @brief Metropolis acceptance probability for the random walk step
    /// @details The proposal is a Gaussian with the block's own marginal covariance, so it is
    /// symmetric and the ratio is just the posterior one. See
    /// ParameterHandlerBase::SetProposalMask()
    double RandomWalkAcceptanceProbability() const;

    /// @brief Point every handler's proposal mask at one of the two blocks
    /// @param RandomWalk Whether to hand the handlers the gradient free block or its complement
    void ApplyProposalMask(const bool RandomWalk);

    /// @brief Add NUTS's own diagnostics, plus a branch saying which kernel ran each step
    void RegisterDiagnosticBranches() override;

    /// @brief Report each kernel's acceptance rate separately
    void PrintKernelRates() const;

    /// @brief Scale the random walk proposal to the block it actually moves
    /// @details Sets the global step scale of every handler holding random walk parameters to
    /// #rw_step_scale and their per parameter step scales to 1, so the proposal width is the throw
    /// matrix times the optimal scaling for #n_rw_pars dimensions and nothing else. Cheap and
    /// idempotent, so it runs every step: the adaptive handler resets both of these the moment it
    /// starts throwing, and that reset has to be undone before the next sweep rather than noticed
    /// afterwards
    void EnforceRandomWalkStepScale();

    /// @brief Tune the proposal width on the random walk's own acceptance, a batch at a time
    /// @details 2.38/sqrt(d) is optimal against the *posterior* covariance. During warmup the
    /// throw matrix is usually still the prior one, which for a well constrained block is far
    /// too wide, so a fixed multiple of it cannot be right. This walks the scale to where the
    /// measured acceptance matches #rw_target_acc, by Robbins-Monro on the log of the scale:
    /// \f$\log s \leftarrow \log s + \gamma_n (\alpha_n - \alpha^*)\f$. Working in the log
    /// keeps the scale positive without a clamp, and the gain decays as
    /// \f$1/\sqrt{n}\f$ so the moves get smaller as the estimate gets better.
    ///
    /// Tuning stops at the end of warmup, exactly as dual averaging freezes the step size. After
    /// that the kernel is fixed, so it leaves the posterior invariant outright and none of this
    /// needs a diminishing adaptation argument to be correct. It also warns, whether or not it is
    /// tuning, because the end of run rate is a whole chain average in which a block that froze
    /// for the first half and then recovered looks healthy
    void TuneRandomWalkScale();

    /// Multiplies the \f$2.38/\sqrt{n_{rw}}\f$ optimal scaling. 1 is the theoretical optimum
    /// against the posterior covariance, and the starting point when #rw_tune is on
    double rw_scale_factor;

    /// Whether to walk #rw_scale_factor onto #rw_target_acc during warmup
    bool rw_tune;

    /// Acceptance the tuning aims for. 0.234 is the high dimensional optimum for a random walk
    double rw_target_acc;

    /// The exponent \f$p\f$ in \f$\Delta\log L \propto s^p\f$, used to size the jump that
    /// gets a stalled block back into range. 2 is a smooth quadratic log posterior; 1 is a
    /// V shaped one, which is what a parameter shifting the binning variable gives, since the
    /// binned likelihood is then piecewise constant and its envelope is \f$|x|\f$ rather than
    /// \f$x^2\f$. Measured at 1.18 over 13 batches for the FD energy scale block.
    ///
    /// 1 is the safe default whatever the truth: assuming 1 where the surface is really quadratic
    /// asks for a smaller scale than needed, which overshoots the acceptance upwards. The block
    /// still moves and the next batch pushes back down. Assuming 2 where it is really 1 would
    /// leave it stalled, which is the failure this is here to escape
    double rw_dlogl_exponent;

    /// How many tuning updates have been done, which sets the Robbins-Monro gain
    int rw_tune_updates;

    /// Whether tuning has been frozen at the end of warmup, so it is only announced once
    bool rw_tune_frozen;

    /// The global step scale the random walk block wants, \f$2.38/\sqrt{n_{rw}}\f$ times
    /// #rw_scale_factor. Known only once SetupExternalBlock() has counted the block
    double rw_step_scale;

    /// How often to check the random walk is still moving, in chain steps. 0 turns the check off
    int rw_check_every;

    /// #rw_steps and #rw_accepted as of the last progress check, to get a batch rate out of them
    unsigned int rw_steps_at_check;
    /// @copydoc rw_steps_at_check
    unsigned int rw_accepted_at_check;

    /// Sweeps thrown out before the Metropolis test because something was outside its bounds.
    /// Separating these matters: a bounds rejection does not care how big the proposal was, so a
    /// scale the tuning keeps shrinking will never fix one
    unsigned int rw_out_of_bounds;
    /// @copydoc rw_out_of_bounds
    unsigned int rw_out_of_bounds_at_check;

    /// Running sum of logLProp - logLCurr over the sweeps that got as far as the Metropolis test,
    /// and how many there were. The mean says whether rejections are a smooth likelihood penalty,
    /// which a smaller proposal would cure, or something a step size cannot reach
    double rw_dlogl_sum;
    /// @copydoc rw_dlogl_sum
    unsigned int rw_dlogl_n;
    /// @copydoc rw_dlogl_sum
    double rw_dlogl_sum_at_check;
    /// @copydoc rw_dlogl_sum
    unsigned int rw_dlogl_n_at_check;

    /// Whether to compose the two kernels every step instead of tossing a coin for one
    bool compose;

    /// Probability of choosing the random walk kernel on any given step, mixture scheme only
    double rw_prob;

    /// Random walk sweeps per step, composition scheme only. Each is a full update of the block
    int rw_sweeps;

    /// Systematic types which go to the random walk because no gradient is available for them
    std::vector<SystType> rw_types;

    /// Parameter names, wildcards allowed, forced into the random walk whatever their type
    std::vector<std::string> rw_names;

    /// Per handler, per parameter: non-zero for the parameters the random walk owns
    std::vector<std::vector<char>> rw_mask;

    /// Per handler, per parameter: the complement of #rw_mask, what NUTS owns
    std::vector<std::vector<char>> hmc_mask;

    /// How many parameters each kernel is responsible for, fixed ones excluded
    int n_rw_pars;
    /// @copydoc n_rw_pars
    int n_hmc_pars;
    /// How many parameters are fixed, and so belong to neither kernel's working set
    int n_fixed_pars;

    /// Which kernel ran this step, for the output tree: 0 NUTS only, 1 random walk only, 2 both.
    /// Under composition it is always 2, so NUTS's own diagnostics are valid on every step
    int kernel_step;

    /// Steps and acceptances of each kernel, for the end of run summary
    unsigned int rw_steps;
    /// @copydoc rw_steps
    unsigned int rw_accepted;
    /// @copydoc rw_steps
    unsigned int hmc_steps;
    /// @copydoc rw_steps
    unsigned int hmc_accepted;
};
