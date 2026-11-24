// HybridMCMC's random walk kernel: the proposal, its width, and the tuning of the width and of the
// number of sweeps during warmup. See HybridMCMC.h for the class.
#include "Algorithms/Hybrid/HybridMCMC.h"

#include <algorithm>
#include <cmath>

// *******************
double HybridMCMC::RandomWalkStepScale() const {
// *******************
    // 2.38/sqrt(d) is the optimal scaling for a random walk on a d dimensional Gaussian, and d is the
    // number of parameters this kernel moves in one proposal, not the number in the handler they
    // happen to live in
    return rw_tuner->GetFactor() * 2.38 / std::sqrt(static_cast<double>(blocks.NRandomWalk()));
}

// *******************
void HybridMCMC::SetupRandomWalkScale() {
// *******************
    rw_step_scale = RandomWalkStepScale();

    // EnforceRandomWalkStepScale() at the end of this function resets every random walk
    // parameter's own step scale to 1, since the throw matrix already carries the prior widths.
    // Measure first what that will override, and by how much: it came from someone's config
    double ScaleLo = 1e30, ScaleHi = -1e30;
    int NOverridden = 0;
    for (size_t s = 0; s < systematics.size(); ++s) {
        for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
            if (!blocks.InRandomWalk(s, i)) continue;
            const double Scale = systematics[s]->GetIndivStepScale(i);
            if (Scale == 1.) continue;
            ScaleLo = std::min(ScaleLo, Scale);
            ScaleHi = std::max(ScaleHi, Scale);
            ++NOverridden;
        }
    }
    MACH3LOG_INFO("Random walk proposal: throw matrix times {:.4f} = 2.38/sqrt({}) on every "
                  "parameter of the block", rw_step_scale, blocks.NRandomWalk());
    if (settings.TuneScale) {
        MACH3LOG_INFO("  that is a starting point: the scale is tuned onto an acceptance of {:.3f} "
                      "every {} steps, and frozen at the end of warmup. 2.38/sqrt(d) is optimal "
                      "against the posterior covariance, and until adaption has one the throw "
                      "matrix is still the prior, which for a well constrained block is much too "
                      "wide", rw_tuner->GetTarget(), settings.CheckEvery);
    } else {
        MACH3LOG_WARN("  RandomWalkTuneScale is off, so this scale is used as is. It is only the "
                      "right one if the throw matrix is already a posterior covariance");
    }
    if (NOverridden > 0) {
        MACH3LOG_INFO("  that replaces {} per parameter step scale(s) from the config, which ran "
                      "from {:.3g} to {:.3g}. The throw matrix already carries the prior widths "
                      "and correlations, so the block needs no second per parameter factor",
                      NOverridden, ScaleLo, ScaleHi);
    }
    if (rw_tuner->GetFactor() != 1.) {
        MACH3LOG_INFO("  RandomWalkScale is {:.3g}, so this is {:.3g} of the optimal scaling",
                      rw_tuner->GetFactor(), rw_tuner->GetFactor());
    }
    EnforceRandomWalkStepScale();
}

// *******************
void HybridMCMC::EnforceRandomWalkStepScale() {
// *******************
    for (size_t s = 0; s < systematics.size(); ++s) {
        if (!blocks.HandlerHasRandomWalk(s)) continue;

        // Changing the global scale of a handler which also holds NUTS parameters is safe: they are
        // never proposed from the throw matrix. It is also what the adaptive handler's own
        // SetIndivStepScaleForSkippedAdaptParams() expects, since it cancels a global change out of
        // the parameters it skipped by rescaling their per parameter factor
        if (std::fabs(systematics[s]->GetGlobalStepScale() - rw_step_scale) > 1e-12) {
            // Quietly: SetStepScale's own warning compares against 2.38/sqrt(every parameter in the
            // handler), which is the number being corrected here, so it would fire on every call
            systematics[s]->SetStepScale(rw_step_scale, false);
        }
        for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
            if (!blocks.InRandomWalk(s, i)) continue;
            if (systematics[s]->GetIndivStepScale(i) != 1.) systematics[s]->SetIndivStepScale(i, 1.);
        }
    }
}

// *******************
void HybridMCMC::DoRandomWalkStep() {
// *******************
    // The masks make every handler propose only its random walk parameters; NUTS's stay put
    if (da.On) {
        // Delayed acceptance: the surrogate is centred on the current point, so (re)anchor it first
        if (!da.Anchored) AnchorSurrogateAtCurrent();
        ApplyProposalMask(true);
        DoDelayedRandomWalkStep();
        ClearProposalMasks();
        return;
    }

    // Plain Metropolis: propose, evaluate the whole posterior there, accept or reject
    ApplyProposalMask(true);
    ProposeRandomWalkStep();
    accProb = RandomWalkAcceptanceProbability();
    if (out_of_bounds) {
        // Nothing was evaluated, so there is no likelihood difference to record
        rw_tuner->RecordOutOfBounds();
    } else {
        const bool Accepted = IsStepAccepted(accProb);
        if (Accepted) AcceptStep();
        rw_tuner->RecordEvaluated(logLProp - logLCurr, Accepted);
    }
    ClearProposalMasks();
}

// *******************
void HybridMCMC::ProposeRandomWalkStep() {
// *******************
    out_of_bounds = false;
    // The mask leaves NUTS's block at its current value, so this is a conditional update of the
    // gradient free block and nothing else
    for (size_t s = 0; s < systematics.size(); ++s) systematics[s]->ProposeStep();
    const double Prior = EvaluatePriors();

    // A proposal outside the bounds costs nothing to reject, so don't pay for the reweight
    if (Prior >= M3::_LARGE_LOGL_) {
        out_of_bounds = true;
        SetSampleLikelihoodsOutOfBounds();
        logLProp = Prior;
        return;
    }
    logLProp = AddSampleLikelihoods(Prior);
}

// *******************
double HybridMCMC::RandomWalkAcceptanceProbability() const {
// *******************
    // With simulated annealing the likelihood ratio is tempered by a temperature falling with the step
    if (anneal) {
        return std::min(1., std::exp(-(logLProp - logLCurr) / std::exp(-static_cast<double>(step) / AnnealTemp)));
    }
    return std::min(1., std::exp(logLCurr - logLProp));
}

// *******************
void HybridMCMC::TuneRandomWalkScale() {
// *******************
    // Only at the end of each batch of RandomWalkCheckEvery steps
    if (settings.CheckEvery <= 0) return;
    if (step == 0 || step % static_cast<unsigned int>(settings.CheckEvery) != 0) return;
    const RandomWalkTuner::Batch Batch = rw_tuner->CloseBatch();
    if (Batch.Steps == 0) return;

    // Shrinking the scale only helps if the rejections are a likelihood penalty. Say which they
    // are on every batch, so a tuning run that cannot work is obvious while it is happening
    MACH3LOG_INFO("Random walk batch: acceptance {:.4f}, out of bounds {:.4f}, mean dlogL over the "
                  "{} evaluated {:+.4g}", Batch.Rate, Batch.OutOfBounds, Batch.NEvaluated, Batch.MeanDLogL);

    // A fit started from a previous file carries on that file's step count, so the warmup
    // window is measured from where this run started, as NoUTurnMC does for the step size
    if (settings.TuneScale && step < stepStart + warmup_steps) {
        const RandomWalkTuner::Update Kind = rw_tuner->Tune(Batch);
        rw_step_scale = RandomWalkStepScale();
        if (Kind == RandomWalkTuner::Update::Jump) {
            MACH3LOG_INFO("Random walk accepted {:.4f} of the last {} steps with a mean dlogL of "
                          "{:.4g}, too wide to learn anything from the acceptance. Jumping the "
                          "scale by {:.4g} to {:.4g}, proposal {:.4g} times the throw matrix",
                          Batch.Rate, settings.CheckEvery, Batch.MeanDLogL, rw_tuner->GetLastChange(),
                          rw_tuner->GetFactor(), rw_step_scale);
        } else {
            MACH3LOG_INFO("Random walk acceptance {:.4f} over the last {} steps, target {:.3f}. "
                          "Scale now {:.4g}, proposal {:.4g} times the throw matrix", Batch.Rate,
                          settings.CheckEvery, rw_tuner->GetTarget(), rw_tuner->GetFactor(), rw_step_scale);
        }
        // EnforceRandomWalkStepScale() pushes the new value into the handlers on the next step
        return;
    }

    // After warmup: announce the freeze once, then only watch for a stalled or overly timid kernel
    if (settings.TuneScale && !rw_tune_frozen) {
        rw_tune_frozen = true;
        MACH3LOG_INFO("Warmup is over, freezing the random walk scale at {:.4g}, proposal {:.4g} "
                      "times the throw matrix. From here the kernel is fixed", rw_tuner->GetFactor(),
                      rw_step_scale);
    }
    WarnAboutBatch(Batch);
}

// *******************
void HybridMCMC::WarnAboutBatch(const RandomWalkTuner::Batch& Batch) const {
// *******************
    if (Batch.Rate < 0.01) {
        MACH3LOG_WARN("Random walk acceptance over the last {} steps is {:.4f}, so the {} gradient "
                      "free parameters have effectively stopped moving", settings.CheckEvery, Batch.Rate,
                      blocks.NRandomWalk());
        if (settings.TuneScale) {
            MACH3LOG_WARN("  the scale was tuned to {:.4g} during warmup and is now frozen, so "
                          "either the posterior moved after warmup or the throw matrix changed "
                          "under it. Lengthen WarmupSteps, or check what adaption did to the "
                          "throw matrix", rw_tuner->GetFactor());
        } else {
            MACH3LOG_WARN("  the proposal is almost certainly too wide for the posterior. Turn "
                          "RandomWalkTuneScale on, or lower RandomWalkScale from {:.3g} by hand",
                          rw_tuner->GetFactor());
        }
    } else if (Batch.Rate > 0.95) {
        MACH3LOG_WARN("Random walk acceptance over the last {} steps is {:.4f}. Almost nothing is "
                      "being rejected, which means the proposal is far too narrow to explore with. "
                      "Scale is {:.3g}", settings.CheckEvery, Batch.Rate, rw_tuner->GetFactor());
    } else {
        MACH3LOG_DEBUG("Random walk acceptance over the last {} steps: {:.4f}", settings.CheckEvery, Batch.Rate);
    }
}

// *******************
void HybridMCMC::TuneSweepCount() {
// *******************
    // Retuned once per batch during warmup, from that batch's measured costs, then fixed
    if (!rw_sweeps_tuning || settings.CheckEvery <= 0) return;
    if (step == 0 || step % static_cast<unsigned int>(settings.CheckEvery) != 0) return;
    if (step >= stepStart + warmup_steps) {
        rw_sweeps_tuning = false;
        MACH3LOG_INFO("Warmup is over, fixing the random walk at {} sweeps per step", rw_sweeps);
        return;
    }
    if (sweep_timing.Sweeps == 0 || sweep_timing.NutsSteps == 0) return;
    // Sweeps cost so the random walk takes TimeShare of each step: n t_sweep = s/(1-s) t_nuts
    const double PerSweep = sweep_timing.SweepSeconds / static_cast<double>(sweep_timing.Sweeps);
    const double PerNuts = sweep_timing.NutsSeconds / static_cast<double>(sweep_timing.NutsSteps);
    const double Wanted = settings.TimeShare / (1. - settings.TimeShare) * PerNuts / PerSweep;
    // At least one sweep, so the block always moves, and at most 500, so a mismeasured
    // batch can't stall a step
    rw_sweeps = std::max(1, std::min(500, static_cast<int>(std::lround(Wanted))));
    MACH3LOG_INFO("Random walk sweeps: {:.3g} ms each against {:.3g} ms per NUTS step, so {} per step for "
                  "{:.0f}% of the time on the random walk", 1e3 * PerSweep, 1e3 * PerNuts, rw_sweeps,
                  100. * settings.TimeShare);
    sweep_timing = SweepTiming{};
}
