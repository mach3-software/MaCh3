// HybridMCMC's flip moves: a deterministic reflection between separated modes (the theta_23 octant,
// say), accepted or rejected on the full posterior. See HybridMCMC::DoFlipStep().
#include "Algorithms/Hybrid/HybridMCMC.h"
#include "Manager/Random.h"

#include <algorithm>
#include <cmath>

// *******************
void HybridMCMC::SetupFlipMoves() {
// *******************
    flips.Handlers.clear();
    for (size_t s = 0; s < systematics.size(); ++s) {
        systematics[s]->SetFlipsInProposeStep(false);
        if (systematics[s]->HasFlips()) flips.Handlers.push_back(s);
    }
    // FlipMoves defaults to one per step when there is something to flip
    flips.PerStep = settings.FlipMoves;
    if (flips.PerStep < 0) flips.PerStep = flips.Handlers.empty() ? 0 : 1;
    if (flips.PerStep > 0 && flips.Handlers.empty()) {
        MACH3LOG_WARN("FlipMoves is {} but no parameter has a FlipParameter, so there is nothing to flip", flips.PerStep);
        flips.PerStep = 0;
    }
    if (flips.PerStep > 0) {
        MACH3LOG_INFO("{} flip move{} per step: one FlipParameter group's reflection at a time, accepted or "
                      "rejected on the full posterior with every other parameter held",
                      flips.PerStep, flips.PerStep == 1 ? "" : "s");
    } else if (!flips.Handlers.empty()) {
        MACH3LOG_WARN("FlipMoves is 0, so the FlipParameter settings are never proposed");
    }
}

// *******************
void HybridMCMC::DoFlipStep() {
// *******************
    // Every handler back on its current values, so the flip is the only thing that moves
    RestoreCurrentPosition();
    // One handler with flips, chosen uniformly; it in turn picks one of its flip groups
    const size_t NHandlers = flips.Handlers.size();
    const size_t Pick = flips.Handlers[std::min(static_cast<size_t>(M3::rand::Uniform() * static_cast<double>(NHandlers)), NHandlers - 1)];
    const double LogJacobian = systematics[Pick]->ProposeFlips();
    ++flips.Proposed;
    ++flips.BatchProposed;
    // A functional flip with no solution at this point proposes nothing, which is a rejection
    if (!std::isfinite(LogJacobian)) {
        ++flips.OutOfBounds;
        RestoreCurrentPosition();
        return;
    }

    // A reflection can land outside the bounds (sin2th_23 below 0.0224 flips past 1 about 0.5112),
    // which costs nothing to reject
    const double Prior = EvaluatePriors();
    if (Prior >= M3::_LARGE_LOGL_) {
        ++flips.OutOfBounds;
        RestoreCurrentPosition();
        return;
    }
    logLProp = AddSampleLikelihoods(Prior);

    // As for the sweeps, the global accepted step counter isn't this move's to keep
    const int AccBefore = accCount;
    // The deterministic move's volume change, |det J|, multiplies the posterior ratio
    if (IsStepAccepted(std::min(1., std::exp(logLCurr - logLProp + LogJacobian)))) {
        AcceptStep();
        ++flips.Accepted;
        ++flips.BatchAccepted;
    }
    accCount = AccBefore;
    // A rejected flip leaves the proposal on the reflected values. Put it back, so nothing after this
    // reads them
    RestoreCurrentPosition();
}

// *******************
void HybridMCMC::ReportFlipBatch() {
// *******************
    if (flips.PerStep <= 0 || settings.CheckEvery <= 0 || flips.BatchProposed == 0) return;
    if (step == 0 || step % static_cast<unsigned int>(settings.CheckEvery) != 0) return;
    MACH3LOG_INFO("Flip moves: {} of {} accepted over the last {} steps ({} of {} over the chain)",
                  flips.BatchAccepted, flips.BatchProposed, settings.CheckEvery, flips.Accepted, flips.Proposed);
    flips.BatchProposed = 0;
    flips.BatchAccepted = 0;
}
