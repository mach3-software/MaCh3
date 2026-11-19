// HybridMCMC's random walk by delayed acceptance: each proposal is screened on a cheap surrogate of
// the samples' likelihood, and only one that passes pays for the exact reweight. See
// HybridMCMC::DoDelayedRandomWalkStep() for the acceptance rule.
#include "Algorithms/HybridMCMC.h"
#include "Parameters/ParameterHandlerGeneric.h"

#include <algorithm>
#include <cmath>

// *******************
void HybridMCMC::SetupSurrogate() {
// *******************
    da.On = false;
    if (settings.SurrogateStride == 0 || nuts_only) return;
    if (anneal) {
        MACH3LOG_ERROR("RandomWalkSurrogateStride and annealing can't be combined: the delayed "
                       "acceptance ratio is written for the untempered posterior");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    //HW The reverse surrogate reuses the forward migrations with the anchor's weights, which is only
    //the surrogate centred on the proposal if the weights didn't move: shifts alone
    for (size_t s = 0; s < systematics.size(); ++s) {
        const auto* Generic = dynamic_cast<const ParameterHandlerGeneric*>(systematics[s]);
        for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
            if (!blocks.InRandomWalk(s, i) || systematics[s]->IsParameterFixed(i)) continue;
            if (Generic == nullptr || Generic->GetParamType(i) != SystType::kFunc) {
                MACH3LOG_ERROR("RandomWalkSurrogateStride needs a random walk block of shift (functional) "
                               "parameters only, and {} isn't one", systematics[s]->GetParFancyName(i));
                throw MaCh3Exception(__FILE__, __LINE__);
            }
        }
    }
    for (auto* Sample : samples) {
        if (!Sample->EnableSurrogate(settings.SurrogateStride)) {
            MACH3LOG_WARN("{} can't give a surrogate, so the random walk stays plain Metropolis", Sample->GetName());
            return;
        }
    }
    da.SampleLLHAtAnchor.assign(samples.size(), 0.);
    da.On = true;
    MACH3LOG_INFO("Random walk by delayed acceptance: each proposal screened on a surrogate from every "
                  "{}th event's migrations, and reweighted exactly only if it passes", settings.SurrogateStride);
}

// *******************
void HybridMCMC::AnchorSurrogateAtCurrent() {
// *******************
    if (!da.Checked) CheckSurrogate();
    RestoreCurrentPosition();
    // An exact reweight at the current point: a cached pass, as nothing since the last sweeps moved the shifts
    da.PriorAtAnchor = 0.;
    for (size_t s = 0; s < systematics.size(); ++s) da.PriorAtAnchor += systematics[s]->GetLikelihood();
    ReweightSamplesWithoutGradient();
    for (size_t i = 0; i < samples.size(); ++i) {
        da.SampleLLHAtAnchor[i] = samples[i]->GetLikelihood();
        samples[i]->AnchorSurrogate();
    }
    da.Anchored = true;
}

// *******************
void HybridMCMC::CheckSurrogate() {
// *******************
    //HW Sampling every event, the surrogate is the exact prediction, forwards and backwards, so any
    //difference is a bug: the migrations, the anchor or the reverse. Proposals are drawn as the
    //random walk would, so they are typical moves of the block
    da.Checked = true;
    for (auto* Sample : samples) Sample->EnableSurrogate(1);
    double Worst = 0.;
    constexpr int NChecks = 5;
    for (int c = 0; c < NChecks; ++c) {
        da.Anchored = false;
        AnchorSurrogateAtCurrent();
        ApplyProposalMask(true);
        for (size_t s = 0; s < systematics.size(); ++s) systematics[s]->ProposeStep();
        ClearProposalMasks();
        for (size_t i = 0; i < samples.size(); ++i) {
            samples[i]->ReweightSurrogate();
            const double Surrogate = samples[i]->GetLikelihood();
            samples[i]->SetGradientsPaused(true);
            samples[i]->Reweight();
            samples[i]->SetGradientsPaused(false);
            const double Exact = samples[i]->GetLikelihood();
            const double Back = samples[i]->ReverseSurrogateLikelihood();
            const double Diff = std::max(std::fabs(Surrogate - Exact), std::fabs(Back - da.SampleLLHAtAnchor[i]));
            MACH3LOG_INFO("Surrogate check {}, {}: -lnL exact {:.10g}, surrogate {:.10g}; back at the "
                          "anchor {:.10g} against {:.10g}", c + 1, samples[i]->GetName(), Exact, Surrogate,
                          Back, da.SampleLLHAtAnchor[i]);
            Worst = std::max(Worst, Diff);
        }
    }
    RestoreCurrentPosition();
    for (auto* Sample : samples) Sample->EnableSurrogate(settings.SurrogateStride);
    da.Anchored = false;
    if (Worst > 1e-6) {
        MACH3LOG_ERROR("The surrogate over every event differs from the exact -lnL by up to {:.3g}, so its "
                       "migrations are wrong. Turn RandomWalkSurrogateStride off and report this", Worst);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    MACH3LOG_INFO("Surrogate over every event matches the exact -lnL to {:.3g}, forwards and back", Worst);
}

// *******************
void HybridMCMC::DoDelayedRandomWalkStep() {
// *******************
    out_of_bounds = false;
    for (size_t s = 0; s < systematics.size(); ++s) systematics[s]->ProposeStep();
    const double Prior = EvaluatePriors();
    if (Prior >= M3::_LARGE_LOGL_) {
        out_of_bounds = true;
        SetSampleLikelihoodsOutOfBounds();
        logLProp = Prior;
        accProb = 0.;
        rw_tuner->RecordOutOfBounds();
        return;
    }

    // Stage one, on the surrogate centred on the current point, which is exact there
    for (size_t i = 0; i < samples.size(); ++i) samples[i]->ReweightSurrogate();
    double Forward = Prior - da.PriorAtAnchor;
    for (size_t i = 0; i < samples.size(); ++i) Forward += samples[i]->GetLikelihood() - da.SampleLLHAtAnchor[i];
    const double LogA1Forward = std::min(0., -Forward);
    ++da.Proposed;
    if (!IsStepAccepted(std::exp(LogA1Forward))) {
        accProb = std::exp(LogA1Forward);
        rw_tuner->RecordEvaluated(Forward, false);
        return;
    }
    ++da.Passed;

    // Stage two, exact
    logLProp = AddSampleLikelihoods(Prior);
    // The reverse move's stage one, from the proposal back to here, on the surrogate centred there
    double Reverse = da.PriorAtAnchor - Prior;
    for (size_t i = 0; i < samples.size(); ++i) Reverse += samples[i]->ReverseSurrogateLikelihood() - sample_llh[i];
    const double LogA1Reverse = std::min(0., -Reverse);
    const double DLogL = logLProp - logLCurr;
    const double LogA2 = std::min(0., -DLogL + LogA1Reverse - LogA1Forward);
    accProb = std::exp(LogA1Forward + LogA2);

    const bool Accepted = IsStepAccepted(std::exp(LogA2));
    if (Accepted) {
        AcceptStep();
        ++da.Accepted;
        // The proposal is the current point now, and the samples hold its exact prediction
        da.PriorAtAnchor = Prior;
        for (size_t i = 0; i < samples.size(); ++i) {
            da.SampleLLHAtAnchor[i] = sample_llh[i];
            samples[i]->AnchorSurrogate();
        }
    }
    rw_tuner->RecordEvaluated(DLogL, Accepted);
}
