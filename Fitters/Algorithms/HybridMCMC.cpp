#include "Algorithms/HybridMCMC.h"
#include "Parameters/ParameterHandlerGeneric.h"
#include "Manager/Random.h"

#include <algorithm>
#include <cmath>
#include <regex>

namespace {
  // ***************
  /// @brief Shell style wildcard match, so a config can ask for "TotalEScale*"
  bool NameMatches(const std::string& Pattern, const std::string& Name) {
  // ***************
    if (Pattern == Name) return true;
    if (Pattern.find('*') == std::string::npos) return false;

    std::string Expr;
    Expr.reserve(Pattern.size() * 2);
    for (const char C : Pattern) {
      if (C == '*') {
        Expr += ".*";
      } else if (std::isalnum(static_cast<unsigned char>(C)) || C == '_') {
        Expr += C;
      } else {
        // Anything else could be a regex metacharacter, so escape it
        Expr += '\\';
        Expr += C;
      }
    }
    return std::regex_match(Name, std::regex(Expr));
  }
}

// ***************
HybridMCMC::HybridMCMC(Manager *const FitManager) : NoUTurnMC(FitManager) {
// ***************
    AlgorithmName = "HybridMCMC";

    const auto Hybrid = fitMan->raw()["General"]["MCMC"]["Hybrid"];

    const auto Scheme = GetFromManager<std::string>(Hybrid["Scheme"], "composition", __FILE__, __LINE__);
    if (Scheme == "composition") compose = true;
    else if (Scheme == "mixture") compose = false;
    else {
        MACH3LOG_ERROR("Hybrid Scheme is '{}', which isn't a scheme. Use 'composition' to run both "
                       "kernels every step, or 'mixture' to toss a coin for one of them", Scheme);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    rw_prob = GetFromManager<double>(Hybrid["RandomWalkProb"], 0.5, __FILE__, __LINE__);
    rw_sweeps = GetFromManager<int>(Hybrid["RandomWalkSweeps"], 1, __FILE__, __LINE__);
    rw_names = GetFromManager<std::vector<std::string>>(Hybrid["RandomWalkNames"], {}, __FILE__, __LINE__);
    rw_scale_factor = GetFromManager<double>(Hybrid["RandomWalkScale"], 1.0, __FILE__, __LINE__);
    rw_check_every = GetFromManager<int>(Hybrid["RandomWalkCheckEvery"], 1000, __FILE__, __LINE__);
    rw_tune = GetFromManager<bool>(Hybrid["RandomWalkTuneScale"], true, __FILE__, __LINE__);
    rw_target_acc = GetFromManager<double>(Hybrid["RandomWalkTargetAcceptance"], 0.234, __FILE__, __LINE__);
    rw_dlogl_exponent = GetFromManager<double>(Hybrid["RandomWalkPenaltyExponent"], 1.0, __FILE__, __LINE__);

    if (rw_dlogl_exponent <= 0.) {
        MACH3LOG_ERROR("RandomWalkPenaltyExponent has to be positive, you asked for {}. It is the "
                       "p in dlogL proportional to scale^p: 1 for a V shaped log posterior, 2 for "
                       "a smooth quadratic one", rw_dlogl_exponent);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    if (rw_target_acc <= 0. || rw_target_acc >= 1.) {
        MACH3LOG_ERROR("RandomWalkTargetAcceptance has to be strictly between 0 and 1, you asked "
                       "for {}", rw_target_acc);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (rw_tune && rw_check_every <= 0) {
        MACH3LOG_ERROR("RandomWalkTuneScale is on but RandomWalkCheckEvery is {}, so there is no "
                       "batch to measure an acceptance over. Give it a positive number of steps, "
                       "or turn the tuning off", rw_check_every);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    if (rw_scale_factor <= 0.) {
        MACH3LOG_ERROR("RandomWalkScale has to be positive, you asked for {}", rw_scale_factor);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Functional parameters are the ones this whole class exists for. Oscillation parameters are in
    // by default for the same reason: NuOscillator hands back no derivative, so SampleHandlerBase
    // drops their gradient and HMC would sample them from their prior
    const auto TypeNames = GetFromManager<std::vector<std::string>>(
        Hybrid["RandomWalkTypes"], {"Functional", "Oscillation"}, __FILE__, __LINE__);

    for (const auto& Name : TypeNames) {
        bool Found = false;
        for (int t = 0; t < SystType::kSystTypes; ++t) {
            const auto Type = static_cast<SystType>(t);
            if (SystType_ToString(Type) != Name) continue;
            rw_types.push_back(Type);
            Found = true;
            break;
        }
        if (!Found) {
            MACH3LOG_ERROR("RandomWalkTypes has '{}', which isn't a systematic type. Valid ones are:", Name);
            for (int t = 0; t < SystType::kSystTypes; ++t) {
                MACH3LOG_ERROR("  {}", SystType_ToString(static_cast<SystType>(t)));
            }
            throw MaCh3Exception(__FILE__, __LINE__);
        }
    }

    n_rw_pars = 0;
    n_hmc_pars = 0;
    kernel_step = 0;
    rw_steps = 0;
    rw_accepted = 0;
    hmc_steps = 0;
    hmc_accepted = 0;
    rw_step_scale = 0.;
    rw_steps_at_check = 0;
    rw_accepted_at_check = 0;
    rw_tune_updates = 0;
    rw_tune_frozen = false;
    rw_out_of_bounds = 0;
    rw_out_of_bounds_at_check = 0;
    rw_dlogl_sum = 0.;
    rw_dlogl_n = 0;
    rw_dlogl_sum_at_check = 0.;
    rw_dlogl_n_at_check = 0;

    if (compose) {
        if (rw_sweeps < 1) {
            MACH3LOG_ERROR("RandomWalkSweeps has to be at least 1, you asked for {}. At 0 the "
                           "gradient free parameters would never move", rw_sweeps);
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        if (Hybrid["RandomWalkProb"]) {
            MACH3LOG_WARN("RandomWalkProb is set but the scheme is composition, which runs both "
                          "kernels every step and so never consults it. Ignoring it");
        }
        MACH3LOG_INFO("Hybrid sampler, composition: {} random walk sweep{} on the gradient free "
                      "block then one NUTS step on the rest, every step",
                      rw_sweeps, rw_sweeps == 1 ? "" : "s");
    } else {
        if (rw_prob <= 0. || rw_prob >= 1.) {
            MACH3LOG_ERROR("RandomWalkProb has to be strictly between 0 and 1, you asked for {}. At 0 "
                           "the gradient free parameters never move and at 1 nothing else does", rw_prob);
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        if (rw_sweeps != 1) {
            MACH3LOG_ERROR("RandomWalkSweeps is {} but the scheme is mixture, which runs one kernel "
                           "per step. Set Scheme to composition to use sweeps", rw_sweeps);
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        MACH3LOG_INFO("Hybrid sampler, mixture: {:.0f}% random walk on the gradient free block, "
                      "{:.0f}% NUTS on the rest", 100. * rw_prob, 100. * (1. - rw_prob));
        MACH3LOG_INFO("Each block is frozen on the steps the other one got, so expect an "
                      "autocorrelation time in chain steps inflated by about {:.1f}x. Scheme "
                      "'composition' avoids that", 1. / (1. - rw_prob));
    }
}

// *******************
void HybridMCMC::SetupExternalBlock() {
// *******************
    rw_mask.resize(systematics.size());
    hmc_mask.resize(systematics.size());
    par_external.assign(velocity.size(), 0);

    n_rw_pars = 0;
    n_hmc_pars = 0;
    n_fixed_pars = 0;

    for (size_t s = 0; s < systematics.size(); ++s) {
        const int NPars = systematics[s]->GetNumParams();
        rw_mask[s].assign(static_cast<size_t>(NPars), 0);
        hmc_mask[s].assign(static_cast<size_t>(NPars), 0);

        for (int i = 0; i < NPars; ++i) {
            const bool RandomWalk = InRandomWalkBlock(s, i);
            rw_mask[s][static_cast<size_t>(i)] = RandomWalk ? 1 : 0;
            hmc_mask[s][static_cast<size_t>(i)] = RandomWalk ? 0 : 1;

            // A fixed parameter is in neither block's working set. Counting it would overstate
            // what each kernel is actually moving, which matters when reading the rates later
            const bool Fixed = systematics[s]->IsParameterFixed(i);
            if (RandomWalk) {
                par_external[static_cast<size_t>(syst_offset[s] + i)] = 1;
                if (!Fixed) ++n_rw_pars;
            } else if (!Fixed) {
                ++n_hmc_pars;
            }
            if (Fixed) ++n_fixed_pars;
        }
    }

    if (n_fixed_pars > 0) {
        MACH3LOG_INFO("Hybrid split: {} parameters on the random walk, {} on NUTS, {} fixed and in "
                      "neither", n_rw_pars, n_hmc_pars, n_fixed_pars);
    } else {
        MACH3LOG_INFO("Hybrid split: {} parameters on the random walk, {} on NUTS", n_rw_pars, n_hmc_pars);
    }

    if (n_rw_pars == 0) {
        MACH3LOG_ERROR("No parameter matched RandomWalkTypes or RandomWalkNames, so the random walk "
                       "kernel has nothing to sample and {:.0f}% of steps would do nothing. Run "
                       "plain NUTS instead, or widen the selection", 100. * rw_prob);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (n_hmc_pars == 0) {
        MACH3LOG_ERROR("Every parameter was put on the random walk, so NUTS has nothing to sample. "
                       "Run plain MR2T2 instead, which is the same sampler without the overhead");
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    for (size_t s = 0; s < systematics.size(); ++s) {
        int NThisHandler = 0;
        for (const char Flag : rw_mask[s]) NThisHandler += (Flag != 0);
        if (NThisHandler == 0) continue;
        MACH3LOG_INFO("  {} contributes {} parameters to the random walk", systematics[s]->GetName(),
                      NThisHandler);
        for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
            if (rw_mask[s][static_cast<size_t>(i)] == 0) continue;
            MACH3LOG_DEBUG("    {}", systematics[s]->GetParFancyName(i));
        }
    }

    // Robbins-Monro steers the global step scale on the fraction of steps that were accepted, but
    // that counter is fed by both kernels and NUTS accepts nearly every step it takes. The measured
    // rate can therefore never fall near the target, and the step scale would be driven up without
    // bound until the random walk stops moving altogether. Refuse rather than quietly do that
    for (const auto& Syst : systematics) {
        if (!Syst->GetDoAdaption()) continue;
        const auto* Adaptive = Syst->GetAdaptiveHandler();
        if (Adaptive == nullptr || !Adaptive->GetUseRobbinsMonro()) continue;

        MACH3LOG_ERROR("{} has UseRobbinsMonro on, which the hybrid sampler can't use. Robbins-Monro "
                       "tunes the step scale on the accepted step fraction, and that counter sees "
                       "both kernels. NUTS accepts almost every step, so the rate stays far above "
                       "any sensible target and the step scale is driven up until the random walk "
                       "stops moving", Syst->GetName());
        MACH3LOG_ERROR("Set UseRobbinsMonro to false and tune the step scale by hand against the "
                       "random walk acceptance this class reports at the end of the chain. The "
                       "covariance adaption itself is fine to leave on, it is only the step scale "
                       "feedback which breaks");
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Now the block is counted, the proposal can be scaled to it. 2.38/sqrt(d) is the optimal
    // scaling for a random walk on a d dimensional Gaussian, and d is the number of parameters
    // this kernel moves in one proposal, not the number in the handler they happen to live in
    rw_step_scale = rw_scale_factor * 2.38 / std::sqrt(static_cast<double>(n_rw_pars));

    // Say what is being overridden and by how much, because these came from someone's config
    double ScaleLo = 1e30, ScaleHi = -1e30;
    int NOverridden = 0;
    for (size_t s = 0; s < systematics.size(); ++s) {
        for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
            if (rw_mask[s][static_cast<size_t>(i)] == 0) continue;
            const double Scale = systematics[s]->GetIndivStepScale(i);
            if (Scale == 1.) continue;
            ScaleLo = std::min(ScaleLo, Scale);
            ScaleHi = std::max(ScaleHi, Scale);
            ++NOverridden;
        }
    }
    MACH3LOG_INFO("Random walk proposal: throw matrix times {:.4f} = 2.38/sqrt({}) on every "
                  "parameter of the block", rw_step_scale, n_rw_pars);
    if (rw_tune) {
        MACH3LOG_INFO("  that is a starting point: the scale is tuned onto an acceptance of {:.3f} "
                      "every {} steps, and frozen at the end of warmup. 2.38/sqrt(d) is optimal "
                      "against the posterior covariance, and until adaption has one the throw "
                      "matrix is still the prior, which for a well constrained block is much too "
                      "wide", rw_target_acc, rw_check_every);
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
    if (rw_scale_factor != 1.) {
        MACH3LOG_INFO("  RandomWalkScale is {:.3g}, so this is {:.3g} of the optimal scaling",
                      rw_scale_factor, rw_scale_factor);
    }
    EnforceRandomWalkStepScale();

    // Dual averaging counts the global step, so only about rw_prob of the warmup window is actually
    // spent tuning the step size. Say so rather than let the step size freeze early and silently
    if (adapt_step_size && !compose) {
        const double EffectiveWarmup = (1. - rw_prob) * static_cast<double>(warmup_steps);
        MACH3LOG_WARN("Dual averaging only updates on a NUTS step, so a warmup of {} steps tunes the "
                      "step size on roughly {:.0f} of them. Scale WarmupSteps up by 1/{:.2f} to get "
                      "the warmup you asked for, or use Scheme 'composition', where every step is a "
                      "NUTS step", warmup_steps, EffectiveWarmup, 1. - rw_prob);
    }
}

// *******************
bool HybridMCMC::InRandomWalkBlock(const size_t s, const int i) const {
// *******************
    const std::string Name = systematics[s]->GetParFancyName(i);
    for (const auto& Pattern : rw_names) {
        if (NameMatches(Pattern, Name)) return true;
    }

    if (rw_types.empty()) return false;

    // Only a generic handler knows what type each of its parameters is. Anything else can still be
    // selected by name, which is why that comes first
    const auto* Generic = dynamic_cast<const ParameterHandlerGeneric*>(systematics[s]);
    if (Generic == nullptr) return false;

    const SystType Type = Generic->GetParamType(i);
    for (const SystType Wanted : rw_types) {
        if (Type == Wanted) return true;
    }
    return false;
}

// *******************
void HybridMCMC::ApplyProposalMask(const bool RandomWalk) {
// *******************
    for (size_t s = 0; s < systematics.size(); ++s) {
        systematics[s]->SetProposalMask(RandomWalk ? rw_mask[s] : hmc_mask[s]);
    }
}

// *******************
void HybridMCMC::ProposeRandomWalkStep() {
// *******************
    out_of_bounds = false;
    double llh = 0.;

    // The mask leaves NUTS's block at its current value, so this is a conditional update of the
    // gradient free block and nothing else
    for (size_t s = 0; s < systematics.size(); ++s) {
        systematics[s]->ProposeStep();

        syst_llh[s] = systematics[s]->GetLikelihood();
        llh += syst_llh[s];
    }

    // A proposal outside the bounds costs nothing to reject, so don't pay for the reweight
    if (llh >= M3::_LARGE_LOGL_) {
        out_of_bounds = true;
        for (size_t i = 0; i < samples.size(); ++i) {
            sample_llh[i] = M3::_LARGE_LOGL_;
        }
        logLProp = llh;
        return;
    }

    for (size_t i = 0; i < samples.size(); ++i) {
        samples[i]->Reweight();
    }
    // Separate loop, so every sample has migrated its events before any likelihood is taken
    for (size_t i = 0; i < samples.size(); ++i) {
        sample_llh[i] = samples[i]->GetLikelihood();
        llh += sample_llh[i];
    }

    logLProp = llh;
}

// *******************
double HybridMCMC::RandomWalkAcceptanceProbability() const {
// *******************
    if (anneal) {
        return std::min(1., std::exp(-(logLProp - logLCurr) / std::exp(-static_cast<double>(step) / AnnealTemp)));
    }
    return std::min(1., std::exp(logLCurr - logLProp));
}

// *******************
void HybridMCMC::DoRandomWalkStep() {
// *******************
    ApplyProposalMask(true);

    ProposeRandomWalkStep();
    accProb = RandomWalkAcceptanceProbability();
    ++rw_steps;

    if (out_of_bounds) {
        // Nothing was evaluated, so there is no likelihood difference to record
        ++rw_out_of_bounds;
    } else {
        rw_dlogl_sum += logLProp - logLCurr;
        ++rw_dlogl_n;
        if (IsStepAccepted(accProb)) {
            AcceptStep();
            ++rw_accepted;
        }
    }

    // Leave nothing masked behind us, so anything which proposes outside this kernel is unaffected
    for (size_t s = 0; s < systematics.size(); ++s) {
        systematics[s]->ClearProposalMask();
    }
}

// *******************
void HybridMCMC::EnforceRandomWalkStepScale() {
// *******************
    for (size_t s = 0; s < systematics.size(); ++s) {
        bool HasRandomWalkPars = false;
        for (const char Flag : rw_mask[s]) {
            if (Flag != 0) { HasRandomWalkPars = true; break; }
        }
        if (!HasRandomWalkPars) continue;

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
            if (rw_mask[s][static_cast<size_t>(i)] == 0) continue;
            if (systematics[s]->GetIndivStepScale(i) != 1.) {
                systematics[s]->SetIndivStepScale(i, 1.);
            }
        }
    }
}

// *******************
void HybridMCMC::TuneRandomWalkScale() {
// *******************
    if (rw_check_every <= 0) return;
    if (step == 0 || step % static_cast<unsigned int>(rw_check_every) != 0) return;

    const unsigned int Steps = rw_steps - rw_steps_at_check;
    if (Steps == 0) return;
    const double Rate = static_cast<double>(rw_accepted - rw_accepted_at_check) /
                        static_cast<double>(Steps);
    const double OobFrac = static_cast<double>(rw_out_of_bounds - rw_out_of_bounds_at_check) /
                           static_cast<double>(Steps);
    const unsigned int NEval = rw_dlogl_n - rw_dlogl_n_at_check;
    const double MeanDLogL = (NEval > 0)
        ? (rw_dlogl_sum - rw_dlogl_sum_at_check) / static_cast<double>(NEval) : 0.;
    rw_steps_at_check = rw_steps;
    rw_accepted_at_check = rw_accepted;
    rw_out_of_bounds_at_check = rw_out_of_bounds;
    rw_dlogl_sum_at_check = rw_dlogl_sum;
    rw_dlogl_n_at_check = rw_dlogl_n;

    // Shrinking the scale only helps if the rejections are a likelihood penalty. Say which they
    // are on every batch, so a tuning run that cannot work is obvious while it is happening
    MACH3LOG_INFO("Random walk batch: acceptance {:.4f}, out of bounds {:.4f}, mean dlogL over the "
                  "{} evaluated {:+.4g}", Rate, OobFrac, NEval, MeanDLogL);

    const bool Tuning = rw_tune && step < warmup_steps;

    if (Tuning) {
        ++rw_tune_updates;
        double LogStep = 0.;
        bool Jumped = false;

        // A batch which accepted nothing says the scale is wrong, but the acceptance error has
        // saturated at -rw_target_acc and carries no magnitude at all. The mean dlogL does, so use
        // that instead: at dlogL proportional to scale^p, the scale which brings the penalty to
        // order one is the current one divided by dlogL^(1/p). That lands in range in an update or
        // two where the Robbins-Monro step, whose travel is the sum of its decaying gains, needs
        // tens of them. It is the difference between escaping a stall and never leaving it
        if (Rate < 0.01 && NEval > 0 && std::isfinite(MeanDLogL) && MeanDLogL > 2.) {
            LogStep = -std::log(MeanDLogL) / rw_dlogl_exponent;
            // One jump is still bounded, so a pathological batch cannot send the scale somewhere
            // absurd. Two updates at this bound cover four orders of magnitude, which is plenty
            LogStep = std::max(-std::log(100.), LogStep);
            Jumped = true;
        } else {
            // Robbins-Monro on the log of the scale. The gain has to decay for the sequence to
            // settle rather than ring around the target, and 1/sqrt(n) is the usual choice
            const double Gain = 2. / std::sqrt(static_cast<double>(rw_tune_updates));
            LogStep = Gain * (Rate - rw_target_acc);

            // Everything accepted says the scale is too small but again not by how much, so cap
            // one ordinary update at a factor of two either way
            const double MaxLogStep = std::log(2.);
            LogStep = std::max(-MaxLogStep, std::min(MaxLogStep, LogStep));
        }

        rw_scale_factor *= std::exp(LogStep);
        rw_step_scale = rw_scale_factor * 2.38 / std::sqrt(static_cast<double>(n_rw_pars));
        if (Jumped) {
            MACH3LOG_INFO("Random walk accepted {:.4f} of the last {} steps with a mean dlogL of "
                          "{:.4g}, too wide to learn anything from the acceptance. Jumping the "
                          "scale by {:.4g} to {:.4g}, proposal {:.4g} times the throw matrix",
                          Rate, rw_check_every, MeanDLogL, std::exp(LogStep), rw_scale_factor,
                          rw_step_scale);
        } else {
            MACH3LOG_INFO("Random walk acceptance {:.4f} over the last {} steps, target {:.3f}. "
                          "Scale now {:.4g}, proposal {:.4g} times the throw matrix", Rate,
                          rw_check_every, rw_target_acc, rw_scale_factor, rw_step_scale);
        }
        // EnforceRandomWalkStepScale() pushes the new value into the handlers on the next step
        return;
    }

    if (rw_tune && !rw_tune_frozen) {
        rw_tune_frozen = true;
        MACH3LOG_INFO("Warmup is over, freezing the random walk scale at {:.4g}, proposal {:.4g} "
                      "times the throw matrix. From here the kernel is fixed", rw_scale_factor,
                      rw_step_scale);
    }

    if (Rate < 0.01) {
        MACH3LOG_WARN("Random walk acceptance over the last {} steps is {:.4f}, so the {} gradient "
                      "free parameters have effectively stopped moving", rw_check_every, Rate,
                      n_rw_pars);
        if (rw_tune) {
            MACH3LOG_WARN("  the scale was tuned to {:.4g} during warmup and is now frozen, so "
                          "either the posterior moved after warmup or the throw matrix changed "
                          "under it. Lengthen WarmupSteps, or check what adaption did to the "
                          "throw matrix", rw_scale_factor);
        } else {
            MACH3LOG_WARN("  the proposal is almost certainly too wide for the posterior. Turn "
                          "RandomWalkTuneScale on, or lower RandomWalkScale from {:.3g} by hand",
                          rw_scale_factor);
        }
    } else if (Rate > 0.95) {
        MACH3LOG_WARN("Random walk acceptance over the last {} steps is {:.4f}. Almost nothing is "
                      "being rejected, which means the proposal is far too narrow to explore with. "
                      "Scale is {:.3g}", rw_check_every, Rate, rw_scale_factor);
    } else {
        MACH3LOG_DEBUG("Random walk acceptance over the last {} steps: {:.4f}", rw_check_every, Rate);
    }
}

// *******************
void HybridMCMC::DoStep() {
// *******************
    // Adaption resets the global and per parameter step scales the moment it starts throwing, and
    // it recomputes the global one over every parameter in the handler rather than over this
    // kernel's block. Put both back before anything proposes
    EnforceRandomWalkStepScale();

    if (compose) DoCompositionStep();
    else         DoMixtureStep();

    TuneRandomWalkScale();
}

// *******************
void HybridMCMC::DoCompositionStep() {
// *******************
    kernel_step = 2;

    // Each sweep is a complete, posterior invariant update of the gradient free block, so any
    // number of them in a row is still invariant. They are cheap next to a NUTS trajectory, which
    // is the whole reason for offering more than one
    //
    //HW The global accepted step counter has to come back to where it started. MCMCBase divides it
    //by the number of steps to report an acceptance rate, and with sweeps it could otherwise be
    //incremented several times in one step and run past the chain length. The random walk's own
    //accepts live in rw_accepted. Only the counter is restored: AcceptStep's real work, moving the
    //proposal into the current state, stands
    const int AccBeforeRandomWalk = accCount;
    for (int Sweep = 0; Sweep < rw_sweeps; ++Sweep) {
        DoRandomWalkStep();
    }
    accCount = AccBeforeRandomWalk;

    // A rejected final sweep leaves this set, and NUTS is about to ask whether its own proposal was
    // in bounds. Don't let the random walk answer for it
    out_of_bounds = false;

    // NUTS restores the whole position from the current values, so it starts from whatever the
    // sweeps left behind. Every step is a NUTS step here, so its diagnostics are live throughout
    // and dual averaging gets the entire warmup window
    const int AccBefore = accCount;
    NoUTurnMC::DoStep();
    ++hmc_steps;
    if (accCount != AccBefore) ++hmc_accepted;
}

// *******************
void HybridMCMC::DoMixtureStep() {
// *******************
    const bool RandomWalk = (M3::rand::Uniform() < rw_prob);
    kernel_step = RandomWalk ? 1 : 0;

    if (RandomWalk) {
        //HW None of NUTS's per-step diagnostics mean anything on a step NUTS didn't take. Without
        //this they keep the values from the last trajectory and the tree fills with stale repeats,
        //which quietly biases anything averaging over them. Cut on Kernel == 0 to get the real ones
        tree_depth_step = 0;
        n_leapfrog_step = 0;
        n_divergent_step = 0;
        energy_error_step = 0.;
        hamiltonian_init_step = 0.;

        DoRandomWalkStep();
        return;
    }

    // NUTS restores the whole position from the current values before it does anything, so it picks
    // up whatever the random walk last accepted without any help from us
    const int AccBefore = accCount;
    NoUTurnMC::DoStep();
    ++hmc_steps;
    if (accCount != AccBefore) ++hmc_accepted;
}

// *******************
void HybridMCMC::RegisterDiagnosticBranches() {
// *******************
    NoUTurnMC::RegisterDiagnosticBranches();
    outTree->Branch("Kernel", &kernel_step, "Kernel/I");
}

// *******************
void HybridMCMC::PrintKernelRates() const {
// *******************
    const double RwRate = (rw_steps > 0)
        ? static_cast<double>(rw_accepted) / static_cast<double>(rw_steps) : 0.;
    const double HmcRate = (hmc_steps > 0)
        ? static_cast<double>(hmc_accepted) / static_cast<double>(hmc_steps) : 0.;

    if (compose) MACH3LOG_INFO("Scheme: composition, {} random walk sweep(s) per step", rw_sweeps);
    else         MACH3LOG_INFO("Scheme: mixture, random walk probability {:.2f}", rw_prob);
    MACH3LOG_INFO("Random walk kernel: {} steps on {} parameters, acceptance {:.3f}",
                  rw_steps, n_rw_pars, RwRate);
    if (rw_tune) {
        MACH3LOG_INFO("  scale tuned over {} updates to {:.4g}, proposal {:.4g} times the throw "
                      "matrix. Set RandomWalkScale to that to start there next time",
                      rw_tune_updates, rw_scale_factor, rw_step_scale);
    }
    if (rw_steps > 0) {
        const double OobFrac = static_cast<double>(rw_out_of_bounds) / static_cast<double>(rw_steps);
        const double MeanDLogL = (rw_dlogl_n > 0) ? rw_dlogl_sum / static_cast<double>(rw_dlogl_n) : 0.;
        MACH3LOG_INFO("  {:.3f} of sweeps were rejected on bounds before any likelihood was "
                      "evaluated, mean dlogL over the {} that were {:+.4g}",
                      OobFrac, rw_dlogl_n, MeanDLogL);
        if (OobFrac > 0.5) {
            MACH3LOG_WARN("  most sweeps never reached the Metropolis test, so the proposal width "
                          "is not what is holding this block still. ProposeRandomWalkStep() sums "
                          "GetLikelihood() over the whole handler, so a parameter outside its "
                          "bounds anywhere in xsec_cov rejects the step whatever the 28 shifts did");
        }
    }
    MACH3LOG_INFO("NUTS kernel:        {} steps on {} parameters, acceptance {:.3f}",
                  hmc_steps, n_hmc_pars, HmcRate);

    // The base class divides its max depth count by the chain length, which counts the random walk
    // steps NUTS never saw, so restate it over the steps NUTS actually took
    if (hmc_steps > 0) {
        MACH3LOG_INFO("Of the NUTS steps, {:.1f}% ran out of doublings rather than U-turning",
                      100. * static_cast<double>(n_max_depth_hit) / static_cast<double>(hmc_steps));
    }

    if (rw_steps > 0 && RwRate < 0.05) {
        MACH3LOG_WARN("The random walk accepted {:.1f}% of its steps over the whole chain, so the "
                      "gradient free parameters barely moved. Lower RandomWalkScale, currently "
                      "{:.3g}, or turn adaption on for their handler so the throw matrix tightens "
                      "onto the posterior", 100. * RwRate, rw_scale_factor);
    } else if (rw_steps > 0 && RwRate > 0.95) {
        MACH3LOG_WARN("The random walk accepted {:.1f}% of its steps, which is too many to be "
                      "exploring. Raise RandomWalkScale, currently {:.3g}", 100. * RwRate,
                      rw_scale_factor);
    }
}

// *******************
void HybridMCMC::RunMCMC() {
// *******************
    NoUTurnMC::RunMCMC();
    PrintKernelRates();
}
