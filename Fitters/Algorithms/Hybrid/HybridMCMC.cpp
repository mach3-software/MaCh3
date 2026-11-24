#include "Algorithms/Hybrid/HybridMCMC.h"
#include "Manager/Random.h"

#include <chrono>

namespace {
  // ***************
  /// @brief General.MCMC.Hybrid, or an empty map: FittingAlgorithm Auto runs with every choice made for it
  YAML::Node HybridSection(Manager* const FitManager) {
  // ***************
    if (CheckNodeExists(FitManager->raw(), "General", "MCMC", "Hybrid")) return FitManager->raw()["General"]["MCMC"]["Hybrid"];
    return YAML::Node(YAML::NodeType::Map);
  }

  // ***************
  /// @brief Whether FittingAlgorithm is Auto
  bool IsAutoMode(Manager* const FitManager) {
  // ***************
    return GetFromManager<std::string>(FitManager->raw()["General"]["FittingAlgorithm"], "", __FILE__, __LINE__) == "Auto";
  }
}

// ***************
HybridMCMC::HybridMCMC(Manager *const FitManager)
    : NoUTurnMC(FitManager),
      settings(HybridSettings::Read(HybridSection(FitManager), IsAutoMode(FitManager), warmup_steps)),
      auto_mode(IsAutoMode(FitManager)),
      rw_sweeps(settings.Sweeps),
      rw_sweeps_tuning(settings.AutoSweeps),
      rw_tuner(std::make_unique<RandomWalkTuner>(settings.ScaleFactor, settings.TargetAcceptance, settings.PenaltyExponent)) {
// ***************
    AlgorithmName = "HybridMCMC";
    //HW Auto picks a NUTS depth limit that leaves trajectories room to U-turn, as Stan's default does
    if (auto_mode && !CheckNodeExists(fitMan->raw(), "General", "MCMC", "MaxTreeDepth")) MaxTreeDepth = 10;
    settings.LogScheme();
}

// *******************
void HybridMCMC::SetupExternalBlock() {
// *******************
    //HW Flips are proposed as a move of their own, see DoFlipStep(), never on the random walk sweeps
    SetupFlipMoves();

    // Decide which parameters each kernel owns. The samples' gradients are already on, so the
    // automatic blocks can ask them which parameters they can differentiate
    blocks.Build(systematics, samples, settings);

    // Tell HamiltonianMC which parameters belong to the other kernel, in its flat ordering, so the
    // metric decouples them and NUTS never moves them
    par_external.assign(velocity.size(), 0);
    for (size_t s = 0; s < systematics.size(); ++s) {
        for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
            if (blocks.InRandomWalk(s, i)) par_external[static_cast<size_t>(syst_offset[s] + i)] = 1;
        }
    }

    if (blocks.NRandomWalk() == 0 && auto_mode) {
        //HW Every parameter has a gradient, so there is nothing for the random walk to do
        MACH3LOG_INFO("Every free parameter has a gradient, so Auto runs plain NUTS steps with no random walk");
        nuts_only = true;
        rw_sweeps = 0;
        rw_sweeps_tuning = false;
        return;
    }
    // From here both kernels run, so check the configuration can work with them
    blocks.RequireBothBlocks(settings.RandomWalkProb);
    blocks.LogRandomWalkMembers(systematics);
    RefuseRobbinsMonro();

    // Everything that needs the random walk block's size: its proposal width, the optional
    // delayed acceptance, and the warning about a mixture's shorter effective warmup
    SetupRandomWalkScale();
    SetupSurrogate();
    WarnAboutMixtureWarmup();
}

// *******************
void HybridMCMC::RefuseRobbinsMonro() const {
// *******************
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
}

// *******************
void HybridMCMC::WarnAboutMixtureWarmup() const {
// *******************
    // Dual averaging counts the global step, so only about RandomWalkProb of the warmup window is
    // actually spent tuning the step size. Say so rather than let the step size freeze early and silently
    if (!adapt_step_size || settings.Compose) return;
    const double EffectiveWarmup = (1. - settings.RandomWalkProb) * static_cast<double>(warmup_steps);
    MACH3LOG_WARN("Dual averaging only updates on a NUTS step, so a warmup of {} steps tunes the "
                  "step size on roughly {:.0f} of them. Scale WarmupSteps up by 1/{:.2f} to get "
                  "the warmup you asked for, or use Scheme 'composition', where every step is a "
                  "NUTS step", warmup_steps, EffectiveWarmup, 1. - settings.RandomWalkProb);
}

// *******************
void HybridMCMC::DoStep() {
// *******************
    // Adaption resets the global and per parameter step scales the moment it starts throwing, and
    // it recomputes the global one over every parameter in the handler rather than over this
    // kernel's block. Put both back before anything proposes
    EnforceRandomWalkStepScale();

    // First, so NUTS or the random walk after it leaves the output branches holding their own values
    for (int f = 0; f < flips.PerStep; ++f) DoFlipStep();

    if (settings.Compose) DoCompositionStep();
    else                  DoMixtureStep();

    // Warmup tuning and batch reports. Each does something only every RandomWalkCheckEvery steps
    TuneRandomWalkScale();
    TuneSweepCount();
    ReportFlipBatch();
}

// *******************
void HybridMCMC::DoCompositionStep() {
// *******************
    kernel_step = 2;

    // Each sweep is a complete, posterior invariant update of the gradient free block, so any
    // number of them in a row is still invariant
    //
    //HW The global accepted step counter has to come back to where it started. MCMCBase divides it
    //by the number of steps to report an acceptance rate, and with sweeps it could otherwise be
    //incremented several times in one step and run past the chain length. The random walk's own
    //accepts are counted in rw_tuner. Only the counter is restored: AcceptStep's real work, moving the
    //proposal into the current state, stands
    const int AccBeforeRandomWalk = accCount;
    const auto SweepsStart = std::chrono::steady_clock::now();
    // NUTS and the flips have moved the rest since the last sweeps, so the surrogate is re-anchored
    da.Anchored = false;
    for (int Sweep = 0; Sweep < rw_sweeps; ++Sweep) {
        DoRandomWalkStep();
    }
    const auto SweepsEnd = std::chrono::steady_clock::now();
    accCount = AccBeforeRandomWalk;

    // A rejected final sweep leaves this set, and NUTS is about to ask whether its own proposal was
    // in bounds. Don't let the random walk answer for it
    out_of_bounds = false;

    // NUTS restores the whole position from the current values, so it starts from whatever the
    // sweeps left behind
    DoNutsStep();

    // Time the two halves, so TuneSweepCount() can balance them during warmup
    sweep_timing.SweepSeconds += std::chrono::duration<double>(SweepsEnd - SweepsStart).count();
    sweep_timing.Sweeps += static_cast<unsigned int>(rw_sweeps);
    sweep_timing.NutsSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - SweepsEnd).count();
    ++sweep_timing.NutsSteps;
}

// *******************
void HybridMCMC::DoMixtureStep() {
// *******************
    const bool RandomWalk = (M3::rand::Uniform() < settings.RandomWalkProb);
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

        da.Anchored = false;
        DoRandomWalkStep();
        return;
    }
    // NUTS restores the whole position from the current values before it does anything, so it picks
    // up whatever the random walk last accepted
    DoNutsStep();
}

// *******************
void HybridMCMC::DoNutsStep() {
// *******************
    // NUTS accepts through the base class's counter, so a change in it says this step moved
    const int AccBefore = accCount;
    NoUTurnMC::DoStep();
    ++hmc_steps;
    if (accCount != AccBefore) ++hmc_accepted;
}

// *******************
void HybridMCMC::ApplyProposalMask(const bool RandomWalk) {
// *******************
    for (size_t s = 0; s < systematics.size(); ++s) {
        systematics[s]->SetProposalMask(blocks.Mask(s, RandomWalk));
    }
}

// *******************
void HybridMCMC::ClearProposalMasks() {
// *******************
    for (size_t s = 0; s < systematics.size(); ++s) {
        systematics[s]->ClearProposalMask();
    }
}

// *******************
double HybridMCMC::EvaluatePriors() {
// *******************
    // Summed in handler order, which is the order every other prior sum in the fitter uses
    double LogL = 0.;
    for (size_t s = 0; s < systematics.size(); ++s) {
        syst_llh[s] = systematics[s]->GetLikelihood();
        LogL += syst_llh[s];
    }
    return LogL;
}

// *******************
void HybridMCMC::ReweightSamplesWithoutGradient() {
// *******************
    for (size_t i = 0; i < samples.size(); ++i) {
        samples[i]->SetGradientsPaused(true);
        samples[i]->Reweight();
        samples[i]->SetGradientsPaused(false);
    }
}

// *******************
double HybridMCMC::AddSampleLikelihoods(double LogL) {
// *******************
    ReweightSamplesWithoutGradient();
    // Separate loop, so every sample has migrated its events before any likelihood is taken
    for (size_t i = 0; i < samples.size(); ++i) {
        sample_llh[i] = samples[i]->GetLikelihood();
        LogL += sample_llh[i];
    }
    return LogL;
}

// *******************
void HybridMCMC::SetSampleLikelihoodsOutOfBounds() {
// *******************
    // What the output tree records for a proposal rejected on bounds before any reweight
    for (size_t i = 0; i < samples.size(); ++i) {
        sample_llh[i] = M3::_LARGE_LOGL_;
    }
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
    const unsigned int rw_steps = rw_tuner->GetSteps();
    const double RwRate = (rw_steps > 0)
        ? static_cast<double>(rw_tuner->GetAccepted()) / static_cast<double>(rw_steps) : 0.;
    const double HmcRate = (hmc_steps > 0)
        ? static_cast<double>(hmc_accepted) / static_cast<double>(hmc_steps) : 0.;

    if (settings.Compose) MACH3LOG_INFO("Scheme: composition, {} random walk sweep(s) per step", rw_sweeps);
    else                  MACH3LOG_INFO("Scheme: mixture, random walk probability {:.2f}", settings.RandomWalkProb);
    MACH3LOG_INFO("Random walk kernel: {} steps on {} parameters, acceptance {:.3f}",
                  rw_steps, blocks.NRandomWalk(), RwRate);
    if (settings.TuneScale) {
        MACH3LOG_INFO("  scale tuned over {} updates to {:.4g}, proposal {:.4g} times the throw "
                      "matrix. Set RandomWalkScale to that to start there next time",
                      rw_tuner->GetUpdates(), rw_tuner->GetFactor(), rw_step_scale);
    }
    if (rw_steps > 0) {
        const double OobFrac = static_cast<double>(rw_tuner->GetOutOfBounds()) / static_cast<double>(rw_steps);
        MACH3LOG_INFO("  {:.3f} of sweeps were rejected on bounds before any likelihood was "
                      "evaluated, mean dlogL over the {} that were {:+.4g}",
                      OobFrac, rw_tuner->GetNEvaluated(), rw_tuner->GetMeanDLogL());
        if (OobFrac > 0.5) {
            MACH3LOG_WARN("  most sweeps never reached the Metropolis test, so the proposal width "
                          "is not what is holding this block still. ProposeRandomWalkStep() sums "
                          "the prior over each whole handler, so a parameter outside its bounds "
                          "anywhere in it rejects the step whatever the random walk moved");
        }
    }
    if (da.On && da.Proposed > 0) {
        MACH3LOG_INFO("  delayed acceptance: {} screened on the surrogate, {} ({:.3f}) passed to the exact "
                      "reweight, {} ({:.3f} of those) accepted there", da.Proposed, da.Passed,
                      static_cast<double>(da.Passed) / static_cast<double>(da.Proposed), da.Accepted,
                      da.Passed > 0 ? static_cast<double>(da.Accepted) / static_cast<double>(da.Passed) : 0.);
    }
    MACH3LOG_INFO("NUTS kernel:        {} steps on {} parameters, acceptance {:.3f}",
                  hmc_steps, blocks.NNuts(), HmcRate);
    if (flips.Proposed > 0) {
        MACH3LOG_INFO("Flip moves:         {} proposed, {} accepted ({:.3f}), {} rejected on bounds", flips.Proposed,
                      flips.Accepted, static_cast<double>(flips.Accepted) / static_cast<double>(flips.Proposed),
                      flips.OutOfBounds);
    }

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
                      "onto the posterior", 100. * RwRate, rw_tuner->GetFactor());
    } else if (rw_steps > 0 && RwRate > 0.95) {
        MACH3LOG_WARN("The random walk accepted {:.1f}% of its steps, which is too many to be "
                      "exploring. Raise RandomWalkScale, currently {:.3g}", 100. * RwRate,
                      rw_tuner->GetFactor());
    }
}

// *******************
void HybridMCMC::RunMCMC() {
// *******************
    NoUTurnMC::RunMCMC();
    PrintKernelRates();
}
