#include "Algorithms/NoUTurnMC.h"

#include <limits>
// HW No-U-Turn Sampler, on top of the HMC machinery

// ***************
NoUTurnMC::NoUTurnMC(Manager *const FitManager) : HamiltonianMC(FitManager) {
// ***************
    AlgorithmName = "NoUTurnMC";

    MaxTreeDepth = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["MaxTreeDepth"], 5, __FILE__, __LINE__);
    MaxDeltaH = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["MaxDeltaH"], 1000., __FILE__, __LINE__);

    if (MaxTreeDepth < 1 || MaxTreeDepth > 16) {
        MACH3LOG_ERROR("MaxTreeDepth has to be between 1 and 16, you asked for {}", MaxTreeDepth);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (MaxDeltaH <= 0.) {
        MACH3LOG_ERROR("MaxDeltaH has to be positive, you asked for {}", MaxDeltaH);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (MaxDeltaH < 10.) {
        MACH3LOG_WARN("MaxDeltaH = {} is small enough that almost every trajectory will be called "
                      "divergent", MaxDeltaH);
    }
    if (CheckNodeExists(fitMan->raw()["General"]["MCMC"]["NLeapFrog"])) {
        MACH3LOG_WARN("NUTS picks its own trajectory length, NLeapFrog is ignored");
    }

    adapt_step_size = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["DualAveraging"], true, __FILE__, __LINE__);
    find_initial_step_size = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["FindInitialStepSize"], true, __FILE__, __LINE__);
    restart_on_metric_update = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["RestartDualAveragingOnMetricUpdate"], true, __FILE__, __LINE__);
    dual_averaging.Target = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["TargetAcceptance"], 0.8, __FILE__, __LINE__);
    warmup_steps = GetFromManager<unsigned int>(fitMan->raw()["General"]["MCMC"]["WarmupSteps"], chainLength / 2, __FILE__, __LINE__);

    if (dual_averaging.Target <= 0. || dual_averaging.Target >= 1.) {
        MACH3LOG_ERROR("TargetAcceptance has to be in (0, 1), you asked for {}", dual_averaging.Target);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (dual_averaging.Target < 0.5) {
        MACH3LOG_WARN("TargetAcceptance = {} is low for HMC, which wants about 0.65 at best and 0.8 "
                      "to be safe", dual_averaging.Target);
    }
    if (adapt_step_size && warmup_steps < 10) {
        MACH3LOG_ERROR("WarmupSteps has to be at least 10, you asked for {}", warmup_steps);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (adapt_step_size && warmup_steps >= chainLength) {
        MACH3LOG_ERROR("WarmupSteps ({}) has to leave some of the {} step chain for sampling",
                       warmup_steps, chainLength);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    step_size_frozen = false;
    metric_changed = false;
    sum_metro_prob = 0.;
    max_energy_error = 0.;
    energy_error_step = 0.;
    hamiltonian_init_step = 0.;
    tree_depth_step = 0;
    n_leapfrog_step = 0;
    n_divergent_step = 0;
    n_max_depth_hit = 0;

    MACH3LOG_INFO("Running NUTS with max tree depth {}, so up to {} leapfrog steps per step",
                  MaxTreeDepth, (1 << MaxTreeDepth) - 1);
    MACH3LOG_INFO("Divergence threshold MaxDeltaH = {}", MaxDeltaH);
    if (adapt_step_size) {
        MACH3LOG_INFO("Dual averaging on, targeting acceptance {} over {} warmup steps",
                      dual_averaging.Target, warmup_steps);
        MACH3LOG_WARN("The first {} steps are step size warmup. They land in the output tree but are "
                      "NOT valid posterior samples, cut them when analysing the chain", warmup_steps);
    } else {
        MACH3LOG_INFO("Dual averaging off, step size stays at DeltaT = {}", delta_t);
    }
}

// *******************
void NoUTurnMC::RunMCMC() {
// *******************
    // All of this has to happen before the base class proposes its first step
    SetupGradients();
    SetupMassMatrix();
    SetupTree();
    ValidateCompatibility();

    // PrepareOutput is idempotent, so calling it here to get at outTree leaves the base class's own
    // call a no-op. Order matches MCMCBase::RunMCMC so the output file stays comparable
    SaveSettings();
    PrepareOutput();
    RegisterDiagnosticBranches();
    InitialiseStepSize();

    MCMCBase::RunMCMC();

    const double MaxDepthFraction = (chainLength > 0)
        ? static_cast<double>(n_max_depth_hit) / static_cast<double>(chainLength) : 0.;
    MACH3LOG_INFO("NUTS finished, {:.1f}% of steps ran out of doublings rather than U-turning",
                  100. * MaxDepthFraction);
    if (adapt_step_size) {
        MACH3LOG_INFO("Final step size {}, which dual averaging clamped {} times",
                      delta_t, dual_averaging.NClamped);
    }
    if (MaxDepthFraction > 0.2) {
        MACH3LOG_WARN("Trajectories are being cut off by MaxTreeDepth = {}, so the chain is closer "
                      "to a random walk than it should be. Raise it, or find out why the step size "
                      "is so small", MaxTreeDepth);
    }
}

// *******************
void NoUTurnMC::SetupTree() {
// *******************
    const size_t NPars = velocity.size();
    const size_t NSyst = systematics.size();
    const size_t NSamp = samples.size();

    auto SetupCandidate = [&](Candidate& Cand) {
        Cand.Position.assign(NPars, 0.);
        Cand.SystLlh.assign(NSyst, 0.);
        Cand.SampleLlh.assign(NSamp, 0.);
    };
    auto SetupEnd = [&](TrajectoryEnd& End) {
        End.Position.assign(NPars, 0.);
        End.Velocity.assign(NPars, 0.);
        End.Gradient.assign(NPars, 0.);
    };

    SetupEnd(end_minus);
    SetupEnd(end_plus);
    SetupCandidate(selected);

    rho.assign(NPars, 0.);
    momentum.assign(NPars, 0.);

    // One level per doubling, so the recursion never allocates
    tree_level.resize(static_cast<size_t>(MaxTreeDepth) + 1);
    for (size_t Depth = 0; Depth < tree_level.size(); ++Depth) {
        tree_level[Depth].Rho.assign(NPars, 0.);
        tree_level[Depth].VelBeg.assign(NPars, 0.);
        tree_level[Depth].VelEnd.assign(NPars, 0.);
        SetupCandidate(tree_level[Depth].Cand);
    }
}

// *******************
void NoUTurnMC::ValidateCompatibility() const {
// *******************
    // FitterBase only sizes syst_llh and sample_llh on the real likelihood path, so the test
    // likelihood would have us writing past the end of both
    if (fTestLikelihood) {
        MACH3LOG_ERROR("NUTS can't be run against the test likelihood, it has no gradient to follow");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (anneal) {
        MACH3LOG_ERROR("NUTS can't be run with simulated annealing, the Hamiltonian assumes a fixed "
                       "target distribution");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    for (const auto& syst : systematics) {
        auto* Adaptive = syst->GetAdaptiveHandler();
        if (Adaptive != nullptr && Adaptive->GetUseRobbinsMonro()) {
            MACH3LOG_WARN("{} has Robbins-Monro enabled, which tunes the random walk step scale that "
                          "NUTS never reads. Only its adapted matrix matters here", syst->GetName());
        }
    }
}

// *******************
void NoUTurnMC::RegisterDiagnosticBranches() {
// *******************
    outTree->Branch("TreeDepth", &tree_depth_step, "TreeDepth/I");
    outTree->Branch("NLeapfrog", &n_leapfrog_step, "NLeapfrog/I");
    outTree->Branch("NDivergent", &n_divergent_step, "NDivergent/I");
    outTree->Branch("StepSize", &delta_t, "StepSize/D");
    outTree->Branch("EnergyError", &energy_error_step, "EnergyError/D");
    outTree->Branch("EnergyErrorMax", &max_energy_error, "EnergyErrorMax/D");
    outTree->Branch("Energy", &hamiltonian_init_step, "Energy/D");
}

// *******************
// Only the starting state, see the header for why the trajectory isn't built here
void NoUTurnMC::ProposeStep() {
// *******************
    out_of_bounds = false;

    // Noticing a rebuilt mass matrix matters: the step size that suited the old metric won't suit
    // the new one, so the dual averaging has to start again
    const int VersionBefore = ThrowMatrixVersionSum();
    RefreshMassMatrix();
    metric_changed = (ThrowMatrixVersionSum() != VersionBefore);

    RestoreCurrentPosition();
    ReweightSamples();
    UpdateGradU();
    logLProp = PotentialEnergy();
}

// *******************
void NoUTurnMC::DoStep() {
// *******************
    // Puts us back at the last accepted state and gives us U and dU/dx there
    ProposeStep();
    logL_init = logLProp;

    ThrowVelocity();
    ke_init = KineticEnergyAndMomentum(momentum);
    const double HamiltonianInit = logL_init + ke_init;

    // Both ends start at the current state, and so does the selected sample
    SaveEnd(end_minus);
    end_plus = end_minus;
    rho = momentum;
    SaveCandidate(selected, logL_init, ke_init);

    // Weights are relative to the starting state, which keeps them from overflowing
    double LogSumWeight = 0.;
    sum_metro_prob = 0.;
    max_energy_error = 0.;
    n_leapfrog_step = 0;
    n_divergent_step = 0;

    bool Moved = false;
    int Depth = 0;

    for (; Depth < MaxTreeDepth; ++Depth) {
        const double SignedStep = (M3::rand::Uniform() < 0.5) ? -delta_t : delta_t;

        // Resume from whichever end we're growing, no reweight needed, see RestoreEnd
        RestoreEnd((SignedStep < 0.) ? end_minus : end_plus);

        // A subtree which U-turns or diverges is thrown away whole, partially using it would break
        // detailed balance
        if (!BuildTree(Depth, SignedStep, HamiltonianInit)) break;

        const TreeLevel& Sub = tree_level[static_cast<size_t>(Depth)];

        // Biased progressive sampling: deliberately favours the fresh subtree over everything
        // already in the tree, which is what makes NUTS mix better than plain multinomial HMC
        if (Sub.LogSumWeight > LogSumWeight ||
            M3::rand::Uniform() < std::exp(Sub.LogSumWeight - LogSumWeight)) {
            selected = Sub.Cand;
            Moved = true;
        }
        LogSumWeight = LogAddExp(LogSumWeight, Sub.LogSumWeight);

        for (size_t i = 0; i < rho.size(); ++i) {
            rho[i] += Sub.Rho[i];
        }
        SaveEnd((SignedStep < 0.) ? end_minus : end_plus);

        if (!NoUTurn(end_minus.Velocity, end_plus.Velocity, rho)) {
            ++Depth;
            break;
        }
    }

    if (Depth == MaxTreeDepth) ++n_max_depth_hit;
    tree_depth_step = Depth;

    WriteBackSelected(HamiltonianInit);

    // Exactly once per step, and only when the trajectory actually moved us. Unconditionally
    // accepting would peg accCount to the chain length and freeze the adaptive accepted-step counter
    if (Moved) AcceptStep();

    UpdateStepSize();
}

// *******************
bool NoUTurnMC::BuildTree(const int Depth, const double SignedStep, const double HamiltonianInit) {
// *******************
    if (Depth == 0) {
        // One leapfrog step: half kick, drift, half kick
        UpdateVelocity(SignedStep);
        UpdatePosition(SignedStep);
        ++n_leapfrog_step;

        // UpdatePosition gives up by teleporting the parameter back to its current value, which
        // would wreck reversibility mid trajectory. Treat it as a divergence instead
        if (out_of_bounds) {
            out_of_bounds = false;
            ++n_divergent_step;
            return false;
        }
        UpdateGradU();
        UpdateVelocity(SignedStep);

        const double Potential = PotentialEnergy();
        const double Kinetic = KineticEnergyAndMomentum(momentum);
        const double Hamiltonian = Potential + Kinetic;
        const double LogWeight = HamiltonianInit - Hamiltonian;
        const bool Finite = std::isfinite(Hamiltonian);

        // Every leaf counts towards the acceptance statistic, discarded ones included, so that
        // divergences pull the step size down
        sum_metro_prob += Finite ? std::min(1., std::exp(LogWeight)) : 0.;
        max_energy_error = std::max(max_energy_error, Finite ? -LogWeight : MaxDeltaH);

        if (!Finite || -LogWeight > MaxDeltaH) {
            ++n_divergent_step;
            return false;
        }

        TreeLevel& Leaf = tree_level[0];
        Leaf.LogSumWeight = LogWeight;
        Leaf.Rho = momentum;
        Leaf.VelBeg = velocity;
        Leaf.VelEnd = velocity;
        SaveCandidate(Leaf.Cand, Potential, Kinetic);
        return true;
    }

    // Left half, which leaves the live state at its far end, exactly where the right half starts
    if (!BuildTree(Depth - 1, SignedStep, HamiltonianInit)) return false;

    TreeLevel& Level = tree_level[static_cast<size_t>(Depth)];
    const TreeLevel& Sub = tree_level[static_cast<size_t>(Depth - 1)];

    Level.LogSumWeight = Sub.LogSumWeight;
    Level.Rho = Sub.Rho;
    Level.VelBeg = Sub.VelBeg;
    Level.Cand = Sub.Cand;

    if (!BuildTree(Depth - 1, SignedStep, HamiltonianInit)) return false;

    // Uniform progressive sampling within a subtree, unlike the biased rule at the top level
    const double LogSumWeightTotal = LogAddExp(Level.LogSumWeight, Sub.LogSumWeight);
    if (M3::rand::Uniform() < std::exp(Sub.LogSumWeight - LogSumWeightTotal)) {
        Level.Cand = Sub.Cand;
    }
    Level.LogSumWeight = LogSumWeightTotal;
    for (size_t i = 0; i < Level.Rho.size(); ++i) {
        Level.Rho[i] += Sub.Rho[i];
    }
    Level.VelEnd = Sub.VelEnd;

    return NoUTurn(Level.VelBeg, Level.VelEnd, Level.Rho);
}

// *******************
bool NoUTurnMC::NoUTurn(const std::vector<double>& VelBeg, const std::vector<double>& VelEnd,
                        const std::vector<double>& Rho) const {
// *******************
    double DotBeg = 0.;
    double DotEnd = 0.;
    for (size_t i = 0; i < Rho.size(); ++i) {
        DotBeg += VelBeg[i] * Rho[i];
        DotEnd += VelEnd[i] * Rho[i];
    }
    // Both ends have to still be heading outwards, which makes the test orientation free
    return DotBeg > 0. && DotEnd > 0.;
}

// *******************
void NoUTurnMC::SaveEnd(TrajectoryEnd& End) {
// *******************
    SavePosition(End.Position);
    End.Velocity = velocity;
    End.Gradient = gradient;
}

// *******************
void NoUTurnMC::RestoreEnd(const TrajectoryEnd& End) {
// *******************
    RestorePosition(End.Position);
    velocity = End.Velocity;
    gradient = End.Gradient;
}

// *******************
void NoUTurnMC::SaveCandidate(Candidate& Cand, const double Potential, const double Kinetic) {
// *******************
    SavePosition(Cand.Position);
    Cand.Potential = Potential;
    Cand.Kinetic = Kinetic;
    // The likelihood breakdown travels with the candidate so the writeback needs no reweight
    Cand.SystLlh = syst_llh;
    Cand.SampleLlh = sample_llh;
}

// *******************
void NoUTurnMC::WriteBackSelected(const double HamiltonianInit) {
// *******************
    // The chain has to record the state we evaluated, not whichever leaf the tree visited last
    RestorePosition(selected.Position);
    logLProp = selected.Potential;
    syst_llh = selected.SystLlh;
    sample_llh = selected.SampleLlh;

    accProb = (n_leapfrog_step > 0)
        ? sum_metro_prob / static_cast<double>(n_leapfrog_step) : 0.;
    energy_error_step = (selected.Potential + selected.Kinetic) - HamiltonianInit;
    hamiltonian_init_step = HamiltonianInit;
}

// **********************
double NoUTurnMC::AcceptanceProbability() {
// **********************
    return accProb;
}

// *******************
int NoUTurnMC::ThrowMatrixVersionSum() const {
// *******************
    int Sum = 0;
    for (const auto& syst : systematics) {
        Sum += syst->GetThrowMatrixUpdates();
    }
    return Sum;
}

// *******************
void NoUTurnMC::InitialiseStepSize() {
// *******************
    if (!adapt_step_size) return;

    if (find_initial_step_size) {
        const double Configured = delta_t;
        delta_t = FindReasonableStepSize();
        MACH3LOG_INFO("Starting step size {} found by doubling from the configured {}", delta_t, Configured);
    }
    dual_averaging.Restart(delta_t);
}

// *******************
double NoUTurnMC::TryStepSize(const TrajectoryEnd& Start, const double StepSize,
                              const double HamiltonianInit) {
// *******************
    RestoreEnd(Start);
    out_of_bounds = false;

    UpdateVelocity(StepSize);
    UpdatePosition(StepSize);
    // A step big enough to bounce off a wall is a step we don't want
    if (out_of_bounds) {
        out_of_bounds = false;
        return -std::numeric_limits<double>::infinity();
    }
    UpdateGradU();
    UpdateVelocity(StepSize);

    const double Hamiltonian = PotentialEnergy() + KineticEnergyAndMomentum(momentum);
    if (!std::isfinite(Hamiltonian)) return -std::numeric_limits<double>::infinity();

    return HamiltonianInit - Hamiltonian;
}

// *******************
// Double or halve until one leapfrog step crosses 0.5 acceptance, so the warmup doesn't have to
// spend itself travelling orders of magnitude
double NoUTurnMC::FindReasonableStepSize() {
// *******************
    RestoreCurrentPosition();
    ReweightSamples();
    UpdateGradU();
    const double PotentialInit = PotentialEnergy();

    ThrowVelocity();
    const double KineticInit = KineticEnergyAndMomentum(momentum);
    const double HamiltonianInit = PotentialInit + KineticInit;

    MACH3LOG_INFO("  [stepsize hunt] U0 {:.6g} KE0 {:.6g} H0 {:.6g}", PotentialInit, KineticInit, HamiltonianInit);

    TrajectoryEnd Start;
    Start.Position.assign(velocity.size(), 0.);
    SaveEnd(Start);

    double StepSize = delta_t;
    int Direction = 0;
    constexpr int MaxTries = 50;
    constexpr double LogHalf = -0.6931471805599453;

    for (int Try = 0; Try < MaxTries; ++Try) {
        const double LogRatio = TryStepSize(Start, StepSize, HamiltonianInit);
        MACH3LOG_INFO("  [stepsize hunt] try {} step {:.6g} logRatio {:.6g}", Try, StepSize, LogRatio);
        const double Want = (LogRatio > LogHalf) ? 1. : -1.;

        // Once we cross the threshold we've bracketed it, so stop
        if (Direction == 0) {
            Direction = static_cast<int>(Want);
        } else if (static_cast<int>(Want) != Direction) {
            break;
        }

        StepSize *= (Direction > 0) ? 2. : 0.5;
        if (StepSize < dual_averaging.MinStepSize || StepSize > dual_averaging.MaxStepSize) {
            MACH3LOG_WARN("Hunting for a starting step size ran off to {}, which usually means the "
                          "gradient or the mass matrix is wrong", StepSize);
            break;
        }
    }

    // Leave the sampler exactly where it started, and consistent
    RestoreEnd(Start);
    ReweightSamples();
    UpdateGradU();
    return StepSize;
}

// *******************
void NoUTurnMC::UpdateStepSize() {
// *******************
    if (!adapt_step_size) return;

    if (step < stepStart + warmup_steps) {
        if (metric_changed && restart_on_metric_update) {
            MACH3LOG_DEBUG("Mass matrix was rebuilt at step {}, restarting dual averaging from {}",
                           step, delta_t);
            dual_averaging.Restart(delta_t);
        }
        // accProb carries the mean acceptance over the trajectory, which is what we steer on
        dual_averaging.Update(accProb);
        delta_t = dual_averaging.StepSize;
    } else if (!step_size_frozen) {
        dual_averaging.Freeze();
        delta_t = dual_averaging.StepSize;
        step_size_frozen = true;
        MACH3LOG_INFO("NUTS warmup done after {} steps, step size frozen at {}", warmup_steps, delta_t);
        if (dual_averaging.NClamped > 0) {
            MACH3LOG_WARN("The step size hit a clamp {} times during warmup", dual_averaging.NClamped);
        }
    } else if (metric_changed) {
        MACH3LOG_WARN("The mass matrix was rebuilt at step {} but the step size is already frozen at "
                      "{}, so the two no longer match. Either warm up for longer than the adaption's "
                      "EndUpdate, or stop the matrix updating sooner", step, delta_t);
    }
}

// **********************
double NoUTurnMC::LogAddExp(const double A, const double B) {
// **********************
    const double Max = std::max(A, B);
    const double Min = std::min(A, B);
    return Max + std::log1p(std::exp(Min - Max));
}
