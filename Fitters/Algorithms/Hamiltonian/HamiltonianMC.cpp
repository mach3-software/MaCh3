#include "Algorithms/Hamiltonian/HamiltonianMC.h"
#include <cmath>
// HW Hamiltonian MCMC!

// ***************
HamiltonianMC::HamiltonianMC(Manager *const FitManager) : MCMCBase(FitManager) {
// ***************
    AlgorithmName = "HamiltonianMC";
    n_leapfrog = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["NLeapFrog"], 10, __FILE__, __LINE__);
    //HW DeltaT is in units of prior sigma, so order one rather than the tiny identity mass values
    delta_t = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["DeltaT"], 0.1, __FILE__, __LINE__);
    //HW Adaption's covariance is a good random walk proposal and a poor metric, so don't take it
    //unless asked. See the adapt_metric docs
    adapt_metric = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["AdaptMetric"], false, __FILE__, __LINE__);
    //HW Damping the correlations by however much adaption hasn't resolved them is what makes an
    //adapted matrix usable as a metric, see the metric_shrinkage docs
    metric_shrinkage = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["MetricShrinkage"], 0., __FILE__, __LINE__);
    metric_from_prior = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["MetricFromPrior"], false, __FILE__, __LINE__);
    if (metric_from_prior && adapt_metric) {
        MACH3LOG_ERROR("MetricFromPrior and AdaptMetric both set: the metric can't stay at the prior and follow adaption");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (metric_from_prior) MACH3LOG_INFO("Metric built from each handler's prior covariance, whatever matrix is loaded or adapted");
    ke_init = 0.;
    logL_init = 0.;

    if (metric_shrinkage < 0. || metric_shrinkage > 1.) {
        MACH3LOG_ERROR("MetricShrinkage has to be in [0, 1], you asked for {}", metric_shrinkage);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (metric_shrinkage > 0.) {
        MACH3LOG_INFO("Adapted metrics shrunk {:.0f}% of the way towards their diagonal; the prior covariance is used as it is", 100. * metric_shrinkage);
        if (!adapt_metric) {
            MACH3LOG_WARN("MetricShrinkage only does something once the metric has correlations "
                          "worth damping. With AdaptMetric false it just rescales the starting "
                          "throw matrix, which is usually already diagonal");
        }
    }

    if (n_leapfrog < 1) {
        MACH3LOG_ERROR("NLeapFrog has to be at least 1, you asked for {}", n_leapfrog);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (delta_t <= 0.) {
        MACH3LOG_ERROR("DeltaT has to be positive, you asked for {}", delta_t);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    MACH3LOG_INFO("Running HMC with {} leapfrog iterations of step size {}", n_leapfrog, delta_t);
}

// *******************
void HamiltonianMC::RunMCMC() {
// *******************
    // All of this has to happen before the base class proposes its first step
    SetupGradients();
    SetupMassMatrix();

    MCMCBase::RunMCMC();
}

// *******************
void HamiltonianMC::SetupGradients() {
// *******************
    if (CheckNodeExists(fitMan->raw()["General"]["MCMC"]["Multicanonical"]) &&
        GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["Multicanonical"]["Enabled"], false, __FILE__, __LINE__)) {
        MACH3LOG_ERROR("HMC can't be run with the multicanonical method, the umbrella penalty isn't "
                       "included in the gradient");
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    const int NPars = BuildFlatOrdering();
    velocity.assign(static_cast<size_t>(NPars), 0.);
    gradient.assign(static_cast<size_t>(NPars), 0.);
    rand_throw.assign(static_cast<size_t>(NPars), 0.);
    MatchSamplesToHandlers();

    MACH3LOG_INFO("HMC will step {} parameters across {} parameter handlers", NPars, systematics.size());
}

// *******************
int HamiltonianMC::BuildFlatOrdering() {
// *******************
    // The flat ordering is every systematic's parameters concatenated, in the order they were added
    syst_offset.resize(systematics.size());
    int NPars = 0;
    for (size_t s = 0; s < systematics.size(); ++s) {
        if (systematics[s]->IsPCA()) {
            MACH3LOG_ERROR("{} is PCA'd, which HMC doesn't support. The gradients are calculated in "
                           "the un-rotated parameter space", systematics[s]->GetName());
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        if (systematics[s]->GetDoAdaption()) {
            if (adapt_metric) {
                MACH3LOG_INFO("{} has adaption enabled, HMC will pick up each adapted matrix as its "
                              "mass matrix", systematics[s]->GetName());
                MACH3LOG_WARN("An adapted covariance is usually too ill conditioned to integrate "
                              "against. If the chain stops exploring once the adaption starts "
                              "throwing, set AdaptMetric to false");
            } else {
                MACH3LOG_INFO("{} has adaption enabled, but AdaptMetric is false, so the mass matrix "
                              "stays at its starting throw matrix. The adaption still accumulates "
                              "and saves its covariance", systematics[s]->GetName());
            }
        }
        syst_offset[s] = NPars;
        NPars += systematics[s]->GetNumParams();
    }

    if (NPars == 0) {
        MACH3LOG_ERROR("HMC needs at least one parameter to step, did you forget AddSystObj?");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    return NPars;
}

// *******************
void HamiltonianMC::MatchSamplesToHandlers() {
// *******************
    // Each sample hands its gradient back indexed by its own parameter handler, so find where that
    // handler's block sits in our flat ordering
    sample_offset.assign(samples.size(), M3::_BAD_INT_);
    for (size_t i = 0; i < samples.size(); ++i) {
        samples[i]->SetComputeGradients(true);

        const ParameterHandlerBase* SampleParHandler = samples[i]->GetGradParHandler();
        for (size_t s = 0; s < systematics.size(); ++s) {
            if (systematics[s] == SampleParHandler) {
                sample_offset[i] = syst_offset[s];
                break;
            }
        }
        if (sample_offset[i] == M3::_BAD_INT_) {
            MACH3LOG_ERROR("{} is driven by a parameter handler which wasn't added to the fit, so "
                           "its gradient can't be matched up with the parameters we step",
                           samples[i]->GetName());
            throw MaCh3Exception(__FILE__, __LINE__);
        }
    }
}

// *******************
void HamiltonianMC::DoStep() {
// *******************
    ProposeStep();
    // Does the MCMC accept this trajectory?
    accProb = AcceptanceProbability();

    if (!out_of_bounds && IsStepAccepted(accProb))
    {
        AcceptStep();
    }
}

// *******************
// Integrate a leapfrog trajectory and evaluate the likelihood where it ends up
void HamiltonianMC::ProposeStep() {
// *******************
    out_of_bounds = false;
    RefreshMassMatrix();

    // Start from the last accepted state. A rejected step leaves the proposed values at the end of
    // the trajectory we threw away, so this has to be done explicitly
    RestoreCurrentPosition();

    // Keep the kinetic energy so we can compare the Hamiltonian at either end of the trajectory
    ThrowVelocity();
    ke_init = KineticEnergy();

    // The first half kick needs dU/dx where the trajectory starts
    ReweightSamples();
    UpdateGradU();

    // Taking U here keeps the Hamiltonian self consistent instead of leaning on logLCurr, which
    // the base class only seeds after this function has already run once
    logL_init = PotentialEnergy();

    for (int f = 0; f < n_leapfrog; ++f) {
        // Half kick, drift, half kick. Adjacent half kicks compose into the usual full kick
        UpdateVelocity(delta_t);
        UpdatePosition(delta_t);
        // A parameter couldn't be brought back inside its bounds, this trajectory is a write off
        if (out_of_bounds) break;
        UpdateGradU();
        UpdateVelocity(delta_t);
    }

    // Finally the potential energy at the end of the trajectory, using the weights UpdatePosition left
    logLProp = PotentialEnergy();
    if (logLProp >= M3::_LARGE_LOGL_) out_of_bounds = true;
}

// *******************
// Refresh the gradient of the potential, U = -logL
void HamiltonianMC::UpdateGradU() {
// *******************
    // Priors come first and overwrite, which resets the gradient for this evaluation. Both loops
    // fill one reused buffer: this runs every leapfrog step, so a fresh vector per handler and per
    // sample was the bulk of the allocation on the hot path
    for (size_t s = 0; s < systematics.size(); ++s) {
        systematics[s]->CalcGradient(grad_block);
        const size_t Offset = static_cast<size_t>(syst_offset[s]);
        for (size_t i = 0; i < grad_block.size(); ++i) {
            gradient[Offset + i] = grad_block[i];
        }
    }

    // Then the samples add on top, already reweighted by UpdatePosition or ProposeStep
    for (size_t i = 0; i < samples.size(); ++i) {
        samples[i]->GetGradLikelihood(grad_block);
        const size_t Offset = static_cast<size_t>(sample_offset[i]);
        for (size_t p = 0; p < grad_block.size(); ++p) {
            gradient[Offset + p] += grad_block[p];
        }
    }

    //HW A frozen parameter has no business contributing a force. Its row of the metric is already
    //decoupled so this changes no dynamics, but it keeps the gradient honest for anything which
    //reads it, and it means a parameter whose likelihood derivative is dropped can't be pushed
    //around by its prior alone
    for (size_t i = 0; i < gradient.size() && i < par_fixed.size(); ++i) {
        if (par_fixed[i]) gradient[i] = 0.;
    }
}

// *******************
void HamiltonianMC::SetupMassMatrix() {
// *******************
    //HW Gives a subclass the chance to say which parameters it is sampling with another kernel,
    //before any of that feeds into par_fixed and the metric below
    SetupExternalBlock();

    par_fixed.assign(velocity.size(), 0);
    mass.assign(systematics.size(), MetricBlock());
    mass_built_at.assign(systematics.size(), -1);
    throw_updates_at_start.resize(systematics.size());
    for (size_t s = 0; s < systematics.size(); ++s) throw_updates_at_start[s] = systematics[s]->GetThrowMatrixUpdates();

    int MaxPars = 0;
    for (size_t s = 0; s < systematics.size(); ++s) {
        const int NPars = systematics[s]->GetNumParams();
        MaxPars = std::max(MaxPars, NPars);
        for (int i = 0; i < NPars; ++i) {
            const size_t Flat = static_cast<size_t>(syst_offset[s] + i);
            //HW A parameter another kernel owns is frozen for us in exactly the same way a fixed
            //one is, so the two go through the same path from here on
            const bool External = !par_external.empty() && par_external[Flat] != 0;
            if (systematics[s]->IsParameterFixed(i) || External) par_fixed[Flat] = 1;
        }
    }
    chol_solve.assign(static_cast<size_t>(MaxPars), 0.);
    SetupDriftBounds(MaxPars);

    for (size_t s = 0; s < systematics.size(); ++s) {
        BuildMassBlock(s);
    }
}

// *******************
void HamiltonianMC::BuildMassBlock(const size_t s) {
// *******************
    // M^-1 is the throw matrix, which is the covariance until adaption replaces it with a better
    // estimate of the target. Either way the prior dynamics become a unit frequency oscillator
    const int NPars = systematics[s]->GetNumParams();
    const TMatrixDSym* Throw = metric_from_prior ? systematics[s]->GetInitialThrowMatrix() : systematics[s]->GetThrowMatrix();
    if (Throw == nullptr || Throw->GetNrows() != NPars) {
        MACH3LOG_ERROR("{} has no usable throw matrix to build a mass matrix from", systematics[s]->GetName());
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    const auto First = par_fixed.begin() + syst_offset[s];
    const std::vector<char> Frozen(First, First + NPars);
    std::vector<double> PriorWidth(static_cast<size_t>(NPars));
    for (int i = 0; i < NPars; ++i) PriorWidth[static_cast<size_t>(i)] = systematics[s]->GetDiagonalError(i);

    //HW Shrinkage is for an adapted covariance's noisy correlations, not the prior's exact ones
    const bool Adapted = !metric_from_prior && systematics[s]->GetThrowMatrixUpdates() > throw_updates_at_start[s];
    mass[s].Build(*Throw, Frozen, Adapted ? metric_shrinkage : 0., PriorWidth, systematics[s]->GetName());
    mass_built_at[s] = systematics[s]->GetThrowMatrixUpdates();
    MACH3LOG_DEBUG("{} mass matrix for {} built from its {}x{} throw matrix", mass[s].IsDiagonal() ? "Diagonal" : "Dense",
                   systematics[s]->GetName(), NPars, NPars);
}

// *******************
void HamiltonianMC::RefreshMassMatrix() {
// *******************
    // Leave the metric exactly as SetupMassMatrix() built it, whatever adaption does to the throw
    // matrix underneath us
    if (!adapt_metric || metric_from_prior) return;

    // Adaption only swaps the throw matrix every AdaptiveUpdateStep steps, so the O(N^3) rebuild is
    // amortised over all the steps in between
    for (size_t s = 0; s < systematics.size(); ++s) {
        if (mass_built_at[s] == systematics[s]->GetThrowMatrixUpdates()) continue;
        MACH3LOG_INFO("Throw matrix of {} was re-adapted, rebuilding its mass matrix", systematics[s]->GetName());
        BuildMassBlock(s);
    }
}

// *******************
// Draw w ~ N(0, M^-1), i.e. correlate standard normal throws with the Cholesky factor of M^-1
void HamiltonianMC::ThrowVelocity() {
// *******************
    for (auto& Throw : rand_throw) {
        Throw = M3::rand::Gaus();
    }
    for (size_t s = 0; s < mass.size(); ++s) {
        const size_t Offset = static_cast<size_t>(syst_offset[s]);
        mass[s].ThrowVelocity(rand_throw.data() + Offset, velocity.data() + Offset);
    }
    ZeroFixedVelocities();
}

// *******************
double HamiltonianMC::KineticEnergy() {
// *******************
    double ke = 0.;
    for (size_t s = 0; s < mass.size(); ++s) {
        ke = mass[s].AddKineticEnergy(velocity.data() + syst_offset[s], chol_solve.data(), ke);
    }
    return ke;
}

// *******************
double HamiltonianMC::KineticEnergyAndMomentum(std::vector<double>& Momentum) {
// *******************
    Momentum.resize(velocity.size());
    double ke = 0.;
    for (size_t s = 0; s < mass.size(); ++s) {
        const size_t Offset = static_cast<size_t>(syst_offset[s]);
        ke = mass[s].AddKineticEnergyAndMomentum(velocity.data() + Offset, chol_solve.data(), Momentum.data() + Offset, ke);
    }
    return ke;
}

// *******************
void HamiltonianMC::RestoreCurrentPosition() {
// *******************
    for (size_t s = 0; s < systematics.size(); ++s) {
        const int NPars = systematics[s]->GetNumParams();
        for (int i = 0; i < NPars; ++i) {
            systematics[s]->SetParProp(i, systematics[s]->GetParCurr(i));
        }
    }
}

// *******************
void HamiltonianMC::ReweightSamples() {
// *******************
    for (auto& samp : samples) {
        samp->Reweight();
    }
}

// *******************
double HamiltonianMC::PotentialEnergy() {
// *******************
    double llh = 0.0;
    for (size_t s = 0; s < systematics.size(); ++s) {
        syst_llh[s] = systematics[s]->GetLikelihood();
        llh += syst_llh[s];
    }

    // A reflected parameter stays in bounds, so only a genuinely mad prior penalty gets us here,
    // and then the samples aren't worth evaluating. The caller decides what that means for the step
    if (llh >= M3::_LARGE_LOGL_) {
        for (size_t i = 0; i < samples.size(); ++i) {
            sample_llh[i] = M3::_LARGE_LOGL_;
        }
        return llh;
    }

    for (size_t i = 0; i < samples.size(); ++i) {
        sample_llh[i] = samples[i]->GetLikelihood();
        llh += sample_llh[i];
    }
    return llh;
}

// *******************
void HamiltonianMC::SavePosition(std::vector<double>& Out) const {
// *******************
    Out.resize(velocity.size());
    for (size_t s = 0; s < systematics.size(); ++s) {
        const int NPars = systematics[s]->GetNumParams();
        const size_t Offset = static_cast<size_t>(syst_offset[s]);
        for (int i = 0; i < NPars; ++i) {
            Out[Offset + static_cast<size_t>(i)] = systematics[s]->GetParProp(i);
        }
    }
}

// *******************
void HamiltonianMC::RestorePosition(const std::vector<double>& In) {
// *******************
    //HW In place, one parameter at a time. The samples and splines hold raw pointers into the
    //proposed values, so the vector itself must never be swapped or resized under them
    for (size_t s = 0; s < systematics.size(); ++s) {
        const int NPars = systematics[s]->GetNumParams();
        const size_t Offset = static_cast<size_t>(syst_offset[s]);
        for (int i = 0; i < NPars; ++i) {
            systematics[s]->SetParProp(i, In[Offset + static_cast<size_t>(i)]);
        }
    }
}

// *******************
void HamiltonianMC::UpdateVelocity(const double Step) {
// *******************
    //HW Half a kick, which in velocity coordinates is w -> w - (h/2) M^-1 grad U
    for (size_t s = 0; s < mass.size(); ++s) {
        const size_t Offset = static_cast<size_t>(syst_offset[s]);
        mass[s].Kick(velocity.data() + Offset, gradient.data() + Offset, Step);
    }
    ZeroFixedVelocities();
}

// *******************
void HamiltonianMC::ZeroFixedVelocities() {
// *******************
    //HW The unit diagonal still hands a fixed parameter a velocity on every kick, so pin it to zero
    for (size_t i = 0; i < velocity.size(); ++i) {
        if (par_fixed[i]) velocity[i] = 0.;
    }
}

// *******************
void HamiltonianMC::UpdatePosition(const double Step) {
// *******************
    for (size_t s = 0; s < systematics.size(); ++s) {
        if (!DriftBlock(s, Step)) {
            out_of_bounds = true;
            break;
        }
    }

    // Both the gradient and the likelihood need the prediction at the new parameter values
    ReweightSamples();
}

// *******************
void HamiltonianMC::SetupDriftBounds(const int MaxPars) {
// *******************
    drift_bounds.assign(systematics.size(), DriftBounds());
    for (size_t s = 0; s < systematics.size(); ++s) {
        const ParameterHandlerBase& Syst = *systematics[s];
        const size_t N = static_cast<size_t>(Syst.GetNumParams());
        DriftBounds& B = drift_bounds[s];
        B.Low.resize(N);
        B.Up.resize(N);
        B.Circular.assign(N, 0);
        for (size_t u = 0; u < N; ++u) {
            const int i = static_cast<int>(u);
            double CircLow = 0., CircUp = 0.;
            //HW A circular parameter (delta_cp) wraps round instead: the drift on a circle is just as
            //volume preserving and reversible, and NUTS's U-turn test only looks at momenta
            if (Syst.GetCircularBounds(i, CircLow, CircUp) && CircUp > CircLow) {
                B.Circular[u] = 1; B.Low[u] = CircLow; B.Up[u] = CircUp;
            } else {
                B.Low[u] = Syst.GetLowerBound(i); B.Up[u] = Syst.GetUpperBound(i);
            }
        }
    }
    // One buffer wide enough for any handler, so a drift allocates nothing
    drift_x.assign(static_cast<size_t>(MaxPars), 0.);
    drift_move.assign(static_cast<size_t>(MaxPars), 0.);
}

// *******************
int HamiltonianMC::FirstWall(const size_t s, const double Remaining, double& HitAt) const {
// *******************
    const DriftBounds& B = drift_bounds[s];
    // The scratch is sized to the widest handler, so stop at this one's parameter count
    const size_t N = B.Low.size();
    int Hit = -1;
    HitAt = Remaining;
    for (size_t u = 0; u < N; ++u) {
        if (B.Circular[u]) continue;
        const double Move = drift_move[u], X = drift_x[u];
        double t = Remaining + 1.;
        if (Move > 0. && X + Remaining * Move > B.Up[u]) t = (B.Up[u] - X) / Move;
        else if (Move < 0. && X + Remaining * Move < B.Low[u]) t = (B.Low[u] - X) / Move;
        if (t < HitAt) { HitAt = std::max(t, 0.); Hit = static_cast<int>(u); }
    }
    return Hit;
}

// *******************
void HamiltonianMC::WrapCircular(const size_t s) {
// *******************
    const DriftBounds& B = drift_bounds[s];
    const size_t N = B.Low.size();
    for (size_t u = 0; u < N; ++u) {
        if (!B.Circular[u] || (drift_x[u] <= B.Up[u] && drift_x[u] >= B.Low[u])) continue;
        const double Width = B.Up[u] - B.Low[u];
        drift_x[u] = B.Low[u] + (drift_x[u] - B.Low[u]) - Width * std::floor((drift_x[u] - B.Low[u]) / Width);
    }
}

// *******************
bool HamiltonianMC::DriftBlock(const size_t s, const double Step) {
// *******************
    //HW A billiard in the metric. Every parameter drifts by Step times its velocity; when one meets a
    //bound, everything moves up to that moment, the velocity is reflected in the metric (which, with
    //a dense metric, also changes the velocities of the parameters correlated with it, see
    //MetricBlock::Reflect), and the rest of the drift carries on from there. Reflecting each
    //parameter on its own after a full drift is only right for a diagonal metric: with a dense one
    //it changes the kinetic energy at every bounce, which shows up as divergences
    ParameterHandlerBase* Syst = systematics[s];
    const int NPars = Syst->GetNumParams();
    double* Velocity = velocity.data() + static_cast<size_t>(syst_offset[s]);
    // The walls are already in drift_bounds[s], so only the values have to be read each drift
    for (int i = 0; i < NPars; ++i) drift_x[static_cast<size_t>(i)] = Syst->GetParProp(i);

    double Remaining = 1.;   // fraction of the drift still to go
    for (int NReflections = 0; ; ++NReflections) {
        // Drift everything up to the first wall (or to the end of the drift if there is none)
        for (int i = 0; i < NPars; ++i) drift_move[static_cast<size_t>(i)] = Step * Velocity[i];
        double HitAt = Remaining;
        const int Hit = FirstWall(s, Remaining, HitAt);
        for (int i = 0; i < NPars; ++i) drift_x[static_cast<size_t>(i)] += HitAt * Step * Velocity[i];
        Remaining -= HitAt;
        if (Hit < 0) break;

        // Put the parameter exactly on the wall it hit, and bounce
        const size_t h = static_cast<size_t>(Hit);
        drift_x[h] = (Step * Velocity[Hit] > 0.) ? drift_bounds[s].Up[h] : drift_bounds[s].Low[h];
        mass[s].Reflect(Velocity, Hit);
        if (NReflections + 1 >= MaxReflections) {
            // DeltaT is far too big for this block's ranges, give up on the trajectory rather than spin
            MACH3LOG_WARN("{} bounced {} times in one drift (last off {}), step {} is too large for its ranges",
                          Syst->GetName(), MaxReflections, Syst->GetParFancyName(Hit), Step);
            return false;
        }
    }

    WrapCircular(s);
    for (int i = 0; i < NPars; ++i) Syst->SetParProp(i, drift_x[static_cast<size_t>(i)]);
    return true;
}

// **********************
// Do we accept the end of the trajectory?
double HamiltonianMC::AcceptanceProbability() {
// **********************
    // The velocity at the end of the trajectory, which carries the proposed kinetic energy
    const double ke_prop = KineticEnergy();

    // U is already the negative log likelihood, so H = -logL + KE at either end
    const double hamiltonian_prop = logLProp + ke_prop;
    const double hamiltonian_init = logL_init + ke_init;

    // Perfect integration would conserve H and accept everything, so this only corrects for the
    // error the finite leapfrog step size introduces
    return std::min(1., std::exp(hamiltonian_init - hamiltonian_prop));
}
