#include "Algorithms/HamiltonianMC.h"
// HW Hamiltonian MCMC!

// ***************
HamiltonianMC::HamiltonianMC(Manager *const FitManager) : MCMCBase(FitManager) {
// ***************
    AlgorithmName = "HamiltonianMC";
    n_leapfrog = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["NLeapFrog"], 10, __FILE__, __LINE__);
    delta_t = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["DeltaT"], 0.05, __FILE__, __LINE__);
    ke_init = 0.;
    logL_init = 0.;

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
    // Lay out the flat parameter ordering and switch the samples over to tracking gradients, both
    // have to happen before the base class proposes its first step
    SetupGradients();

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
            MACH3LOG_WARN("{} has adaption enabled, HMC ignores the adapted step sizes and uses "
                          "DeltaT for every parameter", systematics[s]->GetName());
        }
        syst_offset[s] = NPars;
        NPars += systematics[s]->GetNumParams();
    }

    if (NPars == 0) {
        MACH3LOG_ERROR("HMC needs at least one parameter to step, did you forget AddSystObj?");
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    momentum.assign(static_cast<size_t>(NPars), 0.);
    gradient.assign(static_cast<size_t>(NPars), 0.);

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

    MACH3LOG_INFO("HMC will step {} parameters across {} parameter handlers", NPars, systematics.size());
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

    // Start from the last accepted state. A rejected step leaves the proposed values at the end of
    // the trajectory we threw away, so this has to be done explicitly
    for (size_t s = 0; s < systematics.size(); ++s) {
        const int NPars = systematics[s]->GetNumParams();
        for (int i = 0; i < NPars; ++i) {
            systematics[s]->SetParProp(i, systematics[s]->GetParCurr(i));
        }
    }

    // Fresh momentum, assuming the mass matrix is the identity. Here the KE is just 1/2 p^T p and we
    // keep it so we can compare the Hamiltonian at either end of the trajectory
    ke_init = 0.;
    for (auto& mom : momentum) {
        mom = M3::rand::Gaus();
        ke_init += 0.5 * mom * mom;
    }

    // The first half kick needs dU/dx where the trajectory starts
    for (auto& samp : samples) {
        samp->Reweight();
    }
    UpdateGradU();

    // U where the trajectory starts. We've just reweighted for it anyway, and taking it here keeps
    // the Hamiltonian self consistent instead of leaning on logLCurr, which the base class only
    // seeds after this very function has run for the first time
    logL_init = 0.;
    for (const auto& syst : systematics) {
        logL_init += syst->GetLikelihood();
    }
    for (const auto& samp : samples) {
        logL_init += samp->GetLikelihood();
    }

    for (int f = 0; f < n_leapfrog; ++f) {
        // Half kick, drift, half kick. Back to back half kicks across iterations add up to the full
        // kick of the usual leapfrog, and UpdatePosition leaves the samples reweighted for UpdateGradU
        UpdateMomentum();
        UpdatePosition();
        // A parameter couldn't be brought back inside its bounds, this trajectory is a write off
        if (out_of_bounds) break;
        UpdateGradU();
        UpdateMomentum();
    }

    // Finally the potential energy at the end of the trajectory, using the weights UpdatePosition left
    double llh = 0.0;
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        syst_llh[s] = systematics[s]->GetLikelihood();
        llh += syst_llh[s];
    }

    // A reflected parameter stays in bounds, so only a genuinely mad prior penalty gets us here
    if (llh >= M3::_LARGE_LOGL_)
    {
        out_of_bounds = true;
        for (size_t i = 0; i < samples.size(); ++i) {
            sample_llh[i] = M3::_LARGE_LOGL_;
        }
    }
    else
    {
        for (size_t i = 0; i < samples.size(); ++i)
        {
            sample_llh[i] = samples[i]->GetLikelihood();
            llh += sample_llh[i];
        }
    }

    logLProp = llh;
}

// *******************
// Refresh the gradient of the potential, U = -logL
void HamiltonianMC::UpdateGradU() {
// *******************
    // The prior penalties come first and overwrite, which resets the gradient for this evaluation
    for (size_t s = 0; s < systematics.size(); ++s) {
        const std::vector<double> SystGrad = systematics[s]->CalcGradient();
        const size_t Offset = static_cast<size_t>(syst_offset[s]);
        for (size_t i = 0; i < SystGrad.size(); ++i) {
            gradient[Offset + i] = SystGrad[i];
        }
    }

    // Then the samples add on top. They were reweighted by UpdatePosition, or by ProposeStep for the
    // first evaluation of a trajectory
    for (size_t i = 0; i < samples.size(); ++i) {
        const std::vector<double> SampGrad = samples[i]->GetGradLikelihood();
        const size_t Offset = static_cast<size_t>(sample_offset[i]);
        for (size_t p = 0; p < SampGrad.size(); ++p) {
            gradient[Offset + p] += SampGrad[p];
        }
    }
}

// *******************
void HamiltonianMC::UpdateMomentum() {
// *******************
    for (size_t i = 0; i < momentum.size(); ++i) {
        momentum[i] -= 0.5 * delta_t * gradient[i];
    }
}

// *******************
void HamiltonianMC::UpdatePosition() {
// *******************
    for (size_t s = 0; s < systematics.size(); ++s) {
        const int NPars = systematics[s]->GetNumParams();
        const size_t Offset = static_cast<size_t>(syst_offset[s]);

        for (int i = 0; i < NPars; ++i) {
            const size_t k = Offset + static_cast<size_t>(i);
            double NewVal = systematics[s]->GetParProp(i) + delta_t * momentum[k];

            const double LowBound = systematics[s]->GetLowerBound(i);
            const double UpBound = systematics[s]->GetUpperBound(i);

            // Reflect off the "wall" like a billiard: whatever is left of the drift past the bound
            int NReflections = 0;
            while (NewVal > UpBound || NewVal < LowBound) {
                NewVal = (NewVal > UpBound) ? 2. * UpBound - NewVal : 2. * LowBound - NewVal;
                momentum[k] = -momentum[k];

                if (++NReflections >= MaxReflections) {
                    // DeltaT is far too big for this parameter's range, give up on the trajectory
                    // rather than spin here
                    MACH3LOG_WARN("Parameter {} bounced {} times in one drift, DeltaT = {} is too "
                                  "large for its range [{}, {}]", systematics[s]->GetParFancyName(i),
                                  NReflections, delta_t, LowBound, UpBound);
                    out_of_bounds = true;
                    NewVal = systematics[s]->GetParCurr(i);
                    break;
                }
            }
            systematics[s]->SetParProp(i, NewVal);
        }
    }

    // Both the gradient and the likelihood need the prediction at the new parameter values
    for (auto& samp : samples) {
        samp->Reweight();
    }
}

// **********************
// Do we accept the end of the trajectory?
double HamiltonianMC::AcceptanceProbability() {
// **********************
    // The momentum at the end of the trajectory, which carries the proposed kinetic energy
    double ke_prop = 0.;
    for (const auto& mom : momentum) {
        ke_prop += 0.5 * mom * mom;
    }

    // U is already the negative log likelihood, so H = -logL + KE at either end
    const double hamiltonian_prop = logLProp + ke_prop;
    const double hamiltonian_init = logL_init + ke_init;

    // Perfect integration would conserve H and accept everything, so this only corrects for the
    // error the finite leapfrog step size introduces
    return std::min(1., std::exp(hamiltonian_init - hamiltonian_prop));
}
