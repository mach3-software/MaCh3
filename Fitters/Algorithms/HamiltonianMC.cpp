#include "Algorithms/HamiltonianMC.h"
// HW Hamiltonian MCMC!

// ***************
HamiltonianMC::HamiltonianMC(Manager *const FitManager) : MCMCBase(FitManager) {
// ***************
    AlgorithmName = "HamiltonianMC";
    n_leapfrog = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["NLeapFrog"], 10, __FILE__, __LINE__);
    //HW The mass matrix is the parameter covariance, so DeltaT is in units of prior sigma and a
    //well conditioned problem wants it of order one, not the tiny values an identity mass needs
    delta_t = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["DeltaT"], 0.1, __FILE__, __LINE__);
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
    // Lay out the flat parameter ordering, switch the samples over to tracking gradients and build
    // the mass matrix, all before the base class proposes its first step
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

    // The flat ordering is every systematic's parameters concatenated, in the order they were added
    syst_offset.resize(systematics.size());
    int NPars = 0;
    for (size_t s = 0; s < systematics.size(); ++s) {
        if (systematics[s]->IsPCA()) {
            MACH3LOG_ERROR("{} is PCA'd, which HMC doesn't support. The gradients are calculated in "
                           "the un-rotated parameter space", systematics[s]->GetName());
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        /// @todo HW: the adapted covariance would make a better mass matrix than the prior one,
        /// but it changes mid chain and rebuilding the Cholesky every step is too expensive
        if (systematics[s]->GetDoAdaption()) {
            MACH3LOG_WARN("{} has adaption enabled, HMC takes its mass matrix from the prior "
                          "covariance and won't pick up the adapted one", systematics[s]->GetName());
        }
        syst_offset[s] = NPars;
        NPars += systematics[s]->GetNumParams();
    }

    if (NPars == 0) {
        MACH3LOG_ERROR("HMC needs at least one parameter to step, did you forget AddSystObj?");
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    velocity.assign(static_cast<size_t>(NPars), 0.);
    gradient.assign(static_cast<size_t>(NPars), 0.);
    rand_throw.assign(static_cast<size_t>(NPars), 0.);

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

    // Fresh velocity from the mass matrix. We keep the kinetic energy it carries so we can compare
    // the Hamiltonian at either end of the trajectory
    ThrowVelocity();
    ke_init = KineticEnergy();

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
        UpdateVelocity();
        UpdatePosition();
        // A parameter couldn't be brought back inside its bounds, this trajectory is a write off
        if (out_of_bounds) break;
        UpdateGradU();
        UpdateVelocity();
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
// Build the block diagonal mass matrix, one block per parameter handler
void HamiltonianMC::SetupMassMatrix() {
// *******************
    par_fixed.assign(velocity.size(), 0);
    mass.resize(systematics.size());

    int MaxPars = 0;
    for (size_t s = 0; s < systematics.size(); ++s) {
        MassBlock& Block = mass[s];
        Block.Offset = syst_offset[s];
        Block.NPars = systematics[s]->GetNumParams();
        MaxPars = std::max(MaxPars, Block.NPars);

        // M^-1 is the handler's covariance. That makes the prior part of the dynamics a unit
        // frequency oscillator, since the prior gradient is Sigma^-1 * Delta and M^-1 times it is
        // just Delta, so the leapfrog no longer has to resolve the smallest prior error in the fit
        const TMatrixDSym* Cov = systematics[s]->GetCovMatrix();
        if (Cov == nullptr || Cov->GetNrows() != Block.NPars) {
            MACH3LOG_ERROR("{} has no usable covariance matrix to build a mass matrix from", systematics[s]->GetName());
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        TMatrixDSym InvMass(*Cov);

        // A fixed parameter mustn't move, so decouple it from the rest and leave a unit mass on its
        // diagonal to keep the matrix decomposable. ZeroFixedVelocities does the rest
        for (int i = 0; i < Block.NPars; ++i) {
            if (!systematics[s]->IsParameterFixed(i)) continue;
            par_fixed[static_cast<size_t>(Block.Offset + i)] = 1;
            for (int j = 0; j < Block.NPars; ++j) {
                InvMass(i, j) = 0.;
                InvMass(j, i) = 0.;
            }
            InvMass(i, i) = 1.;
        }

        if (!M3::CanDecomposeMatrix(InvMass)) {
            MACH3LOG_WARN("Covariance of {} isn't positive definite, nudging the diagonal before "
                          "decomposing it", systematics[s]->GetName());
            M3::MakeMatrixPosDef(&InvMass);
        }
        const auto Chol = M3::GetCholeskyDecomposedMatrix(InvMass, systematics[s]->GetName());

        //KS: Flat row major copies, ROOT's indexing is far too slow for the per kick matrix multiply
        const size_t N = static_cast<size_t>(Block.NPars);
        Block.InvMass.resize(N * N);
        Block.Chol.resize(N * N);
        for (int i = 0; i < Block.NPars; ++i) {
            for (int j = 0; j < Block.NPars; ++j) {
                const size_t Element = static_cast<size_t>(i) * N + static_cast<size_t>(j);
                Block.InvMass[Element] = InvMass(i, j);
                Block.Chol[Element] = Chol[static_cast<size_t>(i)][static_cast<size_t>(j)];
            }
        }
        const double MegaBytes = 2. * static_cast<double>(N * N * sizeof(double)) / (1024. * 1024.);
        MACH3LOG_INFO("Mass matrix for {} taken from its {}x{} covariance, {:.1f} MB with its Cholesky factor",
                      systematics[s]->GetName(), Block.NPars, Block.NPars, MegaBytes);
    }
    chol_solve.assign(static_cast<size_t>(MaxPars), 0.);
}

// *******************
// Draw w ~ N(0, M^-1), i.e. correlate standard normal throws with the Cholesky factor of M^-1
void HamiltonianMC::ThrowVelocity() {
// *******************
    for (auto& Throw : rand_throw) {
        Throw = M3::rand::Gaus();
    }

    for (const auto& Block : mass) {
        const int NPars = Block.NPars;
        const double* _restrict_ L = Block.Chol.data();
        const double* _restrict_ Throws = rand_throw.data() + Block.Offset;
        double* _restrict_ Vel = velocity.data() + Block.Offset;

        #ifdef MULTITHREAD
        #pragma omp parallel for
        #endif
        for (int i = 0; i < NPars; ++i) {
            double Sum = 0.;
            //HW L is lower triangular, so row i only touches the first i + 1 throws
            for (int j = 0; j <= i; ++j) {
                Sum += L[i * NPars + j] * Throws[j];
            }
            Vel[i] = Sum;
        }
    }
    ZeroFixedVelocities();
}

// *******************
double HamiltonianMC::KineticEnergy() {
// *******************
    //HW KE = 1/2 w^T M w, and M^-1 = L L^T, so forward substituting L u = w gives KE = 1/2 |u|^2
    //without ever forming M itself
    double ke = 0.;
    for (const auto& Block : mass) {
        const int NPars = Block.NPars;
        const double* _restrict_ L = Block.Chol.data();
        const double* _restrict_ Vel = velocity.data() + Block.Offset;
        double* _restrict_ u = chol_solve.data();

        for (int i = 0; i < NPars; ++i) {
            double Sum = Vel[i];
            for (int j = 0; j < i; ++j) {
                Sum -= L[i * NPars + j] * u[j];
            }
            u[i] = Sum / L[i * NPars + i];
            ke += 0.5 * u[i] * u[i];
        }
    }
    return ke;
}

// *******************
void HamiltonianMC::UpdateVelocity() {
// *******************
    //HW In velocity coordinates the half kick p -> p - (dt/2) grad U becomes
    //w -> w - (dt/2) M^-1 grad U, so the metric enters here rather than in the drift
    for (const auto& Block : mass) {
        const int NPars = Block.NPars;
        const double* _restrict_ InvMass = Block.InvMass.data();
        const double* _restrict_ Grad = gradient.data() + Block.Offset;
        double* _restrict_ Vel = velocity.data() + Block.Offset;

        #ifdef MULTITHREAD
        #pragma omp parallel for
        #endif
        for (int i = 0; i < NPars; ++i) {
            double Sum = 0.;
            #ifdef MULTITHREAD
            #pragma omp simd
            #endif
            for (int j = 0; j < NPars; ++j) {
                Sum += InvMass[i * NPars + j] * Grad[j];
            }
            Vel[i] -= 0.5 * delta_t * Sum;
        }
    }
    ZeroFixedVelocities();
}

// *******************
void HamiltonianMC::ZeroFixedVelocities() {
// *******************
    //HW A fixed parameter is decoupled in the mass matrix but its unit diagonal still hands it a
    //velocity on every kick, so pin it to zero and it can never drift
    for (size_t i = 0; i < velocity.size(); ++i) {
        if (par_fixed[i]) velocity[i] = 0.;
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
            //HW The one place the mass matrix doesn't appear, which is why we carry the velocity and
            //not the momentum: the drift stays axis aligned so the reflection is a simple sign flip
            double NewVal = systematics[s]->GetParProp(i) + delta_t * velocity[k];

            const double LowBound = systematics[s]->GetLowerBound(i);
            const double UpBound = systematics[s]->GetUpperBound(i);

            // Reflect off the "wall" like a billiard: whatever is left of the drift past the bound
            int NReflections = 0;
            while (NewVal > UpBound || NewVal < LowBound) {
                NewVal = (NewVal > UpBound) ? 2. * UpBound - NewVal : 2. * LowBound - NewVal;
                velocity[k] = -velocity[k];

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
    // The velocity at the end of the trajectory, which carries the proposed kinetic energy
    const double ke_prop = KineticEnergy();

    // U is already the negative log likelihood, so H = -logL + KE at either end
    const double hamiltonian_prop = logLProp + ke_prop;
    const double hamiltonian_init = logL_init + ke_init;

    // Perfect integration would conserve H and accept everything, so this only corrects for the
    // error the finite leapfrog step size introduces
    return std::min(1., std::exp(hamiltonian_init - hamiltonian_prop));
}
