#include "Algorithms/HamiltonianMC.h"
// HW Hamiltonian MCMC!

// ***************
HamiltonianMC::HamiltonianMC(Manager *const FitManager) : MCMCBase(FitManager) {
// ***************
    AlgorithmName = "HamiltonianMC";
    n_leapfrog = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["NLeapFrog"], 10, __FILE__, __LINE__);
    //HW DeltaT is in units of prior sigma, so order one rather than the tiny identity mass values
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
            MACH3LOG_INFO("{} has adaption enabled, HMC will pick up each adapted matrix as its mass "
                          "matrix", systematics[s]->GetName());
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
    // Priors come first and overwrite, which resets the gradient for this evaluation
    for (size_t s = 0; s < systematics.size(); ++s) {
        const std::vector<double> SystGrad = systematics[s]->CalcGradient();
        const size_t Offset = static_cast<size_t>(syst_offset[s]);
        for (size_t i = 0; i < SystGrad.size(); ++i) {
            gradient[Offset + i] = SystGrad[i];
        }
    }

    // Then the samples add on top, already reweighted by UpdatePosition or ProposeStep
    for (size_t i = 0; i < samples.size(); ++i) {
        const std::vector<double> SampGrad = samples[i]->GetGradLikelihood();
        const size_t Offset = static_cast<size_t>(sample_offset[i]);
        for (size_t p = 0; p < SampGrad.size(); ++p) {
            gradient[Offset + p] += SampGrad[p];
        }
    }
}

// *******************
void HamiltonianMC::SetupMassMatrix() {
// *******************
    par_fixed.assign(velocity.size(), 0);
    mass.resize(systematics.size());

    int MaxPars = 0;
    for (size_t s = 0; s < systematics.size(); ++s) {
        mass[s].Offset = syst_offset[s];
        mass[s].NPars = systematics[s]->GetNumParams();
        MaxPars = std::max(MaxPars, mass[s].NPars);

        for (int i = 0; i < mass[s].NPars; ++i) {
            if (systematics[s]->IsParameterFixed(i)) par_fixed[static_cast<size_t>(mass[s].Offset + i)] = 1;
        }
    }
    chol_solve.assign(static_cast<size_t>(MaxPars), 0.);

    for (size_t s = 0; s < systematics.size(); ++s) {
        BuildMassBlock(s);
    }
}

// *******************
void HamiltonianMC::BuildMassBlock(const size_t s) {
// *******************
    MassBlock& Block = mass[s];

    // M^-1 is the throw matrix, which is the covariance until adaption replaces it with a better
    // estimate of the target. Either way the prior dynamics become a unit frequency oscillator
    const TMatrixDSym* Throw = systematics[s]->GetThrowMatrix();
    if (Throw == nullptr || Throw->GetNrows() != Block.NPars) {
        MACH3LOG_ERROR("{} has no usable throw matrix to build a mass matrix from", systematics[s]->GetName());
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    TMatrixDSym InvMass(*Throw);

    // A fixed parameter mustn't move, so decouple it and leave a unit mass on its diagonal to keep
    // the matrix decomposable. ZeroFixedVelocities does the rest
    for (int i = 0; i < Block.NPars; ++i) {
        if (!systematics[s]->IsParameterFixed(i)) continue;
        for (int j = 0; j < Block.NPars; ++j) {
            InvMass(i, j) = 0.;
            InvMass(j, i) = 0.;
        }
        InvMass(i, i) = 1.;
    }

    if (!M3::CanDecomposeMatrix(InvMass)) {
        MACH3LOG_WARN("Throw matrix of {} isn't positive definite, nudging the diagonal", systematics[s]->GetName());
        M3::MakeMatrixPosDef(&InvMass);
    }
    const auto Chol = M3::GetCholeskyDecomposedMatrix(InvMass, systematics[s]->GetName());

    //KS: Flat row major copies, ROOT's indexing is too slow for the per kick matrix multiply
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
    Block.BuiltAtUpdate = systematics[s]->GetThrowMatrixUpdates();

    MACH3LOG_DEBUG("Mass matrix for {} built from its {}x{} throw matrix", systematics[s]->GetName(),
                   Block.NPars, Block.NPars);
}

// *******************
void HamiltonianMC::RefreshMassMatrix() {
// *******************
    // Adaption only swaps the throw matrix every AdaptiveUpdateStep steps, so the O(N^3) rebuild is
    // amortised over all the steps in between
    for (size_t s = 0; s < systematics.size(); ++s) {
        if (mass[s].BuiltAtUpdate == systematics[s]->GetThrowMatrixUpdates()) continue;
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
            //HW Lower triangular, so row i only touches the first i + 1 throws
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
    //HW M^-1 = L L^T, so forward substituting L u = w gives 1/2 w^T M w as 1/2 |u|^2
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
double HamiltonianMC::KineticEnergyAndMomentum(std::vector<double>& Momentum) {
// *******************
    //HW M^-1 = L L^T, so M = L^-T L^-1. Forward substituting L u = w gives the kinetic energy as
    //1/2 |u|^2, and back substituting L^T p = u then gives p = M w, both without forming M
    Momentum.resize(velocity.size());
    double ke = 0.;
    for (const auto& Block : mass) {
        const int NPars = Block.NPars;
        const double* _restrict_ L = Block.Chol.data();
        const double* _restrict_ Vel = velocity.data() + Block.Offset;
        double* _restrict_ u = chol_solve.data();
        double* _restrict_ Mom = Momentum.data() + Block.Offset;

        for (int i = 0; i < NPars; ++i) {
            double Sum = Vel[i];
            for (int j = 0; j < i; ++j) {
                Sum -= L[i * NPars + j] * u[j];
            }
            u[i] = Sum / L[i * NPars + i];
            ke += 0.5 * u[i] * u[i];
        }

        //HW L^T is upper triangular, so row i only touches the momenta already solved for
        for (int i = NPars - 1; i >= 0; --i) {
            double Sum = u[i];
            for (int j = i + 1; j < NPars; ++j) {
                Sum -= L[j * NPars + i] * Mom[j];
            }
            Mom[i] = Sum / L[i * NPars + i];
        }
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
            Vel[i] -= 0.5 * Step * Sum;
        }
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
        const int NPars = systematics[s]->GetNumParams();
        const size_t Offset = static_cast<size_t>(syst_offset[s]);

        for (int i = 0; i < NPars; ++i) {
            const size_t k = Offset + static_cast<size_t>(i);
            //HW Carrying the velocity rather than the momentum keeps the drift axis aligned, so the
            //reflection below stays a simple per parameter sign flip
            double NewVal = systematics[s]->GetParProp(i) + Step * velocity[k];

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
                    MACH3LOG_WARN("Parameter {} bounced {} times in one drift, step {} is too "
                                  "large for its range [{}, {}]", systematics[s]->GetParFancyName(i),
                                  NReflections, Step, LowBound, UpBound);
                    out_of_bounds = true;
                    NewVal = systematics[s]->GetParCurr(i);
                    break;
                }
            }
            systematics[s]->SetParProp(i, NewVal);
        }
    }

    // Both the gradient and the likelihood need the prediction at the new parameter values
    ReweightSamples();
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
