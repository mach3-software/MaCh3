#include "Algorithms/HamiltonianMC.h"
// HW Hamiltonian MCMC!

// ***************
HamiltonianMC::HamiltonianMC(Manager *const FitManager) : MCMCBase(FitManager) {
// ***************
    AlgorithmName = "HamiltonianMC";
    n_leapfrog = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["NLeapFrog"], 10, __FILE__, __LINE__);
    delta_t = GetFromManager<float>(fitMan->raw()["General"]["MCMC"]["DeltaT"], 0.05, __FILE__, __LINE__);
}

void HamiltonianMC::RunMCMC(){
    // Need to setup initial momentum vectors
    int n_systs = 0;
    for(const auto& syst: systematics){
        n_systs+=syst->GetNumParams();
    }

    momentum = std::vector<double>(n_systs, 0.0);
    gradient = std::vector<double>(n_systs, 0.0);

    MCMCBase::RunMCMC();
}

void HamiltonianMC::DoStep(){
    // We now set the initial momentum vector (assuming the mass matrix is the identity...)
    ke_init = 0 ;

    for(size_t i=0; i<momentum.size(); ++i){
        double mom = M3::rand::Gaus();
        // We store the initial momentum so we can calculate the hamiltonian at the end
        // Here the ke is just 1/2 p^T p
        ke_init += 0.5*std::pow(mom, 2);
        momentum[i] = mom;
    }


    for(int f=0; f<n_leapfrog; f++){
        ProposeStep();
    }
}

void HamiltonianMC::UpdateGradU(){
    int offset = 0;
    for(auto& syst : systematics){
        auto syst_grad = syst->CalcGradient();
        int npars = syst->GetNumParams();

        for(int i=0; i<npars; ++i){
            gradient[offset+i] = syst_grad[i];
        }
        offset+=npars;
    }

    // Reweight the sample handlers
    for(auto& samp : samples){
        samp->Reweight();
        auto samp_grad = samp->GetGradLikelihood();
        for(size_t g=0; g<gradient.size(); ++g){
            gradient[g] = samp_grad[g];
        }
    }
}

void HamiltonianMC::UpdateMomentum(){
    for(size_t i=0; i<momentum.size(); ++i){
        momentum[i] -= (delta_t/2) * gradient[i];
    }
}


void HamiltonianMC::ProposeStep(){
    // Leap frogging!
    // First we update the momentum
    UpdateMomentum();
    // Now we update our steps
    int offset = 0;
    // First we set the systematic values
    for(auto& syst: systematics){
        for(int i = 0; i < syst->GetNumParams(); ++i){
            const size_t k = offset + i;
            const double curr_val = syst->GetParProp(i);

            syst->SetParProp(i, curr_val + delta_t * momentum[k]);

            // If it's oob we reflect off the "wall" but stay put
            if(!syst->IsInBounds(i)){
                syst->SetParProp(i, curr_val);
                momentum[k] = -momentum[k];
            }
            offset += syst->GetNumParams();
        }
    }

    // Now we reweight
    for(auto& samp: samples){
        samp->Reweight();
    }
    

    // Now we update the gradient
    UpdateGradU();
    UpdateMomentum();
}

double HamiltonianMC::AcceptanceProbability(){
    // First we need to calculate the hamiltonian
    logLProp = 0.0;
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        // Could throw the initial value here to do MCMC stability studies
        // Propose the steps for the systematics
        systematics[s]->ProposeStep();

        // Get the likelihood from the systematics
        syst_llh[s] = systematics[s]->GetLikelihood();
        logLProp += syst_llh[s];

    }
    if (logLProp >= M3::_LARGE_LOGL_)
    {
        return 0.0;
    }



    for(size_t s=0; s<samples.size(); ++s){
        sample_llh[s] = samples[s]->GetLikelihood();
        logLProp += sample_llh[s];
    }


    // Now we have the LLH we can grab the kinetic energy
    double ke_prop = 0;
    for(auto& m : momentum){
        ke_prop+=0.5*std::pow(m, 2);
    }

    // Finally we can get our proposed hamiltonian
    double hamiltonian_prop = logLProp+ke_prop;

    // and our initial hamiltonian
    double hamiltonian_init = logLCurr+ke_init;

    // Now we can get the acceptance probability
    return std::min(1.0, std::exp(hamiltonian_init-hamiltonian_prop));
}