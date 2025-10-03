#include "Fitters/DEMCMC.h"


DEMCMC::DEMCMC(manager *const manager) : MR2T2(manager) {
    #ifdef MPIENABLED
    MACH3LOG_INFO("Using DEMCMC fitter");

    // Get the gamma parameter
    gamma = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["Lambda"], -1.0);
    
    // Gaussian with tiny variance
    scaling_matrix = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["ScalingMatrix"], 1e-5);
    #else
    MACH3LOG_ERROR("Using DEMCMC fitter without MPI... this is not implemented!");
    throw MaCh3Exception(__FILE__, __LINE__);
    #endif
}

#ifdef MPIENABLED
void DEMCMC::RunMCMC()
{
    curr_step_idx = std::vector<int>(systematics.size());
    // We can now loop over systs to assign the right amount of memory
    int idx = 0;
    for(size_t s=0; s < systematics.size(); ++s){
        curr_step_idx[s] = idx;
        for(int p = 0; p < systematics[s]->GetNumParams(); ++p){
            transfer_vec.push_back(0.0);
        }

        idx += systematics[s]->GetNumParams();
    }
    n_params = static_cast<int>(transfer_vec.size());

    all_params = std::vector<double>(n_params * n_procs);

    if(gamma<0){
        gamma = 2.38 / std::sqrt(2*n_params);
    }

    MACH3LOG_INFO("Running DEMCMC with g={} and g={}", gamma, scaling_matrix);

    MR2T2::RunMCMC();
}

void DEMCMC::ProposeStep(){
    // Propose MR2T2 step    
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        // Could throw the initial value here to do MCMC stability studies
        // Propose the steps for the systematics
        for (int p = 0; p < systematics[s]->GetNumParams(); ++p)
        {
            transfer_vec[curr_step_idx[s] + p] = systematics[s]->GetParCurr(p);
        }
    }

    // Gather all parameters from all processes
    MPI_Allgather(transfer_vec.data(), n_params, MPI_DOUBLE,
                  all_params.data(), n_params, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    // NOW we select two fitters from our fitter vector
    int random_a = random->Integer(n_procs);
    while (random_a == mpi_rank)
    {
        random_a = random->Integer(n_procs);
    }
    

    int random_b = random->Integer(n_procs);
    // A bit dodgy but we want to make sure we don't pick the same fitter twice
    while(random_b == random_a or random_b == mpi_rank){
        random_b = random->Integer(n_procs);
    }


    for (size_t s = 0; s < systematics.size(); ++s){
        for (int p = 0; p < systematics[s]->GetNumParams(); ++p)
        {
            // Get the two other steps
            // Now we check if it's fixed
            if (systematics[s]->IsParameterFixed(p)) continue;

            double step_a = all_params[random_a * n_params + curr_step_idx[s] + p];
            double step_b = all_params[random_b * n_params + curr_step_idx[s] + p];

            // Now do the DEMCMC step
            double perturbation = random->Gaus(0, scaling_matrix);
            transfer_vec[curr_step_idx[s] + p] += gamma * (step_a - step_b) + perturbation;
        }
        // systematics[s]->SpecialStepProposal();
    }
    /// We now have everything we need to do the DEMCMC step so we need to wait
    MPI_Barrier(MPI_COMM_WORLD);

    // Final loop
    double llh = 0.0;
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        // Set the proposed step
        for (int p = 0; p < systematics[s]->GetNumParams(); ++p)
        {
            systematics[s]->SetParProp(p, transfer_vec[curr_step_idx[s] + p]);
        }

        // Get the likelihood from the systematics
        syst_llh[s] = systematics[s]->GetLikelihood();
        llh += syst_llh[s];
#ifdef DEBUG
                if (debug)
                    debugFile << "LLH after " << systematics[s]->GetName() << " " << llh << std::endl;
#endif
            
    }


    // Check if we've hit a boundary in the systematics
    // In this case we can save time by not having to reconfigure the simulation
    if (llh >= M3::_LARGE_LOGL_)
    {
        out_of_bounds = true;
#ifdef DEBUG
        if (debug)
            debugFile << "Rejecting based on boundary" << std::endl;
#endif
    }

    // Only reweight when we have a good parameter configuration
    // This speeds things up considerably because for every bad parameter configuration we don't have to reweight the MC
    if (!out_of_bounds)
    {
        // Could multi-thread this
        // But since sample reweight is multi-threaded it's probably better to do that
        for (size_t i = 0; i < samples.size(); ++i)
        {   
            try{
                samples[i]->Reweight();
            }
            catch(...){
                out_of_bounds = true;
                llh = M3::_LARGE_LOGL_;
                break;
            }
        }

        // DB for atmospheric event by event sample migration, need to fully reweight all samples to allow event passing prior to likelihood evaluation
        for (size_t i = 0; i < samples.size(); ++i)
        {
            if(out_of_bounds) break;
            // Get the sample likelihoods and add them
            sample_llh[i] = samples[i]->GetLikelihood();
            llh += sample_llh[i];
#ifdef DEBUG
            if (debug)
                debugFile << "LLH after sample " << i << " " << llh << std::endl;
#endif
        }

        // For when we don't have to reweight, set sample to madness
    }
    else
    {
        for (size_t i = 0; i < samples.size(); ++i)
        {
            // Set the sample_llh[i] to be madly high also to signify a step out of bounds
            sample_llh[i] = M3::_LARGE_LOGL_;
#ifdef DEBUG
            if (debug)
                debugFile << "LLH after REJECT sample " << i << " " << llh << std::endl;
#endif
        }
    }
    // Save the proposed likelihood (class member)
    logLProp = llh;
}
#else
void DEMCMC::ProposeStep(){
    MACH3LOG_ERROR("Using DEMCMC fitter without MPI... this is not implemented!");
    throw MaCh3Exception(__FILE__, __LINE__);
}
#endif