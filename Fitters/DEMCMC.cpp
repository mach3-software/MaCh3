#include "Fitters/DEMCMC.h"

DEMCMC::DEMCMC(manager *manager) : MR2T2(manager) {
    MACH3LOG_ERROR("Using DEMCMC fitter without MPI... this is not implemented!");
    throw MaCh3Exception(__FILE__, __LINE__);
}


#ifdef MPIENABLED
DEMCMC::DEMCMC(manager *const manager, int mpi_rank_) : MR2T2(manager, mpi_rank_) {
    MACH3LOG_INFO("Using DEMCMC fitter");

    // Get the gamma parameter
    gamma = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["Lambda"], 1.0);

    // Gaussian with tiny variance
    scaling_matrix = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["ScalingMatrix"], 1e-6);
}
#endif

void DEMCMC::RunMCMC()
{
    curr_step = std::vector<std::vector<double>>(systematics.size());
    // We can now loop over systs to assign the right amount of memory
    for(size_t s=0; s < systematics.size(); ++s){
        curr_step[s] = std::vector<double>(systematics[s]->GetNumParams(), 0.0);
    }

    MR2T2::RunMCMC();
}

void DEMCMC::ProposeStep(){
    // Propose MR2T2 step    
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        // Could throw the initial value here to do MCMC stability studies
        // Propose the steps for the systematics
        systematics[s]->ProposeStep();
        for (int p = 0; p < systematics[s]->GetNumParams(); ++p)
        {
            curr_step[s][p] = systematics[s]->GetParProp(p);
        }
    }

    // We send a signal to all other fitters that we are ready to do a DEMCMC step
    // We then wait for all other fitters to be ready
    #ifdef MPIENABLED
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    // NOW we select two fitters from our fitter vector
    int random_a = random->Integer(GetNFitters());
    int random_b = random->Integer(GetNFitters());
    // A bit dodgy but we want to make sure we don't pick the same fitter twice
    while(random_b == random_a){
        random_b = random->Integer(GetNFitters());
    }


    for (size_t s = 0; s < systematics.size(); ++s){
        auto a_vec = connectedFitters[random_a]->GetSystObjs()[s]->GetParPropVec();
        auto b_vec = connectedFitters[random_b]->GetSystObjs()[s]->GetParPropVec();
        if(a_vec.size() != b_vec.size()){
            MACH3LOG_ERROR("Connected fitters have different number of parameters this should not happen");
            throw MaCh3Exception(__FILE__, __LINE__);
        }


        for (int p = 0; p < systematics[s]->GetNumParams(); ++p)
        {
            // Get the two other steps
            double step_a = connectedFitters[random_a]->GetSystObjs()[s]->GetParProp(p);
            double step_b = connectedFitters[random_b]->GetSystObjs()[s]->GetParProp(p);

            // Now do the DEMCMC step
            double perturbation = random->Gaus(0, scaling_matrix);
            curr_step[s][p] += gamma * (step_a - step_b) + perturbation;
        }
    }
    /// We now have everything we need to do the DEMCMC step so we need to wait
    #ifdef MPIENABLED
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    // Final loop
    double llh = 0.0;
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        // Set the proposed step
        for (int p = 0; p < systematics[s]->GetNumParams(); ++p)
        {
            systematics[s]->SetParProp(p, curr_step[s][p]);
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
            samples[i]->Reweight();
        }

        // DB for atmospheric event by event sample migration, need to fully reweight all samples to allow event passing prior to likelihood evaluation
        for (size_t i = 0; i < samples.size(); ++i)
        {
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