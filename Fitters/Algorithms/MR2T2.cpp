#include "Fitters/Algorithms/MR2T2.h"


// ***************
MR2T2::MR2T2(Manager *man) : MCMCBase(man) {
// ***************
    AlgorithmName = "MR2T2";
}

// *******************
void MR2T2::DoStep() {
// *******************
    ProposeStep();
    // Does the MCMC accept this step?
    accProb = AcceptanceProbability();

    if (!out_of_bounds && IsStepAccepted(accProb))
    {
        AcceptStep();
    }
}

// *******************
// Do the initial reconfigure of the MCMC
void MR2T2::ProposeStep() {
// *******************
    // Initial likelihood
    out_of_bounds = false;
    
    double llh = 0.0;

    // Loop over the systematics and propose the initial step
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        // Could throw the initial value here to do MCMC stability studies
        // Propose the steps for the systematics
        systematics[s]->ProposeStep();

        // Get the likelihood from the systematics
        syst_llh[s] = systematics[s]->GetLikelihood();
        llh += syst_llh[s];

#ifdef MACH3_DEBUG
        if (debug)
            debugFile << "LLH after " << systematics[s]->GetName() << " " << llh << std::endl;
#endif
    }

    // if we're using the multicanonical method, we need to add the penalty to the
    // likelihood now prior to the Large LLH check
    if (multicanonical) {
    // get the proposed value of delta_cp and apply the multicanonical penalty,
    // weighting it using the beta value to increase or decrease the strength of
    // the penalty
        double delta_cp_value = systematics[multicanonicalHandler->oscCovVar]->GetParProp(multicanonicalHandler->multicanonicalVar);
        double delm23_value = systematics[multicanonicalHandler->oscCovVar]->GetParProp(multicanonicalHandler->multicanonicalVar_dm23);
        double multicanonical_penalty = multicanonicalHandler->GetMulticanonicalWeight(delta_cp_value, delm23_value);

#ifdef MACH3_DEBUG
        // Print the multicanonical penalty and the delta_cp and delm23 values to the debug file
        if (debug) debugFile << " delta_cp: " << delta_cp_value << " delm23: " << delm23_value << " multicanonical_penalty: " << multicanonical_penalty << std::endl;
#endif

        llh += multicanonical_penalty;

        MACH3LOG_DEBUG("Delta CP value: {}", delta_cp_value);
        MACH3LOG_DEBUG("Multicanonical penalty: {}", multicanonical_penalty);
        MACH3LOG_DEBUG("LLH after multicanonical penalty: {}", llh);
    }

    // Check if we've hit a boundary in the systematics
    // In this case we can save time by not having to reconfigure the simulation
    if (llh >= M3::_LARGE_LOGL_)
    {
        out_of_bounds = true;
#ifdef MACH3_DEBUG
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
#ifdef MACH3_DEBUG
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
#ifdef MACH3_DEBUG
            if (debug)
                debugFile << "LLH after REJECT sample " << i << " " << llh << std::endl;
#endif
        }
    }
    // Save the proposed likelihood (class member)
    logLProp = llh;
}

// **********************
// Do we accept the proposed step for all the parameters?
double MR2T2::AcceptanceProbability() {
// **********************
    // Set the acceptance probability to zero
    double acc_prob = 0.0;

    // Calculate acceptance probability
    if (anneal)
        acc_prob = std::min(1., std::exp(-(logLProp - logLCurr) / (std::exp(-step / AnnealTemp))));
    else
        acc_prob = std::min(1., std::exp(logLCurr - logLProp));

    return acc_prob;
}
