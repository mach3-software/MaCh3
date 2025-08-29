#include "Fitters/MR2T2.h"

// *******************
void MR2T2::DoStep() {
// *******************
    ProposeStep();
    // Does the MCMC accept this step?
    accProb = AcceptanceProbability();

    if (IsStepAccepted(accProb) and !out_of_bounds)
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
