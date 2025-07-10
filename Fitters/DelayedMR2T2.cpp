#include "Fitters/DelayedMR2T2.h"

DelayedMR2T2::DelayedMR2T2(manager * const manager) : MR2T2(manager) {
    
    decay_rate = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["DecayRate"], 0.1);
    number_of_iterations = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["NRejections"], 2);
    
    MinLogLikelihood = M3::_LARGE_LOGL_;
    MACH3LOG_INFO("Using Delayed MCMC with decay rate: {} and {} allowed rejections", decay_rate, number_of_iterations);
}

void DelayedMR2T2::ScaleSystematics(double scale)
{
    for (auto &syst : systematics)
    {
        syst->SetStepScale(scale*syst->GetGlobalStepScale(), false);
    }
}

void DelayedMR2T2::ResetSystScale()
{
    // if(step_syst_scale.size() != systematics.size()){
    //     return;
    // }
    for (int i = 0; i < static_cast<int>(systematics.size()); ++i)
    {
        systematics[i]->SetStepScale(start_step_scale[i], false);
    }
}

//**********************************************
void DelayedMR2T2::StoreCurrentStep()
// *********************************************
{
    /*
    Stores values from  step
    */
    // Ensure vector is correctly sized
    current_step_vals.resize(systematics.size());
    start_step_scale.resize(systematics.size());

    // Fill with old proposed step, need to copy to avoid using the same step twice
    for(int i=0; i<static_cast<int>(systematics.size()); ++i){
        current_step_vals[i] = std::vector<double>(systematics[i]->GetParCurrVec());
        start_step_scale[i] = systematics[i]->GetGlobalStepScale();
    }
}

double DelayedMR2T2::AcceptanceProbability()
{
    // From https://arxiv.org/pdf/2010.04190 and DRAM Matlab implementation

    double numerator = std::max(0.0, std::exp(MinLogLikelihood - logLProp) - 1.0);
    double denominator = std::exp(MinLogLikelihood - logLCurr) - 1.0;
    if (denominator <= 0.0) return 1.0;

    // Large enough that the the minus 1 will make NO difference so the max term cancels
    if(std::isinf(numerator) or std::isinf(denominator) ){
        return MR2T2::AcceptanceProbability();
    }


    return std::min(numerator / denominator, 1.0);
}

void DelayedMR2T2::DoStep()
{
    // Set the initial rejection to false
    StoreCurrentStep();
    
    bool accepted = false;
    MinLogLikelihood = M3::_LARGE_LOGL_;

    for(int i=0; i<number_of_iterations; ++i)
    {
        ProposeStep();

        // Small hack to ensure we can "leapfrog"
        for (auto &syst : systematics)
        {
            syst->AcceptStep();
        }

        if(out_of_bounds or logLProp>MinLogLikelihood){
            continue;
        }


        if(i==0){
            accProb = MR2T2::AcceptanceProbability();
        }
        else{
            accProb = AcceptanceProbability();
        }
        if(IsStepAccepted(accProb)){
            AcceptStep();
            accepted = true;
            // Update
            break;
        }
        // Temporary, means we can "leapfrog"

        ScaleSystematics(decay_rate);
        MinLogLikelihood = logLProp;
    }

    // If we reject we need to reset everything [this is probably a tad too slow...
    if(!accepted){
        for(int i=0; i<static_cast<int>(systematics.size()); ++i){
            for(int j=0; j<static_cast<int>(systematics[i]->GetParCurrVec().size()); ++j){
                systematics[i]->SetParCurrProp(j, current_step_vals[i][j]);
            }
        }
    }
    ResetSystScale();
}
