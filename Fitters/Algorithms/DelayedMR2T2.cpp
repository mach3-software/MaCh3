#include "Fitters/Algorithms/DelayedMR2T2.h"

// *************************
DelayedMR2T2::DelayedMR2T2(Manager* const FitManager) : MR2T2(FitManager) {
// *************************
    AlgorithmName = "DelayedMR2T2";
    auto Config = fitMan->raw();
    // Step scale for delayed step = step scale prev delay * decay_rate
    decay_rate = GetFromManager<double>(Config["General"]["MCMC"]["DecayRate"], 0.1, __FILE__, __LINE__);
    // Maximum number of rejections
    max_rejections = GetFromManager<int>(Config["General"]["MCMC"]["MaxRejections"], 1, __FILE__, __LINE__);
    // How big is the first step scale. This scales parameter step by initial_scale * <step scale in ParameterHandler>
    initial_scale = GetFromManager<int>(Config["General"]["MCMC"]["InitialScale"], 1, __FILE__, __LINE__);
    // On delayed on out of bounds steps. This saves some compute time and possibly regains efficiency
    delay_on_oob_only = GetFromManager<bool>(Config["General"]["MCMC"]["DelayOnlyOutBounds"], false, __FILE__, __LINE__);
    // Allow to intrude another delayed rejection based on defined here probability
    delay_probability = GetFromManager<double>(Config["General"]["MCMC"]["DelayProbability"], 1.0, __FILE__, __LINE__);

    if(delay_probability > 1.0){
        MACH3LOG_WARN("Probability of delaying steps > 1, setting to 1");
        delay_probability = 1.0;
    }
    if(delay_probability <= 0.0){
        MACH3LOG_CRITICAL("Probability of delaying steps <= 0, setting to 0, delay will no longer have an impact!");
        delay_probability = 0.0;
        // Since we're not delaying at all set max_rejections = 0 to save time
        max_rejections = 0;
    }
    MACH3LOG_INFO("Using Delayed MCMC with decay rate: {} and {} allowed rejections", decay_rate, max_rejections);
}

// *************************
void DelayedMR2T2::ScaleSystematics(const double scale) {
// *************************
    for (auto &syst : systematics)
    {
        syst->SetStepScale(scale*syst->GetGlobalStepScale(), false);
    }
}

// *************************
void DelayedMR2T2::ResetSystScale() {
// *************************
    for (int i = 0; i < static_cast<int>(systematics.size()); ++i)
    {
        systematics[i]->SetStepScale(start_step_scale[i], false);
    }
}

// *************************
void DelayedMR2T2::PrepareOutput() {
// *************************
    FitterBase::PrepareOutput();

    for(size_t i = 0; i < systematics.size(); ++i){
        if(systematics[i]->GetDoAdaption()){
            if(systematics[i]->GetAdaptiveHandler()->GetUseRobbinsMonro()){
                MACH3LOG_ERROR("Right now Robbins-Monro doesn't work with Delayed MCMC");
                MACH3LOG_ERROR("Usual adaptive works fine though");
                MACH3LOG_ERROR("Ask Dr Wallace for details");
                throw MaCh3Exception(__FILE__ , __LINE__ );
            }
        }
    }
    // Store delayed specific settings
    outTree->Branch("DelayedStep", &accepted_delayed, "DelayedStep/B");
}

//**********************************************
void DelayedMR2T2::StoreCurrentStep() {
// *********************************************
    // Stores values from  step
    // Ensure vector is correctly sized
    current_step_vals.resize(systematics.size());
    start_step_scale.resize(systematics.size());

    // Fill with old proposed step, need to copy to avoid using the same step twice
    for(int i=0; i<static_cast<int>(systematics.size()); ++i){
        current_step_vals[i] = std::vector<double>(systematics[i]->GetParCurrVec());
        start_step_scale[i] = systematics[i]->GetGlobalStepScale();
    }
}

// *************************
double DelayedMR2T2::AcceptanceProbability() {
// *************************
    // From https://arxiv.org/pdf/2010.04190 and DRAM Matlab implementation
    const double numerator = std::max(0.0, std::exp(MinLogLikelihood - logLProp) - 1.0);
    const double denominator = std::exp(MinLogLikelihood - logLCurr) - 1.0;
    if (denominator <= 0.0) return 1.0;

    // Large enough that the the minus 1 will make NO difference so the max term cancels
    if(std::isinf(numerator) || std::isinf(denominator) ){
        return MR2T2::AcceptanceProbability();
    }

    return std::min(numerator / denominator, 1.0);
}

// *************************
void DelayedMR2T2::DoStep() {
// *************************
    // Set the initial rejection to false
    StoreCurrentStep();

    // Apply initial scale
    ScaleSystematics(initial_scale);

    bool accepted = false;
    bool delay = false;
    accepted_delayed = false;

    MinLogLikelihood = M3::_LARGE_LOGL_;

    for(int i = 0; i<= max_rejections; ++i)
    {
        ProposeStep();

        // Small hack to ensure we can "leapfrog"
        for (auto &syst : systematics)
        {
            syst->AcceptStep();
        }

        if(out_of_bounds || logLProp > MinLogLikelihood){
            continue;
        }

        if(i == 0){
            accProb = MR2T2::AcceptanceProbability();
        }
        else{
            accProb = AcceptanceProbability();
            delay = true;
        }
        if(IsStepAccepted(accProb)){
            AcceptStep();
            accepted = true;

            // Keep track of if step is delayed + accepted for storage
            accepted_delayed = delay;

            // Update
            break;
        }
        /// OOB Check
        if(!out_of_bounds && delay_on_oob_only) {
            // If we are not out of bounds and only delaying on out of bounds steps
            // We can skip the delay
            break;
        }
        
        /// Probabilistic condition
        if (!ProbabilisticDelay()) {
            break;
        }


        // Temporary, means we can "leapfrog"
        ScaleSystematics(decay_rate);
        MinLogLikelihood = logLProp;
    }

    // If we reject we need to reset everything [this is probably a tad too slow...
    // especially in PCA...
    if(!accepted){
        for(int i=0; i<static_cast<int>(systematics.size()); ++i){
            for(int j=0; j<static_cast<int>(systematics[i]->GetParCurrVec().size()); ++j){
                systematics[i]->SetParCurrProp(j, current_step_vals[i][j]);
            }
        }
    }
    ResetSystScale();
}

// *************************
bool DelayedMR2T2::ProbabilisticDelay() const {
// *************************
    // We can delay probabilistically
    return (M3::rand::Uniform() > delay_probability);
}
