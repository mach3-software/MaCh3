#include "Fitters/MCMCBase.h"

// *************************
// Initialise the Manager and make it an object of mcmc class
// Now we can dump Manager settings to the output file
MCMCBase::MCMCBase(Manager *man) : FitterBase(man) {
// *************************
    // Beginning step number
    stepStart = 0;

    // Starting parameters should be thrown
    out_of_bounds = false;
    chainLength = Get<unsigned>(fitMan->raw()["General"]["MCMC"]["NSteps"], __FILE__, __LINE__);
    if (chainLength < 10){
        MACH3LOG_ERROR("MCMC chain length must be at least 10 steps, otherwise this will result in a floating point exception.");
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    AnnealTemp = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["AnnealTemp"], -999, __FILE__ , __LINE__);
    if (AnnealTemp < 0)
        anneal = false;
    else
    {
        MACH3LOG_INFO("Enabling simulated annealing with T = {}", AnnealTemp);
        anneal = true;
    }
}


// *******************
// Run the Markov chain with all the systematic objects added
void MCMCBase::RunMCMC() {
// *******************
    // Multicanonical method toggle from yaml config
    multicanonical = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["Multicanonical"]["Enabled"], false, __FILE__, __LINE__);
    MACH3LOG_INFO("Multicanonical Method: {}", multicanonical);

    if (multicanonical) {
        /// Initialise the multicanonical handler
        multicanonicalHandler = std::make_unique<MulticanonicalMCMCHandler>();
    
        // Initialize the multicanonical handler with the systematics
        multicanonicalHandler->InitializeMulticanonicalHandlerConfig(fitMan, systematics);
        AlgorithmName += "_UmbrellaSampling"; // Append to the algorithm name
#ifdef MACH3_DEBUG
    // Enable debug output stream for multicanonical handler if debug is enabled
    multicanonicalHandler->setDebugStream(&debugFile, debug);
#endif
    }

    // Save the settings into the output file
    SaveSettings();

    // Prepare the output branches
    PrepareOutput();

    // Remove obsolete memory and make other checks before fit starts
    SanitiseInputs();

    // Print Progress before Propose Step
    PrintProgress(false);

    // Only propose step if we are running fresh chain. If we are running from previous chain then it's not needed and we can use usual pipeline
    if(stepStart == 0) {
        // Reconfigure the samples, systematics and oscillation for first weight
        // ProposeStep sets logLProp
        ProposeStep();

        // Initialise the value of the multicanonical parameter to the centres of the umbrellas
        if (multicanonical){
            multicanonicalHandler->InitializeMulticanonicalParams(systematics);
        }
        // Set the current logL to the proposed logL for the 0th step
        // Accept the first step to set logLCurr: this shouldn't affect the MCMC because we ignore the first N steps in burn-in
        logLCurr = logLProp;
    }

    // Begin MCMC
    const auto StepEnd = stepStart + chainLength;
    for (step = stepStart; step < StepEnd; ++step)
    {
        DoMCMCStep();
    }
    // Save all the MCMC output
    SaveOutput();

    // Process MCMC
    ProcessMCMC();

    // Save the adaptive MCMC
    for (const auto &syst : systematics)
    {
        if (syst->GetDoAdaption())
        {
            auto adaptive_handler = syst->GetAdaptiveHandler();
            adaptive_handler->SaveAdaptiveToFile(adaptive_handler->GetOutFileName(), syst->GetName(), true);
        }
    }
}

// *******************
void MCMCBase::DoMCMCStep() {
// *******************
    /// Starts step timer, prints progress
    PreStepProcess();
    /// Step proposal, acceptance etc
    DoStep();
    /// Tree filling etc.
    PostStepProcess();
}

// *******************
void MCMCBase::PreStepProcess() {
// *******************
    stepClock->Start();
    out_of_bounds = false;

    // Print 10 steps in total
    if ((step - stepStart) % (chainLength / 10) == 0)
    {
        CheckAcceptanceRates();
        PrintProgress();
    }
}

// *************************
void MCMCBase::PostStepProcess() {
// *************************
    //KS: Some version of ROOT keep spamming about accessing already deleted object which is wrong and not helpful...
    int originalErrorLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;

    stepClock->Stop();
    stepTime = stepClock->RealTime();

    // Write step to output tree
    outTree->Fill();

    // Do Adaptive MCMC
    AdaptiveStep();

    if (step % auto_save == 0){
        outTree->AutoSave();
    }
    gErrorIgnoreLevel = originalErrorLevel;
}

// *******************
// Print the fit output progress
void MCMCBase::PrintProgress(bool StepsPrint) {
// *******************
    if(StepsPrint) MACH3LOG_INFO("Step:\t{}/{}, current: {:.2f}, proposed: {:.2f}", step - stepStart, chainLength, logLCurr, logLProp);
    if(StepsPrint) MACH3LOG_INFO("Accepted/Total steps: {}/{} = {:.2f}", accCount, step - stepStart, static_cast<double>(accCount) / static_cast<double>(step - stepStart));

    for (size_t i = 0; i < samples.size(); ++i) {
        samples[i]->PrintRates();
    }

    for (ParameterHandlerBase *cov : systematics) {
        cov->PrintPreFitCurrPropValues();
    }
#ifdef MACH3_DEBUG
    if (debug)
    {
        debugFile << "\n-------------------------------------------------------" << std::endl;
        debugFile << "Step:\t" << step + 1 << "/" << chainLength << "  |  current: " << logLCurr << " proposed: " << logLProp << std::endl;
    }
#endif
}

// *******************
// Handles warnings for low acceptance rates
void MCMCBase::CheckAcceptanceRates() {
// *******************
    const auto StepEnd = stepStart + chainLength;

    // Skip check for reallllly early steps
    if(step - stepStart < std::min(10.0, static_cast<double>(StepEnd)/10)) {
        return;
    }
    // KS: Do not add "throw", this can crash CI or some short dummy chains for testing
    if(accCount==0 && step - stepStart > StepEnd/10) {
        MACH3LOG_CRITICAL("No steps were accepted in the MCMC chain after {} steps. Please check your  step sizes.", step - stepStart);
    }

    double acc_rate = static_cast<double>(accCount) / static_cast<double>(step - stepStart);

    if (acc_rate < 0.01) {
        MACH3LOG_WARN("Acceptance rate is very low: {:.4f}. This may indicate that the proposal distribution is not well-tuned. Consider reducing the step size.", acc_rate);
    } else if (acc_rate > 0.5) {
        MACH3LOG_WARN("Acceptance rate is very high: {:.4f}. This may indicate that the proposal distribution is too narrow. Consider increasing the step sizes.", acc_rate);
    }
}

// *******************
void MCMCBase::StartFromPreviousFit(const std::string &FitName) {
// *******************
    // Use base class
    FitterBase::StartFromPreviousFit(FitName);

    // For MCMC we also need to set stepStart
    TFile *infile = M3::Open(FitName, "READ", __FILE__, __LINE__);
    TTree *posts = infile->Get<TTree>("posteriors");
    unsigned int step_val = 0;

    posts->SetBranchAddress("step", &step_val);
    posts->GetEntry(posts->GetEntries() - 1);

    stepStart = step_val;
    // KS: Also update number of steps if using adaption
    for (unsigned int i = 0; i < systematics.size(); ++i) {
        if (systematics[i]->GetDoAdaption()) {
            systematics[i]->SetNumberOfSteps(stepStart);
        }
    }
    infile->Close();
    delete infile;
}

// *************************
void MCMCBase::AdaptiveStep() {
// *************************
    // Save the Adaptive output
    for (const auto &syst : systematics)
    {
        if (syst->GetDoAdaption()){
            syst->UpdateAdaptiveCovariance();
        }
    }
}

// *************************
bool MCMCBase::IsStepAccepted(const double acc_prob) {
// *************************
    // Get the random number
    const double fRandom = random->Rndm();
    // Do the accept/reject
    #ifdef MACH3_DEBUG
    debugFile << " logLProp: " << logLProp << " logLCurr: " << logLCurr << " acc_prob: " << acc_prob << " fRandom: " << fRandom << std::endl;
    #endif

    if (fRandom > acc_prob)
    {
        // Reject
        return false;
    }
    return true;
}

// *************************
void MCMCBase::AcceptStep() {
// *************************
    ++accCount;
    logLCurr = logLProp;

    // Loop over systematics and accept
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        systematics[s]->AcceptStep();
    }
}
