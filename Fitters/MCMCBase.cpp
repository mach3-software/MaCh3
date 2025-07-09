#include "MCMCBase.h"

#include "mcmc.h"

// *************************
// Initialise the manager and make it an object of mcmc class
// Now we can dump manager settings to the output file
MCMCBase::MCMCBase(manager *man) : FitterBase(man)
{
    // *************************
    // Beginning step number
    stepStart = 0;

    // Starting parameters should be thrown
    reject = false;
    chainLength = Get<unsigned>(fitMan->raw()["General"]["MCMC"]["NSteps"], __FILE__, __LINE__);

    AnnealTemp = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["AnnealTemp"], -999);
    if (AnnealTemp < 0)
        anneal = false;
    else
    {
        MACH3LOG_INFO("Enabling simulated annealing with T = {}", AnnealTemp);
        anneal = true;
    }
}

bool MCMCBase::AcceptStep(const double acc_prob){
    // Get the random number
    const double fRandom = random->Rndm();
    // Do the accept/reject
    if (fRandom > accProb)
    {
        // Reject
        return false;
    }
    ++accCount;
    accept = true;

    logLCurr = logLProp;

    // Loop over systematics and accept
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        systematics[s]->AcceptStep();
    }

    return accept;
}

// *******************
// Run the Markov chain with all the systematic objects added
void MCMCBase::RunMCMC()
{
    // *******************
    // Save the settings into the output file
    SaveSettings();

    // Prepare the output branches
    PrepareOutput();

    // Remove obsolete memory and make other checks before fit starts
    SanitiseInputs();

    // Reconfigure the samples, systematics and oscillation for first weight
    // ProposeStep sets logLProp
    ProposeStep();
    // Set the current logL to the proposed logL for the 0th step
    // Accept the first step to set logLCurr: this shouldn't affect the MCMC because we ignore the first N steps in burn-in
    logLCurr = logLProp;

    // Begin MCMC
    const auto StepEnd = stepStart + chainLength;
    for (step = stepStart; step < StepEnd; ++step)
    {
        stepClock->Start();
        DoMCMCStep()
        // Auto save the output
        if (step % auto_save == 0)
            outTree->AutoSave();
    }

    // Save all the MCMC output
    SaveOutput();

    // Process MCMC
    ProcessMCMC();
}


void MCMCBase::DoMCMCStep(){
    // Set the initial rejection to false
    reject = false;

    // Print 10 steps in total
    if ((step - stepStart) % (chainLength / 10) == 0)
    {
        PrintProgress();
    }

    // Propose current step variation and save the systematic likelihood that results in this step being taken
    // Updates logLProp
    ProposeStep();
    // Does the MCMC accept this step?
    if (!reject)
    {
        acc_prob = CheckStep();
        AcceptStep(acc_prob);
    }

    stepClock->Stop();
    stepTime = stepClock->RealTime();

    // Write step to output tree
    outTree->Fill();

    // Do Adaptive MCMC
    AdaptiveStep()
}

// *******************
// Print the fit output progress
void MCMCBase::PrintProgress()
{
    // *******************
    MACH3LOG_INFO("Step:\t{}/{}, current: {:.2f}, proposed: {:.2f}", step - stepStart, chainLength, logLCurr, logLProp);
    MACH3LOG_INFO("Accepted/Total steps: {}/{} = {:.2f}", accCount, step - stepStart, static_cast<double>(accCount) / static_cast<double>(step - stepStart));

    for (ParameterHandlerBase *cov : systematics)
    {
        cov->PrintNominalCurrProp();
    }
#ifdef DEBUG
    if (debug)
    {
        debugFile << "\n-------------------------------------------------------" << std::endl;
        debugFile << "Step:\t" << step + 1 << "/" << chainLength << "  |  current: " << logLCurr << " proposed: " << logLProp << std::endl;
    }
#endif
}

// *******************
void MCMCBase::StartFromPreviousFit(const std::string &FitName)
{
    // *******************
    // Use base class
    FitterBase::StartFromPreviousFit(FitName);

    // For MCMC we also need to set stepStart
    TFile *infile = new TFile(FitName.c_str(), "READ");
    TTree *posts = infile->Get<TTree>("posteriors");
    int step_val = -1;

    posts->SetBranchAddress("step", &step_val);
    posts->GetEntry(posts->GetEntries() - 1);

    stepStart = step_val + 1;
    // KS: Also update number of steps if using adaption
    for (unsigned int i = 0; i < systematics.size(); ++i)
    {
        if (systematics[i]->GetDoAdaption())
        {
            systematics[i]->SetNumberOfSteps(step_val);
        }
    }
    infile->Close();
    delete infile;
}


void MCMCBase::AdaptiveStep(){
    // Save the Adaptive output
    for (const auto &syst : systematics)
    {
        if (syst->GetDoAdaption())
        {
            syst->GetAdaptiveHandler()->SaveAdaptiveToFile(syst->GetAdaptiveHandler()->GetOutFileName(), syst->GetName(), true);
        }
    }
}