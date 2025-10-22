#include "Fitters/MCMCBase.h"

// *************************
// Initialise the manager and make it an object of mcmc class
// Now we can dump manager settings to the output file
MCMCBase::MCMCBase(manager *man) : FitterBase(man) {
// *************************
    Init();
}


void MCMCBase::Init(){
    // Beginning step number
    stepStart = 0;

    // Starting parameters should be thrown
    out_of_bounds = false;
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

void MCMCBase::PrepareOutput()
{
    FitterBase::PrepareOutput();

    // We also need to resize our proposal vectors
    sample_llh_prop.resize(samples.size());
    syst_llh_prop.resize(systematics.size());
}

// *******************
// Run the Markov chain with all the systematic objects added
void MCMCBase::RunMCMC() {
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
    for(size_t s=0; s < systematics.size(); ++s){    
        sample_llh[s] = sample_llh_prop[s];
    }
    for(size_t s=0; s < systematics.size(); ++s){
        syst_llh = syst_llh_prop;
    }

    // Begin MCMC
    const auto StepEnd = stepStart + chainLength;

#ifdef MPIENABLED
    // We send a signal to all other fitters that we are ready to start MCMC
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    for (step = stepStart; step < StepEnd; ++step)
    {
        DoMCMCStep();
    }
    // Save all the MCMC output
    SaveOutput();

    // Process MCMC
    ProcessMCMC();
}



// *******************
void MCMCBase::DoMCMCStep() {
// *******************
    /// Starts step timer, prints progress
    PreStepProcess();
    /// Step proposal, acceptance etc

    #ifdef MPIENABLED
    // We send a signal to all other fitters that we are ready to propose a step
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

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
        PrintProgress();
        /// Prevent thread conflcits
        #ifdef MPIENABLED
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
    }
}

// *************************
void MCMCBase::PostStepProcess() {
// *************************
    stepClock->Stop();
    stepTime = stepClock->RealTime();

    // Write step to output tree
    outTree->Fill();

    // Do Adaptive MCMC
    AdaptiveStep();

    if (step % auto_save == 0){
        outTree->AutoSave();
    }
}

// *******************
// Print the fit output progress
void MCMCBase::PrintProgress() {
// *******************
    MACH3LOG_INFO("Step:\t{}/{}, current: {:.2f}, proposed: {:.2f}", step - stepStart, chainLength, logLCurr, logLProp);
    MACH3LOG_INFO("Accepted/Total steps: {}/{} = {:.2f}", accCount, step - stepStart, static_cast<double>(accCount) / static_cast<double>(step - stepStart));
    #ifdef MPIENABLED
    // Print this info for all other MPI ranks
    if(mpi_rank==0){
        for(int i = 1; i < n_procs; ++i){
            double mpi_info[3] = {0.0, 0.0, 0.0};
            MPI_Recv(mpi_info, 3, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int rank = static_cast<int>(mpi_info[0]);
            int acc_steps = static_cast<int>(mpi_info[1]);
            int total_steps = static_cast<int>(mpi_info[2]);
            double acc_rate = 0.0;
            if(total_steps > 0){
                acc_rate = static_cast<double>(acc_steps) / static_cast<double>(total_steps);
            }
            MACH3LOG_INFO("Rank {}: Accepted/Total steps: {}/{} = {:.2e}", rank, acc_steps, total_steps, acc_rate);
        }
    }
    else{
        double mpi_info[3] = {static_cast<double>(mpi_rank), static_cast<double>(accCount), static_cast<double>(step - stepStart)};
        MPI_Send(mpi_info, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        return;
    }
    #endif

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
void MCMCBase::StartFromPreviousFit(const std::string &FitName) {
// *******************
    // Use base class
    FitterBase::StartFromPreviousFit(FitName);

    // For MCMC we also need to set stepStart
    TFile *infile = new TFile(FitName.c_str(), "READ");
    TTree *posts = infile->Get<TTree>("posteriors");
    unsigned int step_val = 0;

    posts->SetBranchAddress("step", &step_val);
    posts->GetEntry(posts->GetEntries() - 1);

    stepStart = step_val;
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
    #ifdef DEBUG
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
        syst_llh[s] = syst_llh_prop[s];
    }
    // Loop over samples and swap
    for( size_t i = 0; i < samples.size(); ++i)
    {
        sample_llh[i] = sample_llh_prop[i];
    }
}
