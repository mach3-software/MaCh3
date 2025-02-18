#include "mcmc.h"

// *************************
// Initialise the manager and make it an object of mcmc class
// Now we can dump manager settings to the output file
mcmc::mcmc(manager *man) : FitterBase(man) {
// *************************
  // Beginning step number
  stepStart = 0;

  // Starting parameters should be thrown
  reject = false;
  chainLength = fitMan->raw()["General"]["MCMC"]["NSteps"].as<unsigned>();

  AnnealTemp = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["AnnealTemp"], -999);
  if(AnnealTemp < 0) anneal = false;
  else
  {
    MACH3LOG_INFO("Enabling simulated annealing with T = {}", AnnealTemp);
    anneal = true;
  }
}

// *************************
// Destructor: close the logger and output file
mcmc::~mcmc() {
// *************************

}

// **********************
// Do we accept the proposed step for all the parameters?
void mcmc::CheckStep() {
// **********************
  bool accept = false;

  // Set the acceptance probability to zero
  accProb = 0.0;

  // Calculate acceptance probability
  if (anneal) accProb = TMath::Min(1.,TMath::Exp( -(logLProp-logLCurr) / (TMath::Exp(-step/AnnealTemp)))); 
  else accProb = TMath::Min(1., TMath::Exp(logLCurr-logLProp));

  // Get the random number
  double fRandom = random->Rndm();

  // Do the accept/reject
  if (fRandom <= accProb) {
    accept = true;
    ++accCount;
  } else {
    accept = false;
  }

  #ifdef DEBUG
  if (debug) debugFile << " logLProp: " << logLProp << " logLCurr: " << logLCurr << " accProb: " << accProb << " fRandom: " << fRandom << std::endl;
  #endif

  // Update all the handlers to accept the step
  if (accept && !reject) {
    logLCurr = logLProp;

    // Loop over systematics and accept
    for (size_t s = 0; s < systematics.size(); ++s) {
      systematics[s]->acceptStep();
    }
  }

  stepClock->Stop();
  stepTime = stepClock->RealTime();

  // Write step to output tree
  outTree->Fill();
}

// *******************
// Run the Markov chain with all the systematic objects added
void mcmc::runMCMC() {

  // *******************

  // Set up the multicanonical weight on the delta_cp parameter
  multicanonical = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["Multicanonical"],false);
  MACH3LOG_INFO("Multicanonical Method: {}", multicanonical);
  if (multicanonical) {
    multicanonicalBeta = fitMan->raw()["General"]["MCMC"]["MulticanonicalBeta"].as<double>();
    MACH3LOG_INFO("Setting multicanonical beta to {}", multicanonicalBeta);
    MACH3LOG_INFO("Looping over systematics to find delta_cp parameter");
    for (size_t s = 0; s < systematics.size(); s++) {
      MACH3LOG_INFO("Systematic: {}", systematics[static_cast<int>(s)]->getName());
      if (systematics[static_cast<int>(s)]->getName() == "osc_cov") {
        oscCovVar = static_cast<int>(s);
        MACH3LOG_INFO("Found osc_cov systematic saving in variable {}", oscCovVar);
        for (int i = 0; i < systematics[static_cast<int>(s)]->GetNumParams(); i++) {
          MACH3LOG_INFO("Parameter: {}", systematics[static_cast<int>(s)]->GetParName(i));
          if (systematics[static_cast<int>(s)]->GetParName(i) == "delta_cp") {
            multicanonicalVar = i;
            MACH3LOG_INFO("Setting multicanonical weight on delta_cp parameter int {}",i);
          }
        }
      }
    }
  }

  // Save the settings into the output file
  SaveSettings();

  // Prepare the output branches
  PrepareOutput();

  // Reconfigure the samples, systematics and oscillation for first weight
  // ProposeStep sets logLProp
  ProposeStep();
  // Set the current logL to the proposed logL for the 0th step
  // Accept the first step to set logLCurr: this shouldn't affect the MCMC because we ignore the first N steps in burn-in
  logLCurr = logLProp;


  // Begin MCMC
  for (step = stepStart; step < stepStart+chainLength; ++step)
  {
    stepClock->Start();
    // Set the initial rejection to false
    reject = false;

    // Print 10 steps in total
    if ( (step-stepStart) % (chainLength/10) == 0) {
      PrintProgress();
    }

    // Propose current step variation and save the systematic likelihood that results in this step being taken
    // Updates logLProp
    ProposeStep();

    // Does the MCMC accept this step?
    CheckStep();

    // Auto save the output
    if (step % auto_save == 0) outTree->AutoSave();
  }

  // Save all the MCMC output
  SaveOutput();

  // Process MCMC
  ProcessMCMC();
}

// *******************
// Do the initial reconfigure of the MCMC
void mcmc::ProposeStep() {
// *******************
  // Initial likelihood
  double llh = 0.0;

  // Initiate to false
  reject = false;

  // Loop over the systematics and propose the initial step
  for (size_t s = 0; s < systematics.size(); ++s) {
    // Could throw the initial value here to do MCMC stability studies
    // Propose the steps for the systematics
    systematics[s]->proposeStep();

    // Get the likelihood from the systematics
    syst_llh[s] = systematics[s]->GetLikelihood();
    llh += syst_llh[s];

    #ifdef DEBUG
    if (debug) debugFile << "LLH after " << systematics[s]->getName() << " " << llh << std::endl;
    #endif
  }

  // if we're using the multicanonical method, we need to add the penalty to the likelihood now prior to the Large LLH check
  if (multicanonical){
    // MACH3LOG_INFO("LLH before multicanonical penalty: {}", llh);
    delta_cp_value = systematics[oscCovVar]->getParProp(multicanonicalVar);
    // MACH3LOG_INFO("Delta CP value: {}", delta_cp_value);
    multicanonical_penalty = pow(GetMulticanonicalWeight(delta_cp_value),1-multicanonicalBeta);
    // MACH3LOG_INFO("Multicanonical penalty: {}", multicanonical_penalty);
    llh += multicanonical_penalty;
    // MACH3LOG_INFO("LLH after multicanonical penalty: {}", llh);
  }

  // Check if we've hit a boundary in the systematics
  // In this case we can save time by not having to reconfigure the simulation
  if (llh >= M3::_LARGE_LOGL_) {
    reject = true;
    #ifdef DEBUG
    if (debug) debugFile << "Rejecting based on boundary" << std::endl;
    #endif
  }

  // Only reweight when we have a good parameter configuration
  // This speeds things up considerably because for every bad parameter configuration we don't have to reweight the MC
  if (!reject)
  {
    // Could multi-thread this
    // But since sample reweight is multi-threaded it's probably better to do that
    for (size_t i = 0; i < samples.size(); ++i)
    {
      samples[i]->reweight();
    }

    //DB for atmospheric event by event sample migration, need to fully reweight all samples to allow event passing prior to likelihood evaluation
    for (size_t i = 0; i < samples.size(); ++i) {
      // Get the sample likelihoods and add them
      sample_llh[i] = samples[i]->GetLikelihood();
      llh += sample_llh[i];
      #ifdef DEBUG
      if (debug) debugFile << "LLH after sample " << i << " " << llh << std::endl;
      #endif
    }

  // For when we don't have to reweight, set sample to madness
  } else {
    for (size_t i = 0; i < samples.size(); ++i) {
      // Set the sample_llh[i] to be madly high also to signify a step out of bounds
      sample_llh[i] = M3::_LARGE_LOGL_;
      #ifdef DEBUG
      if (debug) debugFile << "LLH after REJECT sample " << i << " " << llh << std::endl;
      #endif
    }
  }

  // Save the proposed likelihood (class member)
  logLProp = llh;
}

inline double mcmc::GetMulticanonicalWeight(double deltacp){

  constexpr double inv_sqrt_2pi = 0.3989422804014337;
  constexpr double neg_half_sigma = -0.5;

  double exp1 = std::exp(neg_half_sigma * (deltacp - TMath::Pi()) * (deltacp - TMath::Pi()));
  double exp2 = std::exp(neg_half_sigma * (deltacp) * (deltacp));
  double exp3 = std::exp(neg_half_sigma * (deltacp + TMath::Pi()) * (deltacp + TMath::Pi()));
  ///delta_cp_log_likelihood = -TMath::Log(TMath::Gaus(deltacp,TMath::Pi(),1,kTRUE)+TMath::Gaus(deltacp,0,1,kTRUE)+TMath::Gaus(deltacp,-TMath::Pi(),1,kTRUE));

  return -std::log(inv_sqrt_2pi * (exp1 + exp2 + exp3));
}

// *******************
// Print the fit output progress
void mcmc::PrintProgress() {
// *******************
  MACH3LOG_INFO("Step:\t{}/{}, current: {:.2f}, proposed: {:.2f}", step - stepStart, chainLength, logLCurr, logLProp);
  MACH3LOG_INFO("Accepted/Total steps: {}/{} = {:.2f}", accCount, step - stepStart, static_cast<double>(accCount) / static_cast<double>(step - stepStart));

  for (covarianceBase *cov : systematics) {
    if (cov->getName() == "xsec_cov") {
      MACH3LOG_INFO("Cross-section parameters: ");
      cov->printNominalCurrProp();
    }
  }
  #ifdef DEBUG
  if (debug) {
    debugFile << "\n-------------------------------------------------------" << std::endl;
    debugFile << "Step:\t" << step + 1 << "/" << chainLength << "  |  current: " << logLCurr << " proposed: " << logLProp << std::endl;
  }
  #endif
}

// *******************
void mcmc::StartFromPreviousFit(const std::string& FitName) {
// *******************
  // Use base class
  FitterBase::StartFromPreviousFit(FitName);

  // For MCMC we also need to set stepStart
  TFile *infile = new TFile(FitName.c_str(), "READ");
  TTree *posts = infile->Get<TTree>("posteriors");
  int step_val = -1;

  posts->SetBranchAddress("step",&step_val);
  posts->GetEntry(posts->GetEntries()-1);

  stepStart = step_val + 1;
  infile->Close();
  delete infile;
}
