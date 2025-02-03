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
  chainLength = fitMan->raw()["General"]["MCMC"]["NSteps"].as<double>();

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

// *************************
// Load starting positions from the end of a previous chain
void mcmc::ReadParsFromFile(std::string file) {
// *************************
  MACH3LOG_INFO("MCMC getting starting position from {}", file);

  TFile *infile = new TFile(file.c_str(), "READ");
  TTree *posts = (TTree*)infile->Get("posteriors");
  TObjArray* brlis = (TObjArray*)posts->GetListOfBranches();
  int nbr = brlis->GetEntries();
  TString* branch_names = new TString[nbr];
  double* branch_vals = new double[nbr];

  for (int i = 0; i < nbr; ++i) {
    TBranch *br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();
    branch_names[i] = bname;
    std::cout << " * Loading " << bname << std::endl;
    posts->SetBranchAddress(branch_names[i], &branch_vals[i]);
  }

  posts->GetEntry(posts->GetEntries()-1);

  delete[] branch_names;
  delete[] branch_vals;
  infile->Close();
  delete infile;
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

    if (osc) {
      osc->acceptStep();
    }
    
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

  bool multicanonical = fitMan->raw()["General"]["MCMC"]["Multicanonical"].as<bool>();
  double multicanonicalBeta = fitMan->raw()["General"]["MCMC"]["MulticanonicalBeta"].as<double>();

  // Propose steps for the oscillation handlers
  if (osc) {
    osc->proposeStep();

    // Now get the likelihoods for the oscillation
    osc_llh = osc->GetLikelihood();

    if (multicanonical){
      std::string deltaCPParamNumber = osc->GetParName(1);

      double delta_cp_value = osc->getParProp(1);

      std::cout << "Delta CP value: " << delta_cp_value << std::endl;

      osc_llh = osc_llh * pow(GetMulticanonicalWeight(delta_cp_value),1-multicanonicalBeta);

    }
    std::cout << "Oscillation likelihood: " << osc_llh << std::endl;

    // Add the oscillation likelihoods to the reconfigure likelihoods
    llh += osc_llh;

    #ifdef DEBUG
    if (debug) debugFile << "LLH for oscillation handler: " << llh << std::endl;
    #endif
  }

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

  // Check if we've hit a boundary in the systematics
  // In this case we can save time by not having to reconfigure the simulation
  if (llh >= _LARGE_LOGL_) {
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
      sample_llh[i] = _LARGE_LOGL_;
      #ifdef DEBUG
      if (debug) debugFile << "LLH after REJECT sample " << i << " " << llh << std::endl;
      #endif
    }
  }

  // Save the proposed likelihood (class member)
  logLProp = llh;
}

double mcmc::GetMulticanonicalWeight(double deltacp){
  double delta_cp_value = deltacp;
  double delta_cp_log_likelihood;

  // calculate log likelihood according to gaussian around -1/2pi
  // #####proof of concept#####, this need to be changed to a sensible prior for dcp
  delta_cp_log_likelihood = -TMath::Log(TMath::Exp(-pow(delta_cp_value + 0.5*TMath::Pi(),2)));

  return delta_cp_log_likelihood;
}

// *******************
// Print the fit output progress
void mcmc::PrintProgress() {
// *******************

  MACH3LOG_INFO("Step:\t{}/{}, current: {:.2f}, proposed: {:.2f}", step - stepStart, chainLength, logLCurr, logLProp);
  MACH3LOG_INFO("Accepted/Total steps: {}/{} = {:.2f}", accCount, step - stepStart, static_cast<double>(accCount) / static_cast<double>(step - stepStart));

  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it) {
    if (std::string((*it)->getName()) == "xsec_cov") {
      MACH3LOG_INFO("Cross-section parameters: ");
      (*it)->printNominalCurrProp();
    }
  }
  #ifdef DEBUF
  if (debug) {
    debugFile << "\n-------------------------------------------------------" << std::endl;
    debugFile << "Step:\t" << step + 1 << "/" << chainLength << "  |  current: " << logLCurr << " proposed: " << logLProp << std::endl;
  }
  #endif
}
