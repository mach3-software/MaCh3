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
  chainLength = Get<unsigned>(fitMan->raw()["General"]["MCMC"]["NSteps"], __FILE__, __LINE__);

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
  // Multicanonical method toggle from spline
  multicanonical = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["Multicanonical"],false);
  MACH3LOG_INFO("Multicanonical Method: {}", multicanonical);

  multicanonicalSpline = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["MulticanonicalSpline"],false);

  if (multicanonical) {

    if (multicanonicalSpline){
      std::string splineFile = GetFromManager<std::string>(fitMan->raw()["General"]["MCMC"]["MulticanonicalSplineFile"],"nofile");
      TFile *file = new TFile(splineFile.c_str(), "READ");
      dcp_spline_IO = static_cast<TSpline3*>(file->Get("dcp_spline_IO"));
      MACH3LOG_INFO("Using multicanonical spline from file {}", splineFile);
      dcp_spline_IO->Eval(0.0); // check that the spline is valid
      std::cout << "Spline evaluated at 0.0 gives value " << dcp_spline_IO->Eval(0.0) << std::endl;

      dcp_spline_NO = static_cast<TSpline3*>(file->Get("dcp_spline_NO"));
      MACH3LOG_INFO("Using multicanonical spline from file {}", splineFile);
      dcp_spline_NO->Eval(0.0); // check that the spline is valid
      std::cout << "Spline evaluated at 0.0 gives value " << dcp_spline_NO->Eval(0.0) << std::endl;

      // // check for empty TSpline
      // if (dcp_spline == nullptr) {
      //   MACH3LOG_ERROR("Spline not found in file {}", splineFile);
      //   throw std::runtime_error("Spline not found in file");
      // }

    } else {
      // Get the multicanonical sigma values from the configuration file
      multicanonicalSigma = fitMan->raw()["General"]["MCMC"]["MulticanonicalSigma"].as<double>();
    }

    // Get the multicanonical beta value from the configuration file
    multicanonicalBeta = fitMan->raw()["General"]["MCMC"]["MulticanonicalBeta"].as<double>();
    MACH3LOG_INFO("Setting multicanonical beta to {}", multicanonicalBeta);


    MACH3LOG_INFO("Looping over systematics to find delta_cp parameter");
    // Loop over the systematics and find the osc_cov systematic and the delta_cp parameter number
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
          if (systematics[static_cast<int>(s)]->GetParName(i) == "delm2_23") {
            multicanonicalVar_dm23 = i;
            MACH3LOG_INFO("Setting delm2_23 parameter int {}",i);
          } else {
            MACH3LOG_ERROR("No dm23 parameter found in osc_cov systematic");
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

  //Save the Adaptive output
  for(const auto& syst : systematics){
    if(syst->getUseAdaptive()){
      syst->getAdaptiveHandler()->SaveAdaptiveToFile(syst->getAdaptiveHandler()->output_file_name, syst->getName(), true);
    }
  }

  // Process MCMC
  ProcessMCMC();
}
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

    // get the proposed value of delta_cp and apply the multicanonical pentalty, weighting it using the beta value to increase or decrease the strenght of the penalty
    delta_cp_value = systematics[oscCovVar]->getParProp(multicanonicalVar);
    delm23_value = systematics[oscCovVar]->getParProp(multicanonicalVar_dm23);
    
    if (multicanonicalSpline) {
      // Get the multicanonical weight from the spline
      multicanonical_penalty = GetMulticanonicalWeightSpline(delta_cp_value, dcp_spline_IO, dcp_spline_NO, delm23_value)*(multicanonicalBeta);
    } else {
      // Get the multicanonical weight from the Gaussian
      multicanonical_penalty = GetMulticanonicalWeightGaussian(delta_cp_value)*(multicanonicalBeta);
    }
    llh += multicanonical_penalty;
    
    // MACH3LOG_INFO("Delta CP value: {}", delta_cp_value);
    // MACH3LOG_INFO("Multicanonical penalty: {}", multicanonical_penalty);
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

inline double mcmc::GetMulticanonicalWeightGaussian(double deltacp){
  // precalculate constants
  constexpr double inv_sqrt_2pi = 0.3989422804014337;
  double sigma = multicanonicalSigma;
  double neg_half_sigma_sq = -1/(2*sigma*sigma); // sigma = 1 => -0.5; sigma = 0.5 => -2

  // three gaussians centered at -pi, 0, pi with sigma pre-defined above
  double exp1 = std::exp(neg_half_sigma_sq * (deltacp - TMath::Pi()) * (deltacp - TMath::Pi()));
  double exp2 = std::exp(neg_half_sigma_sq * (deltacp) * (deltacp));
  double exp3 = std::exp(neg_half_sigma_sq * (deltacp + TMath::Pi()) * (deltacp + TMath::Pi()));
  ///delta_cp_log_likelihood = -TMath::Log(TMath::Gaus(deltacp,TMath::Pi(),1,kTRUE)+TMath::Gaus(deltacp,0,1,kTRUE)+TMath::Gaus(deltacp,-TMath::Pi(),1,kTRUE));

  // return the log likelihood, ie the log of the normalised sum of the gaussians
  return -std::log(inv_sqrt_2pi * (1/sigma) * (exp1 + exp2 + exp3));
}

inline double mcmc::GetMulticanonicalWeightSpline(double deltacp, TSpline3 *spline_IO, TSpline3 *spline_NO, double delm23){
  double dcp_spline_val;

  if (delm23 < 0){
    dcp_spline_val = spline_IO->Eval(deltacp);
    return -2*(-std::log(dcp_spline_val)+std::log(spline_IO->Eval(-TMath::Pi()/2)));
  } else {
    dcp_spline_val = spline_NO->Eval(deltacp);
    return -2*(-std::log(dcp_spline_val)+std::log(spline_NO->Eval(-TMath::Pi()/2)));
  }
  // std::cout << "Evaluating spline at delta_cp = " << deltacp << " gives value " << dcp_spline_val << "with -log lh of :" << -log(dcp_spline_val) << std::endl;
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
