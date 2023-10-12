#include "mcmc.h"

// *************************
// Initialise the manager and make it an object of mcmc class
// Now we can dump manager settings to the output file
mcmc::mcmc(manager *man) : FitterBase(man) {
// *************************

  // Beginning step number
  stepStart = 0;

  // Starting parameters should be thrown
  init_pos = false;
  reject = false;
  chainLength = fitMan->raw()["General.NSteps"].as<double>();

  //ETA - currently don't have this in manager as it needs some love
  //AnnealTemp = fitMan->GetTemp();
  AnnealTemp = -999;
  if(AnnealTemp < 0) anneal = false;
  else
  {
    std::cout << "- Enabling simulated annealing with T = " << AnnealTemp << std::endl;
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
  std::cout << "MCMC getting starting position from " << file << std::endl;

  TFile *infile = new TFile(file.c_str(), "READ");
  TTree *posts = (TTree*)infile->Get("posteriors");
  TObjArray* brlis = (TObjArray*)posts->GetListOfBranches();
  int nbr = brlis->GetEntries();
  TString branch_names[nbr];
  double branch_vals[nbr];

  for (int i = 0; i < nbr; ++i) {
    TBranch *br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();
    branch_names[i] = bname;
    std::cout << " * Loading " << bname << std::endl;
    posts->SetBranchAddress(branch_names[i], &branch_vals[i]);
  }

  posts->GetEntry(posts->GetEntries()-1);
  for (int i = 0; i < nbr; ++i) {
    init_pars.insert( std::pair<TString, double>(branch_names[i], branch_vals[i]));
  }

  init_pos = true;
  infile->Close();
  delete infile;
}

// *************************
// When we've passed a previous MCMC chain, now find the starting value
double mcmc::FindStartingValue(std::string par_name) {
// *************************

  if (!init_pos) {
    std::cout << "- Can't find starting value, file not specified!" << std::endl;
    throw;
  }

  std::map<TString, double>::const_iterator it = init_pars.find(par_name);

  if (it == init_pars.end()) {
    std::cout << "- Can't find parameter of name " << par_name << std::endl;
    throw;
  }

  return it->second;
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

  if (debug) debugFile << " logLProp: " << logLProp << " logLCurr: " << logLCurr << " accProb: " << accProb << " fRandom: " << fRandom << std::endl;

  // Update all the handlers to accept the step
  if (accept && !reject) {
    logLCurr = logLProp;

    if (osc) {
      osc->acceptStep();
    }
    if (osc2) {
      osc2->acceptStep();
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
  for (step = stepStart; step < stepStart+chainLength; step++) {

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
// Process the MCMC output to get postfit etc
void mcmc::ProcessMCMC() {
// *******************

  if (fitMan == NULL) return;

  // Process the MCMC
  if (fitMan->raw()["General.ProcessMCMC"].as<bool>()) {

    // Make the processor
    MCMCProcessor Processor(std::string(outputFile->GetName()), false);
    
    Processor.Initialise();
    // Make the TVectorD pointers which hold the processed output
    TVectorD *Central = NULL;
    TVectorD *Errors = NULL;
    TVectorD *Central_Gauss = NULL;
    TVectorD *Errors_Gauss = NULL;
    TVectorD *Peaks = NULL;
    
    // Make the postfit
    Processor.GetPostfit(Central, Errors, Central_Gauss, Errors_Gauss, Peaks);
    Processor.DrawPostfit();

    // Make the TMatrix pointers which hold the processed output
    TMatrixDSym *Covariance = NULL;
    TMatrixDSym *Correlation = NULL;

    // Make the covariance matrix
    Processor.GetCovariance(Covariance, Correlation);
    Processor.DrawCovariance();

    std::vector<TString> BranchNames = Processor.GetBranchNames();

    
    // Re-open the TFile
    if (!outputFile->IsOpen()) {
      std::cout << "Opening output again to update with means..." << std::endl;
      outputFile = new TFile(fitMan->raw()["General.Output.Filename"].as<std::string>().c_str(), "UPDATE");
    }

    
    Central->Write("PDF_Means");
    Errors->Write("PDF_Errors");
    Central_Gauss->Write("Gauss_Means");
    Errors_Gauss->Write("Errors_Gauss");
    Covariance->Write("Covariance");
    Correlation->Write("Correlation");
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

  // Propose steps for the oscillation handlers
  if (osc) {
    osc->proposeStep();

    // If we have separate neutrino and anti-neutrino oscillation parameters
    if (osc2) {
      osc2->proposeStep();

      // Fix some parameters in osc equal to osc2 (if it's set in equalOscPar)

      // double* returned by getPropPars() has oscpar[0-5] as you would expect (theta12, 
      // theta23, theta13, dm21, dm32, dcp), oscpar[6] = 2, oscpar[7] = L, 
      // and oscpar[8] = density. If beta is present it has oscpar[9] = beta. 
      // We don't want this - the vector<double> given to setParameters() should have
      // only the oscillation parametes (so size = 6 if no beta, or 7 if using beta).
      // Construct the correct vector in the loop (at the same time as setting some
      // parameters equal between osc and osc2 if desired)

      double* osc_proppars = osc->getPropPars();
      double* osc2_proppars = osc2->getPropPars();
      std::vector<double> oscpars;
      std::vector<double> oscpars2;

      for (int par = 0; par < int(osc->getSize()+3); par++) {
        // skip the extra entries
        if (par==6 || par==7 || par==8) {
          continue;

          // normal osc pars
        } else if (par<6) {
          // Set correct parameter values for osc = osc2
          if (equalOscPar[par]==1) {
            osc_proppars[par] = osc2_proppars[par];
          }

          // Set same mass hierarchy if desired
          // dm2
          if (equalMH && par==4) {
            double sign = osc_proppars[par]*osc2_proppars[par]; // +ve if both are same MH, -ve if not
            if (sign<0) {
              osc_proppars[par]*=(-1.0);
            }
          }

          // Now make vectors
          oscpars.push_back(osc_proppars[par]);
          oscpars2.push_back(osc2_proppars[par]);
          // beta
        } else if (par==9) {

          // Set correct parameter values for osc = osc2
          if (equalOscPar[par]==1) {
            osc_proppars[par] = osc2_proppars[par];
          }

          // Now make vectors
          oscpars.push_back(osc_proppars[par]);
          oscpars2.push_back(osc2_proppars[par]);
        }
      }

      // Now set the correct parameters
      osc->setParameters(oscpars);
      osc2->setParameters(oscpars2);
    }

    // Now get the likelihoods for the oscillation
    osc_llh = osc->getLikelihood();
    if (osc2) {
      osc_llh += osc2->getLikelihood();
    }

    // Add the oscillation likelihoods to the reconfigure likelihoods
    llh += osc_llh;

    if (debug) debugFile << "LLH for oscillation handler: " << llh << std::endl;
  }

  int stdIt = 0;
  // Loop over the systematics and propose the initial step
  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it, ++stdIt) {

    // Could throw the initial value here to do MCMC stability studies
    // Propose the steps for the systematics
    if (!osc_only) {
      (*it)->proposeStep();
    }

    // Get the likelihood from the systematics
    syst_llh[stdIt] = (*it)->getLikelihood();
    llh += syst_llh[stdIt];

    if (debug) debugFile << "LLH after " << systematics[stdIt]->getName() << " " << llh << std::endl;
  }

  // Check if we've hit a boundary in the systematics
  // In this case we can save time by not having to reconfigure the simulation
  if (llh >= __LARGE_LOGL__) {
    reject = true;
    if (debug) debugFile << "Rejecting based on boundary" << std::endl;
  }

  // Only reweight when we have a good parameter configuration
  // This speeds things up considerably because for every bad parameter configuration we don't have to reweight the MC
  if (!reject) {

    // Could multi-thread this
    // But since sample reweight is multi-threaded it's probably better to do that
    for (size_t i = 0; i < samples.size(); i++) {

      // If we're running with different oscillation parameters for neutrino and anti-neutrino
      if (osc && osc2) {
        samples[i]->reweight(osc->getPropPars(), osc2->getPropPars());
      } else if (osc) {
        samples[i]->reweight(osc->getPropPars());
        // If we aren't using any oscillation
      } else {
        double* fake = NULL;
        samples[i]->reweight(fake);
      }
    }

    //DB for atmospheric event by event sample migration, need to fully reweight all samples to allow event passing prior to likelihood evaluation
    for (size_t i = 0; i < samples.size(); i++) {
      // Get the sample likelihoods and add them
      sample_llh[i] = samples[i]->getLikelihood();
      llh += sample_llh[i];
      if (debug) debugFile << "LLH after sample " << i << " " << llh << std::endl;
    }

  // For when we don't have to reweight, set sample to madness
  } else {
    for (size_t i = 0; i < samples.size(); i++) {
      // Set the sample_llh[i] to be madly high also to signfiy a step out of bounds
      sample_llh[i] = __LARGE_LOGL__;
      if (debug) debugFile << "LLH after REJECT sample " << i << " " << llh << std::endl;
    }
  }

  // Save the proposed likelihood (class member)
  logLProp = llh;
}

// *******************
// Print the fit output progress
void mcmc::PrintProgress() {
// *******************

  std::cout << "Step:\t" << step-stepStart << "/" << chainLength << "  |  current: " << logLCurr << " proposed: " << logLProp << std::endl;
  std::cout << "Accepted/Total steps: " << accCount << "/" << step-stepStart << " = " << double(accCount)/double(step - stepStart) << std::endl;

  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it) {
    if (std::string((*it)->getName()) == "xsec_cov") {
      std::cout << "Cross-section parameters: " << std::endl;
      (*it)->printNominalCurrProp();
    }
  }

  if (debug) {
    debugFile << "\n-------------------------------------------------------" << std::endl;
    debugFile << "Step:\t" << step + 1 << "/" << chainLength << "  |  current: " << logLCurr << " proposed: " << logLProp << std::endl;
  }
}
