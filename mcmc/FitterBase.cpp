#include "FitterBase.h"


// *************************
// Initialise the manager and make it an object of FitterBase class
// Now we can dump manager settings to the output file
FitterBase::FitterBase(manager * const man) : fitMan(man) {
// *************************
  
  std::cout << "Trying to get information from fitMan" << std::endl;
  std::cout << "Seed is " << fitMan->raw()["General"]["Seed"].as<int>() << std::endl;

  random = new TRandom3(fitMan->raw()["General"]["Seed"].as<int>());

  // Counter of the accepted # of steps
  accCount = 0;
  step = 0;

  clock = new TStopwatch;
  stepClock = new TStopwatch;
  // Fit summary and debug info
  debug = fitMan->raw()["General"]["Debug"].as<bool>();
  std::string outfile = fitMan->raw()["General"]["Output"]["FileName"].as<std::string>();

  // Could set these from manager which is passed!
  osc_only = false;

  // Save output every auto_save steps
  //you don't want this too often https://root.cern/root/html606/TTree_8cxx_source.html#l01229
  auto_save = fitMan->raw()["General"]["MCMC"]["AutoSave"].as<int>();
  // Do we want to save the nominal parameters to output
  save_nominal = true;

  // Set the output file
  outputFile = new TFile(outfile.c_str(), "RECREATE");
  outputFile->cd();
  // Set output tree
  outTree = new TTree("posteriors", "Posterior_Distributions");
  // Auto-save every 200MB, the bigger the better https://root.cern/root/html606/TTree_8cxx_source.html#l01229
  outTree->SetAutoSave(-200E6);

  // Prepare the output log file
  if (debug) debugFile.open((outfile+".log").c_str());

  // Clear the samples and systematics
  samples.clear();
  systematics.clear();
  osc = NULL;

  sample_llh = NULL;
  syst_llh = NULL;
}


// *************************
// Destructor: close the logger and output file
FitterBase::~FitterBase() {
// *************************
  delete random;
  delete[] sample_llh;
  delete[] syst_llh;
  delete clock;
  delete stepClock;
  std::cout << "Done!" << std::endl;
}

// *******************
// Prepare the output tree
void FitterBase::PrepareOutput() {
// *******************

  // Check that we have added samples
  if (!samples.size()) {
    std::cerr << "No samples! Stopping MCMC" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Prepare the output trees
  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it) {
    (*it)->setBranches(*outTree);
  }

  if (osc) {
    osc->setBranches(*outTree);
    outTree->Branch("LogL_osc", &osc_llh, "LogL_osc/D");
  }

  outTree->Branch("LogL", &logLCurr, "LogL/D");
  outTree->Branch("accProb", &accProb, "accProb/D");
  outTree->Branch("step", &step, "step/I");
  outTree->Branch("stepTime", &stepTime, "stepTime/D");

  // Store individual likelihood components
  // Using double means double as large output file!
  sample_llh = new double[samples.size()];
  syst_llh = new double[systematics.size()];

  for (size_t i = 0; i < samples.size(); ++i) {
    std::stringstream oss, oss2;
    oss << "LogL_sample_" << i;
    oss2 << oss.str() << "/D";
    outTree->Branch(oss.str().c_str(), &sample_llh[i], oss2.str().c_str());

    // For adding sample dependent branches to the posteriors tree
    samples[i]->setMCMCBranches(outTree);
  }

  for (size_t i = 0; i < systematics.size(); ++i) {
    std::stringstream oss, oss2;
    oss << "LogL_systematic_" << systematics[i]->getName();
    oss2 << oss.str() << "/D";
    outTree->Branch(oss.str().c_str(), &syst_llh[i], oss2.str().c_str());
  }

  std::cout << "\n-------------------- Starting MCMC --------------------" << std::endl;

  if (debug) {
    PrintInitialState();
    debugFile << "----- Starting MCMC -----" << std::endl;
  }
  // Time the progress
  clock->Start();
}

// *******************
void FitterBase::SaveOutput() {
// *******************

  //Stop Clock
  clock->Stop();

  outputFile->cd();
  outTree->Write();

  std::cout << "\n" << step << " steps took " << clock->RealTime() << " seconds to complete. (" << clock->RealTime() / step << "s / step).\n" << accCount<< " steps were accepted." << std::endl;

  if (debug)
  {
      debugFile << "\n\n" << step << " steps took " << clock->RealTime() << " seconds to complete. (" << clock->RealTime() / step << "s / step).\n" << accCount<< " steps were accepted." << std::endl;
  }

  if(debug) debugFile.close();

  outputFile->Close();
}



// *************************
// Add samplePDF object to the Markov Chain
void FitterBase::addSamplePDF(samplePDFBase * const sample) {
// *************************
  std::cout << "Adding samplePDF object " << std::endl;
  samples.push_back(sample);

  return;
}


// *************************
// Add flux systematics, cross-section systematics, ND systematics to the chain
void FitterBase::addSystObj(covarianceBase * const cov) {
// *************************

  std::cout << "Adding systematic object " << cov->getName() << std::endl;
  systematics.push_back(cov);

  // Save an array of nominal
  if (save_nominal) {
    std::vector<double> vec = cov->getNominalArray();
    size_t n = vec.size();
    double *n_vec = new double[n];
    for (size_t i = 0; i < n; ++i) {
      n_vec[i] = vec[i];
    }
    TVectorT<double> t_vec(n, n_vec);
    TString nameof = TString(cov->getName());
    nameof = nameof.Append("_nom");
    t_vec.Write(nameof);
    delete[] n_vec;
  }

  return;
}

// *************************
// Add an oscillation handler
// Similar to systematic really, but handles the oscillation weights
void FitterBase::addOscHandler(covarianceOsc * const oscf) {
// *************************

  osc = oscf;

  if (save_nominal) {

    std::vector<double> vec = oscf->getNominalArray();
    size_t n = vec.size();
    double *n_vec = new double[n];
    for (size_t i = 0; i < n; ++i) {
      n_vec[i] = vec[i];
    }
    TVectorT<double> t_vec(n, n_vec);
    TString nameof = TString(oscf->getName());
    nameof = nameof.Append("_nom");
    t_vec.Write(nameof);
    delete[] n_vec;
  }

  return;
}

// *************************
// When using separate oscillations for neutrino and anti-neutrino
void FitterBase::addOscHandler(covarianceOsc *oscf, covarianceOsc *oscf2) {
// *************************

  osc = oscf;
  osc2 = oscf2;

  if (save_nominal) {
    std::vector<double> vec = oscf->getNominalArray();
    size_t n = vec.size();
    double *n_vec = new double[n];
    for (size_t i = 0; i < n; ++i) {
      n_vec[i] = vec[i];
    }
    TVectorT<double> t_vec(n, n_vec);
    TString nameof = TString(oscf->getName());
    nameof = nameof.Append("_nom");
    t_vec.Write(nameof);

    std::vector<double> vec2 = oscf2->getNominalArray();
    size_t n2 = vec2.size();
    double *n_vec2 = new double[n];
    for (size_t i = 0; i < n; ++i) {
      n_vec2[i] = vec2[i];
    }
    TVectorT<double> t_vec2(n2, n_vec2);
    TString nameof2 = TString(oscf2->getName());
    nameof2 = nameof2.Append("_2_nom");
    t_vec2.Write(nameof2);
    delete[] n_vec;
    delete[] n_vec2;
  }

  // Set whether osc and osc2 should be equal (default: no for all parameters)
  for (int i=0; i<osc->getSize(); i++) {
    equalOscPar.push_back(0);
  }

  // Set whether to force osc and osc2 to use the same mass hierarchy (default: no)
  equalMH = false;

  return;
}


// **********************
// Print our initial state for the chain
void FitterBase::PrintInitialState() {
// **********************
  if (debug) {
    for (size_t i = 0; i < systematics.size(); ++i) {
      debugFile << "\nnominal values for " << systematics[i]->getName() << std::endl;
      std::vector<double> noms = systematics[i]->getNominalArray();
      for (size_t j = 0; j < noms.size(); ++j) {
        debugFile << noms[j] << " ";
      }
    }
    debugFile << std::endl;
  }
}
