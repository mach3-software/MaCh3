#include "FitterBase.h"

#include "TRandom.h"
#include "TStopwatch.h"

// *************************
// Initialise the manager and make it an object of FitterBase class
// Now we can dump manager settings to the output file
FitterBase::FitterBase(manager * const man) : fitMan(man) {
// *************************

  random = new TRandom3(fitMan->raw()["General"]["Seed"].as<int>());

  // Counter of the accepted # of steps
  accCount = 0;
  step = 0;

  clock = new TStopwatch;
  stepClock = new TStopwatch;
  #ifdef DEBUG
  // Fit summary and debug info
  debug = fitMan->raw()["General"]["Debug"].as<bool>();
  #endif

  std::string outfile = fitMan->raw()["General"]["OutputFile"].as<std::string>();

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

  FileSaved = false;
  SettingsSaved = false;
  OutputPrepared = false;

  //Create TDirectory
  CovFolder = outputFile->mkdir("CovarianceFolder");
  outputFile->cd();

  #ifdef DEBUG
  // Prepare the output log file
  if (debug) debugFile.open((outfile+".log").c_str());
  #endif

  // Clear the samples and systematics
  samples.clear();
  systematics.clear();
  osc = NULL;

  sample_llh = nullptr;
  syst_llh = nullptr;

  fTestLikelihood = false;
  //ETA - No guarantee that "Fitter" field exists so check this first before
  //checking ["Fitter"]["FitTestLikelihood"]
  if(fitMan->raw()["General"]["Fitter"])
  {
    if(fitMan->raw()["General"]["Fitter"]["FitTestLikelihood"]){
      fTestLikelihood = fitMan->raw()["General"]["Fitter"]["FitTestLikelihood"].as<bool>();
    }
  }

}


// *************************
// Destructor: close the logger and output file
FitterBase::~FitterBase() {
// *************************
  SaveOutput();
  if(random != nullptr) delete random;
  if(sample_llh != nullptr) delete[] sample_llh;
  if(syst_llh != nullptr) delete[] syst_llh;
  if(outputFile != nullptr) delete outputFile;
  delete clock;
  delete stepClock;
  MACH3LOG_INFO("Closing MaCh3 Fitter Engine");
}


// *******************
// Prepare the output tree
void FitterBase::SaveSettings() {
// *******************

  if(SettingsSaved) return;
  fitMan->SaveSettings(outputFile);
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  MACH3LOG_INFO("#####Current Setup#####");
  MACH3LOG_INFO("Number of covariances: {}", systematics.size());
  for(unsigned int i = 0; i < systematics.size(); i++)
    MACH3LOG_INFO("{}: Cov name: {}",i,  systematics[i]->getName());
  MACH3LOG_INFO("Number of SamplePDFs: {}", samples.size());
  for(unsigned int i = 0; i < systematics.size(); i++)
    MACH3LOG_INFO("{}: SamplePDF name: {}",i , samples[i]->GetName());

  SettingsSaved = true;
}

// *******************
// Prepare the output tree
void FitterBase::PrepareOutput() {
// *******************

  if(OutputPrepared) return;
  //MS: Check if we are fitting the test likelihood, rather than T2K likelihood, and only setup T2K output if not
  if(!fTestLikelihood)
  {
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
  }
  else
  {
    outTree->Branch("LogL", &logLCurr, "LogL/D");
    outTree->Branch("accProb", &accProb, "accProb/D");
    outTree->Branch("step", &step, "step/I");
    outTree->Branch("stepTime", &stepTime, "stepTime/D");
  }

  MACH3LOG_INFO("-------------------- Starting MCMC --------------------");
  #ifdef DEBUG
  if (debug) {
    PrintInitialState();
    debugFile << "----- Starting MCMC -----" << std::endl;
  }
  #endif
  // Time the progress
  clock->Start();

  OutputPrepared = true;
}

// *******************
void FitterBase::SaveOutput() {
// *******************

  if(FileSaved) return;
  //Stop Clock
  clock->Stop();

  outputFile->cd();
  outTree->Write();

  MACH3LOG_INFO("{} steps took {:.2f} seconds to complete. ({:.2f}s / step).", step, clock->RealTime(), clock->RealTime() / step);
  MACH3LOG_INFO("{} steps were accepted.", accCount);
  #ifdef DEBUG
  if (debug)
  {
    debugFile << "\n\n" << step << " steps took " << clock->RealTime() << " seconds to complete. (" << clock->RealTime() / step << "s / step).\n" << accCount<< " steps were accepted." << std::endl;
    debugFile.close();
  }
  #endif

  outputFile->Close();
  FileSaved = true;
}

// *************************
// Add samplePDF object to the Markov Chain
void FitterBase::addSamplePDF(samplePDFBase * const sample) {
// *************************
  MACH3LOG_INFO("Adding {} object ", sample->GetName());
  samples.push_back(sample);

  return;
}

// *************************
// Add flux systematics, cross-section systematics, ND systematics to the chain
void FitterBase::addSystObj(covarianceBase * const cov) {
// *************************

  MACH3LOG_INFO("Adding systematic object {}", cov->getName());
  systematics.push_back(cov);

  // Save an array of nominal
  if (save_nominal) {
    CovFolder->cd();
    double *n_vec = new double[cov->getSize()];
    for (int i = 0; i < cov->getSize(); ++i)
      n_vec[i] = cov->getParInit(i);

    TVectorT<double> t_vec(cov->getSize(), n_vec);
    t_vec.Write((std::string(cov->getName()) + "_prior").c_str());
    delete[] n_vec;

    cov->getCovMatrix()->Write(cov->getName());
    outputFile->cd();
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
    CovFolder->cd();
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
    outputFile->cd();
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
    CovFolder->cd();
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

    outputFile->cd();
  }

  return;
}


// **********************
// Print our initial state for the chain
void FitterBase::PrintInitialState() {
// **********************

  #ifdef DEBUG
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
  #endif
}

// *************************
// Run LLH scan
void FitterBase::RunLLHScan() {
// *************************

  // Save the settings into the output file
  SaveSettings();

  MACH3LOG_INFO("Starting LLH Scan");
  int TotalNSamples = 0;
  for(unsigned int i = 0; i < samples.size(); i++ )
  {
    TotalNSamples += samples[i]->GetNsamples();
  }

  //KS: Turn it on if you want LLH scan for each ND sample separately, which increase time significantly but can be useful for validating new samples or dials.
  bool PlotAllNDsamplesLLH = false;
  if(fitMan->raw()["General"]["LLHScanBySample"])
    PlotAllNDsamplesLLH = fitMan->raw()["General"]["LLHScanBySample"].as<bool>();

  std::vector<std::string> SkipVector;
  if(fitMan->raw()["General"]["LLHScanSkipVector"])
  {
    SkipVector = fitMan->raw()["General"]["LLHScanSkipVector"].as<std::vector<std::string>>();
    std::cout<<" Found skip vector with "<<SkipVector.size()<<" entries "<<std::endl;
  }

  // Now finally get onto the LLH scan stuff
  // Very similar code to MCMC but never start MCMC; just scan over the parameter space

  std::vector<TDirectory *> Cov_LLH(systematics.size());
  for(unsigned int ivc = 0; ivc < systematics.size(); ivc++ )
  {
    std::string NameTemp = systematics[ivc]->getName();
    NameTemp = NameTemp.substr(0, NameTemp.find("_cov")) + "_LLH";
    Cov_LLH[ivc] = outputFile->mkdir(NameTemp.c_str());
  }

  std::vector<TDirectory *> SampleClass_LLH(samples.size());
  for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
  {
    std::string NameTemp = samples[ivs]->GetName();
    SampleClass_LLH[ivs] = outputFile->mkdir(NameTemp.c_str());
  }

  TDirectory *Sample_LLH = outputFile->mkdir("Sample_LLH");
  TDirectory *Total_LLH = outputFile->mkdir("Total_LLH");

  TDirectory **SampleSplit_LLH = nullptr;
  if(PlotAllNDsamplesLLH)
  {
    SampleSplit_LLH = new TDirectory*[TotalNSamples];
    int iterator = 0;
    for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
    {
      for(__int__ is = 0; is < samples[ivs]->GetNsamples(); is++ )
      {
        SampleSplit_LLH[iterator] = outputFile->mkdir((samples[ivs]->GetSampleName(is)+ "_LLH").c_str());
        iterator++;
      }
    }
  }
  // Number of points we do for each LLH scan
  const int n_points = GetFromManager<int>(fitMan->raw()["General"]["LLHScanPoints"], 100);

  // We print 5 reweights
  const int countwidth = double(n_points)/double(5);

  bool isxsec = false;
  // Loop over the covariance classes
  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    if (std::string((*it)->getName()) == "xsec_cov")
    {
      isxsec = true;
    } else {
      isxsec = false;
    }

    // Scan over all the parameters
    // Get the number of parameters
    int npars = (*it)->getSize();
    bool IsPCA = (*it)->IsPCA();
    if (IsPCA) npars = (*it)->getNpars();
    for (int i = 0; i < npars; ++i)
    {
      // Get the parameter name
      std::string name = (*it)->GetParName(i);
      if (IsPCA) name += "_PCA";
      // For xsec we can get the actual name, hurray for being informative
      if (isxsec) name = (*it)->GetParFancyName(i);
      bool skip = false;
      for(unsigned int is = 0; is < SkipVector.size(); is++)
      {
        if(name.substr(0, SkipVector[is].length()) == SkipVector[is])
        {
          skip = true;
          break;
        }
      }
      if(skip) continue;

      // Get the parameter priors and bounds
      double prior = (*it)->getParInit(i);
      if (IsPCA) prior = (*it)->getParCurr_PCA(i);

      // Get the covariance matrix and do the +/- nSigma
      double nSigma = 1;
      if (IsPCA) nSigma = 0.5;
      // Set lower and upper bounds relative the prior
      double lower = prior - nSigma*(*it)->getDiagonalError(i);
      double upper = prior + nSigma*(*it)->getDiagonalError(i);
      // If PCA, transform these parameter values to the PCA basis
      if (IsPCA) {
        lower = prior - nSigma*sqrt(((*it)->getEigenValues())(i));
        upper = prior + nSigma*sqrt(((*it)->getEigenValues())(i));

        std::cout << "eval " << i << " = " << (*it)->getEigenValues()(i) << std::endl;
        std::cout << "prior " << i << " = " << prior << std::endl;
        std::cout << "lower " << i << " = " << lower << std::endl;
        std::cout << "upper " << i << " = " << upper << std::endl;
        std::cout << "nSigma " << nSigma << std::endl;
      }

      // Cross-section and flux parameters have boundaries that we scan between, check that these are respected in setting lower and upper variables
      if (lower < (*it)->GetLowerBound(i)) {
        lower = (*it)->GetLowerBound(i);
      }
      if (upper > (*it)->GetUpperBound(i)) {
        upper = (*it)->GetUpperBound(i);
      }
      MACH3LOG_INFO("Scanning {} with {} steps, from {} - {}, prior = {}", name, n_points, lower, upper, prior);

      // Make the TH1D
      TH1D *hScan = new TH1D((name+"_full").c_str(), (name+"_full").c_str(), n_points, lower, upper);
      hScan->SetTitle(std::string(std::string("2LLH_full, ") + name + ";" + name + "; -2(ln L_{sample} + ln L_{xsec+flux} + ln L_{det})").c_str());

      TH1D *hScanSam = new TH1D((name+"_sam").c_str(), (name+"_sam").c_str(), n_points, lower, upper);
      hScanSam->SetTitle(std::string(std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());


      TH1D **hScanSample = new TH1D*[samples.size()];
      double *nSamLLH = new double[samples.size()];
      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        std::string NameTemp = samples[ivs]->GetName();
        hScanSample[ivs] = new TH1D((name+"_"+NameTemp).c_str(), (name+"_" + NameTemp).c_str(), n_points, lower, upper);
        hScanSample[ivs]->SetTitle(std::string(std::string("2LLH_" + NameTemp + ", ") + name + ";" + name + "; -2(ln L_{" + NameTemp +"})").c_str());
        nSamLLH[ivs] = 0.;
      }

      TH1D **hScanCov = new TH1D*[systematics.size()];
      double *nCovLLH = new double[systematics.size()];
      for(unsigned int ivc = 0; ivc < systematics.size(); ivc++ )
      {
        std::string NameTemp = systematics[ivc]->getName();
        NameTemp = NameTemp.substr(0, NameTemp.find("_cov"));

        hScanCov[ivc] = new TH1D((name+"_"+NameTemp).c_str(), (name+"_" + NameTemp).c_str(), n_points, lower, upper);
        hScanCov[ivc]->SetTitle(std::string(std::string("2LLH_" + NameTemp + ", ") + name + ";" + name + "; -2(ln L_{" + NameTemp +"})").c_str());
        nCovLLH[ivc] = 0.;
      }

      TH1D **hScanSamSplit = nullptr;
      double *sampleSplitllh = nullptr;
      if(PlotAllNDsamplesLLH)
      {
        int iterator = 0;
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          hScanSamSplit = new TH1D*[TotalNSamples];
          sampleSplitllh = new double[TotalNSamples];
          for(__int__ is = 0; is < samples[ivs]->GetNsamples(); is++ )
          {
            hScanSamSplit[iterator] = new TH1D( (name+samples[ivs]->GetSampleName(is)).c_str(), (name+samples[ivs]->GetSampleName(is)).c_str(), n_points, lower, upper );
            hScanSamSplit[iterator]->SetTitle(std::string(std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());
            iterator++;
          }
        }
      }

      // Scan over the parameter space
      for (int j = 0; j < n_points; j++)
      {
        if (j % countwidth == 0)
          MaCh3Utils::PrintProgressBar(j, n_points);

        // For PCA we have to do it differently
        if (IsPCA) {
          (*it)->setParProp_PCA(i, hScan->GetBinCenter(j+1));
        } else {
          // Set the parameter
          (*it)->setParProp(i, hScan->GetBinCenter(j+1));
        }

        // Reweight the MC
        double *fake = 0;
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          samples[ivs]->reweight(fake);
        }
        //Total LLH
        double totalllh = 0;

        // Get the -log L likelihoods
        double samplellh = 0;

        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          nSamLLH[ivs] = samples[ivs]->GetLikelihood();
          samplellh += nSamLLH[ivs];
        }

        for(unsigned int ivc = 0; ivc < systematics.size(); ivc++ )
        {
          nCovLLH[ivc] = systematics[ivc]->GetLikelihood();
          totalllh += nCovLLH[ivc];
        }

        totalllh += samplellh;

        if(PlotAllNDsamplesLLH)
        {
          int iterator = 0;
          for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
          {
            for(__int__ is = 0; is < samples[ivs]->GetNsamples(); is++)
            {
              sampleSplitllh[iterator] = samples[ivs]->getSampleLikelihood(is);
              iterator++;
            }
          }
        }

        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          hScanSample[ivs]->SetBinContent(j+1, 2*nSamLLH[ivs]);
        }
        for(unsigned int ivc = 0; ivc < systematics.size(); ivc++ )
        {
          hScanCov[ivc]->SetBinContent(j+1, 2*nCovLLH[ivc]);
        }

        hScanSam->SetBinContent(j+1, 2*samplellh);
        hScan->SetBinContent(j+1, 2*totalllh);

        if(PlotAllNDsamplesLLH)
        {
          int iterator = 0;
          for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
          {
            for(__int__ is = 0; is < samples[ivs]->GetNsamples(); is++)
            {
              hScanSamSplit[is]->SetBinContent(j+1, 2*sampleSplitllh[is]);
              iterator++;
            }
          }
        }
      }
      for(unsigned int ivc = 0; ivc < systematics.size(); ivc++ )
      {
        Cov_LLH[ivc]->cd();
        hScanCov[ivc]->Write();
        delete hScanCov[ivc];
      }

      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        SampleClass_LLH[ivs]->cd();
        hScanSample[ivs]->Write();
        delete hScanSample[ivs];
      }

      Sample_LLH->cd();
      hScanSam->Write();
      Total_LLH->cd();
      hScan->Write();

      delete[] hScanCov;
      delete[] nCovLLH;
      delete[] hScanSample;
      delete[] nSamLLH;
      delete hScanSam;
      delete hScan;

      hScanCov = nullptr;
      nCovLLH = nullptr;
      hScanSample = nullptr;
      nSamLLH = nullptr;
      hScanSam = nullptr;
      hScan = nullptr;

      if(PlotAllNDsamplesLLH)
      {
        int iterator = 0;
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          for(__int__ is = 0; is < samples[ivs]->GetNsamples(); is++)
          {
            SampleSplit_LLH[iterator]->cd();
            hScanSamSplit[iterator]->Write();
            delete hScanSamSplit[iterator];
            iterator++;
          }
        }
        delete[] hScanSamSplit;
      }

      // Reset the parameters to their prior central values
      if (IsPCA) {
        (*it)->setParProp_PCA(i, prior);
      } else {
        (*it)->setParProp(i, prior);
      }
    }//end loop over systematics
  }//end loop covariance classes

  for(unsigned int ivc = 0; ivc < systematics.size(); ivc++ )
  {
    Cov_LLH[ivc]->Write();
    delete Cov_LLH[ivc];
  }

  for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
  {
    SampleClass_LLH[ivs]->Write();
    delete SampleClass_LLH[ivs];
  }

  Sample_LLH->Write();
  delete Sample_LLH;

  Total_LLH->Write();
  delete Total_LLH;

  if(PlotAllNDsamplesLLH)
  {
    int iterator = 0;
    for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
    {
      for(__int__ is = 0; is < samples[ivs]->GetNsamples(); is++ )
      {
        SampleSplit_LLH[iterator]->Write();
        delete SampleSplit_LLH[iterator];
        iterator++;
      }
    }
  }
  SaveOutput();
}


// *************************
// Run 2D LLH scan
void FitterBase::Run2DLLHScan() {
// *************************

  // Save the settings into the output file
  SaveSettings();

  TDirectory *Sample_2DLLH = outputFile->mkdir("Sample_2DLLH");
  std::vector<std::string> SkipVector;
  if(fitMan->raw()["General"]["LLHScanSkipVector"])
  {
    SkipVector = fitMan->raw()["General"]["LLHScanSkipVector"].as<std::vector<std::string>>();
    std::cout<<" Found skip vector with "<<SkipVector.size()<<" entries "<<std::endl;
  }

  // Number of points we do for each LLH scan
  const int n_points = 20;
  // We print 5 reweights
  const int countwidth = double(n_points)/double(5);

  bool isxsec = false;
  // Loop over the covariance classes
  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    if (std::string((*it)->getName()) == "xsec_cov")
    {
      isxsec = true;
    } else {
      isxsec = false;
    }
    // Scan over all the parameters
    // Get the number of parameters
    int npars = (*it)->getSize();
    bool IsPCA = (*it)->IsPCA();
    if (IsPCA) npars = (*it)->getNpars();

    for (int i = 0; i < npars; ++i)
    {

      std::string name_x = (*it)->GetParName(i);
      if (IsPCA) name_x += "_PCA";
      // For xsec we can get the actual name, hurray for being informative
      if (isxsec) name_x = (*it)->GetParFancyName(i);

      // Get the parameter priors and bounds
      double prior_x = (*it)->getParInit(i);
      if (IsPCA) prior_x = (*it)->getParCurr_PCA(i);

      // Get the covariance matrix and do the +/- nSigma
      double nSigma = 1;
      if (IsPCA) nSigma = 0.5;
      // Set lower and upper bounds relative the prior
      double lower_x = prior_x - nSigma*(*it)->getDiagonalError(i);
      double upper_x = prior_x + nSigma*(*it)->getDiagonalError(i);
      // If PCA, transform these parameter values to the PCA basis
      if (IsPCA) {
        lower_x = prior_x - nSigma*sqrt(((*it)->getEigenValues())(i));
        upper_x = prior_x + nSigma*sqrt(((*it)->getEigenValues())(i));

        std::cout << "eval " << i << " = " << (*it)->getEigenValues()(i) << std::endl;
        std::cout << "prior " << i << " = " << prior_x << std::endl;
        std::cout << "lower " << i << " = " << lower_x << std::endl;
        std::cout << "upper " << i << " = " << upper_x << std::endl;
        std::cout << "nSigma " << nSigma << std::endl;
      }

      // Cross-section and flux parameters have boundaries that we scan between, check that these are respected in setting lower and upper variables
      if (lower_x < (*it)->GetLowerBound(i)) {
        lower_x = (*it)->GetLowerBound(i);
      }
      if (upper_x > (*it)->GetUpperBound(i)) {
        upper_x = (*it)->GetUpperBound(i);
      }

      bool skip = false;
      for(unsigned int is = 0; is < SkipVector.size(); is++)
      {
        if(name_x.substr(0, SkipVector[is].length()) == SkipVector[is])
        {
          skip = true;
          break;
        }
      }
      if(skip) continue;

      for (int j = 0; j < i; ++j)
      {
        std::string name_y = (*it)->GetParName(j);
        if (IsPCA) name_y += "_PCA";
        // For xsec we can get the actual name, hurray for being informative
        if (isxsec) name_y = (*it)->GetParFancyName(j);

        bool skip = false;
        for(unsigned int is = 0; is < SkipVector.size(); is++)
        {
          if(name_y.substr(0, SkipVector[is].length()) == SkipVector[is])
          {
            skip = true;
            break;
          }
        }
        if(skip) continue;

        // Get the parameter priors and bounds
        double prior_y = (*it)->getParInit(j);
        if (IsPCA) prior_y = (*it)->getParCurr_PCA(j);

        // Set lower and upper bounds relative the prior
        double lower_y = prior_y - nSigma*(*it)->getDiagonalError(j);
        double upper_y = prior_y + nSigma*(*it)->getDiagonalError(j);
        // If PCA, transform these parameter values to the PCA basis
        if (IsPCA) {
          lower_y = prior_y - nSigma*sqrt(((*it)->getEigenValues())(j));
          upper_y = prior_y + nSigma*sqrt(((*it)->getEigenValues())(j));

          std::cout << "eval " << j << " = " << (*it)->getEigenValues()(j) << std::endl;
          std::cout << "prior " << j << " = " << prior_y << std::endl;
          std::cout << "lower " << j << " = " << lower_y << std::endl;
          std::cout << "upper " << j << " = " << upper_y << std::endl;
          std::cout << "nSigma " << nSigma << std::endl;
        }

        // Cross-section and flux parameters have boundaries that we scan between, check that these are respected in setting lower and upper variables
        if (lower_y < (*it)->GetLowerBound(j)) {
          lower_y = (*it)->GetLowerBound(j);
        }
        if (upper_y > (*it)->GetUpperBound(j)) {
          upper_y = (*it)->GetUpperBound(j);
        }
        MACH3LOG_INFO("Scanning X {} with {} steps, from {} - {}, prior = {}", name_x, n_points, lower_x, upper_x, prior_x);
        MACH3LOG_INFO("Scanning Y {} with {} steps, from {} - {}, prior = {}", name_y, n_points, lower_y, upper_y, prior_y);

        TH2D *hScanSam = new TH2D((name_x + "_" + name_y + "_sam").c_str(), (name_x + "_" + name_y + "_sam").c_str(), n_points, lower_x, upper_x, n_points, lower_y, upper_y);
        hScanSam->GetXaxis()->SetTitle(name_x.c_str());
        hScanSam->GetYaxis()->SetTitle(name_y.c_str());
        hScanSam->GetZaxis()->SetTitle("2LLH_sam");

        // Scan over the parameter space
        for (int x = 0; x < n_points; x++)
        {
          if (x % countwidth == 0)
            MaCh3Utils::PrintProgressBar(x, n_points);

          for (int y = 0; y < n_points; y++)
          {
            // For PCA we have to do it differently
            if (IsPCA) {
              (*it)->setParProp_PCA(i, hScanSam->GetXaxis()->GetBinCenter(x+1));
              (*it)->setParProp_PCA(j, hScanSam->GetYaxis()->GetBinCenter(y+1));
            } else {
              // Set the parameter
              (*it)->setParProp(i, hScanSam->GetXaxis()->GetBinCenter(x+1));
              (*it)->setParProp(j, hScanSam->GetYaxis()->GetBinCenter(y+1));
            }

            // Reweight the MC
            double *fake = 0;
            for(unsigned int ivs = 0; ivs < samples.size(); ivs++) {
              samples[ivs]->reweight(fake);
            }

            // Get the -log L likelihoods
            double samplellh = 0;
            for(unsigned int ivs = 0; ivs < samples.size(); ivs++) {
              samplellh += samples[ivs]->GetLikelihood();
            }
            hScanSam->SetBinContent(x+1, y+1, 2*samplellh);
          }// end loop over y points
        } // end loop over x points

        Sample_2DLLH->cd();
        hScanSam->Write();

        delete hScanSam;
        hScanSam = nullptr;

        // Reset the parameters to their prior central values
        if (IsPCA) {
          (*it)->setParProp_PCA(i, prior_x);
          (*it)->setParProp_PCA(j, prior_y);
        } else {
          (*it)->setParProp(i, prior_x);
          (*it)->setParProp(j, prior_y);
        }
      } //end loop over systematics y
    }//end loop over systematics X
  }//end loop covariance classes

  Sample_2DLLH->Write();
  delete Sample_2DLLH;

  SaveOutput();
}
