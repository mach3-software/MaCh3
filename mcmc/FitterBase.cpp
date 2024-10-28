#include "FitterBase.h"

#include "TRandom.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"

// *************************
// Initialise the manager and make it an object of FitterBase class
// Now we can dump manager settings to the output file
FitterBase::FitterBase(manager * const man) : fitMan(man) {
// *************************

  //Get mach3 modes from manager
  Modes = fitMan->GetMaCh3Modes();
  random = std::make_unique<TRandom3>(fitMan->raw()["General"]["Seed"].as<int>());

  // Counter of the accepted # of steps
  accCount = 0;
  step = 0;

  clock = std::make_unique<TStopwatch>();
  stepClock = std::make_unique<TStopwatch>();
  #ifdef DEBUG
  // Fit summary and debug info
  debug = GetFromManager<bool>(fitMan->raw()["General"]["Debug"], false);
  #endif

  auto outfile = fitMan->raw()["General"]["OutputFile"].as<std::string>();
  // Save output every auto_save steps
  //you don't want this too often https://root.cern/root/html606/TTree_8cxx_source.html#l01229
  auto_save = fitMan->raw()["General"]["MCMC"]["AutoSave"].as<int>();
  // Do we want to save the nominal parameters to output
  save_nominal = true;

  // Do we want to save proposal? This will break plotting scripts and is heave for disk space and step time. Only use when debugging
  SaveProposal = false;

  #ifdef MULTITHREAD
  //KS: TODO This should help with performance when saving entries to ROOT file. I didn't have time to validate hence commented out
  //Based on other tests it is really helpful
  //ROOT::EnableImplicitMT();
  #endif
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

  TotalNSamples = 0;
  fTestLikelihood = GetFromManager<bool>(fitMan->raw()["General"]["Fitter"]["FitTestLikelihood"], false);
}

// *************************
// Destructor: close the logger and output file
FitterBase::~FitterBase() {
// *************************
  SaveOutput();

  if(outputFile != nullptr) delete outputFile;
  outputFile = nullptr;
  MACH3LOG_DEBUG("Closing MaCh3 Fitter Engine");
}

// *******************
// Prepare the output tree
void FitterBase::SaveSettings() {
// *******************
  if(SettingsSaved) return;

  outputFile->cd();

  TDirectory* MaCh3Version = outputFile->mkdir("MaCh3Engine");
  MaCh3Version->cd();

  if (std::getenv("MaCh3_ROOT") == NULL) {
    MACH3LOG_ERROR("Need MaCh3_ROOT environment variable");
    MACH3LOG_ERROR("Please remember about source bin/setup.MaCh3.sh");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (std::getenv("MACH3") == NULL) {
    MACH3LOG_ERROR("Need MACH3 environment variable");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  std::string header_path = std::string(std::getenv("MACH3"));
  header_path += "/version.h";
  FILE* file = fopen(header_path.c_str(), "r");
  //KS: It is better to use experiment specific header file. If given experiment didn't provide it we gonna use one given by Core MaCh3.
  if (!file)
  {
    header_path = std::string(std::getenv("MaCh3_ROOT"));
    header_path += "/version.h";
  }
  else
  {
    fclose(file);
  }

  // EM: embed the cmake generated version.h file
  TMacro versionHeader("version_header", "version_header");
  versionHeader.ReadFile(header_path.c_str());
  versionHeader.Write();

  TNamed Engine(GetName(), GetName());
  Engine.Write(GetName().c_str());

  MaCh3Version->Write();
  delete MaCh3Version;

  outputFile->cd();

  fitMan->SaveSettings(outputFile);

  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  MACH3LOG_INFO("#####Current Setup#####");
  MACH3LOG_INFO("Number of covariances: {}", systematics.size());
  for(unsigned int i = 0; i < systematics.size(); ++i)
    MACH3LOG_INFO("{}: Cov name: {}, it has {} params", i, systematics[i]->getName(), systematics[i]->GetNumParams());
  MACH3LOG_INFO("Number of SamplePDFs: {}", samples.size());
  for(unsigned int i = 0; i < samples.size(); ++i)
    MACH3LOG_INFO("{}: SamplePDF name: {}, it has {} samples",i , samples[i]->GetName(), samples[i]->GetNsamples());

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
      MACH3LOG_CRITICAL("No samples Found! Is this really what you wanted?");
      //throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    // Prepare the output trees
    for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it) {
      (*it)->SetBranches(*outTree, SaveProposal);
    }

    outTree->Branch("LogL", &logLCurr, "LogL/D");
    outTree->Branch("accProb", &accProb, "accProb/D");
    outTree->Branch("step", &step, "step/I");
    outTree->Branch("stepTime", &stepTime, "stepTime/D");

    // Store individual likelihood components
    // Using double means double as large output file!
    sample_llh.resize(samples.size());
    syst_llh.resize(systematics.size());

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
  TotalNSamples += sample->GetNsamples();
  MACH3LOG_INFO("Adding {} object, with {} samples", sample->GetName(), sample->GetNsamples());
  samples.push_back(sample);

  return;
}

// *************************
// Add flux systematics, cross-section systematics, ND systematics to the chain
void FitterBase::addSystObj(covarianceBase * const cov) {
// *************************
  MACH3LOG_INFO("Adding systematic object {}, with {} params", cov->getName(), cov->GetNumParams());
  systematics.push_back(cov);

  CovFolder->cd();
  std::vector<double> n_vec(cov->GetNumParams());
  for (int i = 0; i < cov->GetNumParams(); ++i)
    n_vec[i] = cov->getParInit(i);

  TVectorT<double> t_vec(cov->GetNumParams(), n_vec.data());
  t_vec.Write((std::string(cov->getName()) + "_prior").c_str());

  cov->getCovMatrix()->Write(cov->getName().c_str());

  TH2D* CorrMatrix = cov->GetCorrelationMatrix();
  CorrMatrix->Write((cov->getName() + std::string("_Corr")).c_str());
  delete CorrMatrix;

  // If we have yaml config file for covariance let's save it
  YAML::Node Config = cov->GetConfig();
  if(!Config.IsNull())
  {
    TMacro ConfigSave = YAMLtoTMacro(Config, (std::string("Config_") + cov->getName()));
    ConfigSave.Write();
  }

  outputFile->cd();

  return;
}


// *******************
void FitterBase::StartFromPreviousFit(const std::string& FitName) {
// *******************
  MACH3LOG_INFO("Getting starting position from {}", FitName);

  TFile *infile = new TFile(FitName.c_str(), "READ");
  TTree *posts = (TTree*)infile->Get("posteriors");
  int step_val = 0;
  double log_val = _LARGE_LOGL_;
  posts->SetBranchAddress("step",&step_val);
  posts->SetBranchAddress("LogL",&log_val);

  for (size_t s = 0; s < systematics.size(); ++s)
  {
    TDirectory* CovarianceFolder = (TDirectory*)infile->Get("CovarianceFolder");

    std::string ConfigName = "Config_" + std::string(systematics[s]->getName());
    TMacro *ConfigCov = (TMacro*)(CovarianceFolder->Get(ConfigName.c_str()));
    // KS: Not every covariance uses yaml, if it uses yaml make sure they are identical
    if (ConfigCov != nullptr) {
      // Config which was in MCMC from which we are starting
      YAML::Node CovSettings = TMacroToYAML(*ConfigCov);
      // Config from currently used cov object
      YAML::Node ConfigCurrent = systematics[s]->GetConfig();

      if (!compareYAMLNodes(CovSettings, ConfigCurrent))
      {
        MACH3LOG_ERROR("Yaml configs in previous chain and current one are different", FitName);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
    }
    if(ConfigCov == nullptr) delete ConfigCov;

    std::vector<double> branch_vals(systematics[s]->GetNumParams());
    for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
      branch_vals[i] = _BAD_DOUBLE_;
      posts->SetBranchAddress(systematics[s]->GetParName(i).c_str(), &branch_vals[i]);
    }
    posts->GetEntry(posts->GetEntries()-1);

    for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
      if(branch_vals[i] == _BAD_DOUBLE_)
      {
        MACH3LOG_ERROR("Parameter {} is unvitalised with value {}", i, branch_vals[i]);
        MACH3LOG_ERROR("Please check more precisely chain you passed {}", FitName);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
    }
    systematics[s]->setParameters(branch_vals);
    systematics[s]->acceptStep();

    MACH3LOG_INFO("Printing new starting values for: {}", systematics[s]->getName());
    systematics[s]->printNominalCurrProp();

    CovarianceFolder->Close();
    delete CovarianceFolder;
  }
  logLCurr = log_val;

  infile->Close();
  delete infile;

  for (size_t s = 0; s < systematics.size(); ++s) {
    if(systematics[s]->getDoAdaption()){ //Use separate throw matrix for xsec
      systematics[s]->setNumberOfSteps(step_val);
    }
  }
}

// *******************
// Process the MCMC output to get postfit etc
void FitterBase::ProcessMCMC() {
// *******************
  if (fitMan == nullptr) return;

  // Process the MCMC
  if (fitMan->raw()["General"]["ProcessMCMC"].as<bool>()) {

    // Make the processor
    MCMCProcessor Processor(std::string(outputFile->GetName()));

    Processor.Initialise();
    // Make the TVectorD pointers which hold the processed output
    TVectorD *Central = nullptr;
    TVectorD *Errors = nullptr;
    TVectorD *Central_Gauss = nullptr;
    TVectorD *Errors_Gauss = nullptr;
    TVectorD *Peaks = nullptr;

    // Make the postfit
    Processor.GetPostfit(Central, Errors, Central_Gauss, Errors_Gauss, Peaks);
    Processor.DrawPostfit();

    // Make the TMatrix pointers which hold the processed output
    TMatrixDSym *Covariance = nullptr;
    TMatrixDSym *Correlation = nullptr;

    // Make the covariance matrix
    Processor.GetCovariance(Covariance, Correlation);
    Processor.DrawCovariance();

    std::vector<TString> BranchNames = Processor.GetBranchNames();

    // Re-open the TFile
    if (!outputFile->IsOpen()) {
      MACH3LOG_INFO("Opening output again to update with means..");
      outputFile = new TFile(fitMan->raw()["General"]["Output"]["Filename"].as<std::string>().c_str(), "UPDATE");
    }
    Central->Write("PDF_Means");
    Errors->Write("PDF_Errors");
    Central_Gauss->Write("Gauss_Means");
    Errors_Gauss->Write("Errors_Gauss");
    Covariance->Write("Covariance");
    Correlation->Write("Correlation");
  }
}

// *************************
// Run Drag Race
void FitterBase::DragRace(const int NLaps) {
// *************************

  MACH3LOG_INFO("Let the Race Begin!");
  // Reweight the MC
  for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
  {
    TStopwatch clockRace;
    clockRace.Start();
    for(int Lap = 0; Lap < NLaps; Lap++) {
      samples[ivs]->reweight();
    }
    clockRace.Stop();
    MACH3LOG_INFO("It took {:.4f} s to reweights {} times sample: {}", clockRace.RealTime(), NLaps, samples[ivs]->GetName());
    MACH3LOG_INFO("On average {:.6f}", clockRace.RealTime()/NLaps);
  }

  for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
  {
    TStopwatch clockRace;
    clockRace.Start();
    for(int Lap = 0; Lap < NLaps; Lap++) {
      samples[ivs]->GetLikelihood();
    }
    clockRace.Stop();
    MACH3LOG_INFO("It took {:.4f} s to calculate  GetLikelihood {} times sample:  {}", clockRace.RealTime(), NLaps, samples[ivs]->GetName());
    MACH3LOG_INFO("On average {:.6f}", clockRace.RealTime()/NLaps);
  }
  // Get vector of proposed steps. If we want to run LLH scan or something else after we need to revert changes after proposing steps multiple times
  std::vector<std::vector<double>> StepsValuesBefore(systematics.size());
  for (size_t s = 0; s < systematics.size(); ++s) {
    StepsValuesBefore[s] = systematics[s]->getProposed();
  }
  for (size_t s = 0; s < systematics.size(); ++s) {
    TStopwatch clockRace;
    clockRace.Start();
    for(int Lap = 0; Lap < NLaps; Lap++) {
      systematics[s]->proposeStep();
    }
    clockRace.Stop();
    MACH3LOG_INFO("It took {:.4f} s to propose step {} times cov:  {}",  clockRace.RealTime(), NLaps, systematics[s]->getName());
    MACH3LOG_INFO("On average {:.6f}", clockRace.RealTime()/NLaps);
  }
  for (size_t s = 0; s < systematics.size(); ++s) {
    systematics[s]->setParameters(StepsValuesBefore[s]);
  }

  for (size_t s = 0; s < systematics.size(); ++s) {
    TStopwatch clockRace;
    clockRace.Start();
    for(int Lap = 0; Lap < NLaps; Lap++) {
      systematics[s]->GetLikelihood();
    }
    clockRace.Stop();
    MACH3LOG_INFO("It took {:.4f} s to calculate  get likelihood {} times cov:  {}",  clockRace.RealTime(), NLaps, systematics[s]->getName());
    MACH3LOG_INFO("On average {:.6f}", clockRace.RealTime()/NLaps);
  }
  MACH3LOG_INFO("End of race");
}

// *************************
// Run LLH scan
void FitterBase::RunLLHScan() {
// *************************

  // Save the settings into the output file
  SaveSettings();

  MACH3LOG_INFO("Starting LLH Scan");

  //KS: Turn it on if you want LLH scan for each ND sample separately, which increase time significantly but can be useful for validating new samples or dials.
  bool PlotAllNDsamplesLLH = false;
  if(fitMan->raw()["General"]["LLHScanBySample"])
    PlotAllNDsamplesLLH = fitMan->raw()["General"]["LLHScanBySample"].as<bool>();

  std::vector<std::string> SkipVector;
  if(fitMan->raw()["General"]["LLHScanSkipVector"])
  {
    SkipVector = fitMan->raw()["General"]["LLHScanSkipVector"].as<std::vector<std::string>>();
    MACH3LOG_INFO("Found skip vector with {} entries", SkipVector.size());
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

  std::vector<TDirectory *>SampleSplit_LLH;
  if(PlotAllNDsamplesLLH)
  {
    SampleSplit_LLH.resize(TotalNSamples);
    int SampleIterator = 0;
    for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
    {
      for(_int_ is = 0; is < samples[ivs]->GetNsamples(); is++ )
      {
        SampleSplit_LLH[SampleIterator] = outputFile->mkdir((samples[ivs]->GetSampleName(is)+ "_LLH").c_str());
        SampleIterator++;
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
    int npars = (*it)->GetNumParams();
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
        MACH3LOG_INFO("eval {} = {:.2f}", i, (*it)->getEigenValues()(i));
        MACH3LOG_INFO("prior {} = {:.2f}", i, prior);
        MACH3LOG_INFO("lower {} = {:.2f}", i, lower);
        MACH3LOG_INFO("upper {} = {:.2f}", i, upper);
        MACH3LOG_INFO("nSigma = {:.2f}", nSigma);
      }

      // Cross-section and flux parameters have boundaries that we scan between, check that these are respected in setting lower and upper variables
      if (lower < (*it)->GetLowerBound(i)) {
        lower = (*it)->GetLowerBound(i);
      }
      if (upper > (*it)->GetUpperBound(i)) {
        upper = (*it)->GetUpperBound(i);
      }
      MACH3LOG_INFO("Scanning {} with {} steps, from {:.2f} - {:.2f}, prior = {:.2f}", name, n_points, lower, upper, prior);

      // Make the TH1D
      auto hScan = std::make_unique<TH1D>((name + "_full").c_str(), (name + "_full").c_str(), n_points, lower, upper);
      hScan->SetTitle(std::string(std::string("2LLH_full, ") + name + ";" + name + "; -2(ln L_{sample} + ln L_{xsec+flux} + ln L_{det})").c_str());

      auto hScanSam = std::make_unique<TH1D>((name + "_sam").c_str(), (name + "_sam").c_str(), n_points, lower, upper);
      hScanSam->SetTitle(std::string(std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());

      std::vector<TH1D*> hScanSample(samples.size());
      std::vector<double> nSamLLH(samples.size());
      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        std::string NameTemp = samples[ivs]->GetName();
        hScanSample[ivs] = new TH1D((name+"_"+NameTemp).c_str(), (name+"_" + NameTemp).c_str(), n_points, lower, upper);
        hScanSample[ivs]->SetTitle(std::string(std::string("2LLH_" + NameTemp + ", ") + name + ";" + name + "; -2(ln L_{" + NameTemp +"})").c_str());
        nSamLLH[ivs] = 0.;
      }

      std::vector<TH1D*> hScanCov(systematics.size());
      std::vector<double> nCovLLH(systematics.size());
      for(unsigned int ivc = 0; ivc < systematics.size(); ivc++ )
      {
        std::string NameTemp = systematics[ivc]->getName();
        NameTemp = NameTemp.substr(0, NameTemp.find("_cov"));

        hScanCov[ivc] = new TH1D((name+"_"+NameTemp).c_str(), (name+"_" + NameTemp).c_str(), n_points, lower, upper);
        hScanCov[ivc]->SetTitle(std::string(std::string("2LLH_" + NameTemp + ", ") + name + ";" + name + "; -2(ln L_{" + NameTemp +"})").c_str());
        nCovLLH[ivc] = 0.;
      }

      std::vector<TH1D *> hScanSamSplit;
      std::vector<double> sampleSplitllh;
      if(PlotAllNDsamplesLLH)
      {
        int SampleIterator = 0;
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          hScanSamSplit.resize(TotalNSamples);
          sampleSplitllh.resize(TotalNSamples);
          for(_int_ is = 0; is < samples[ivs]->GetNsamples(); is++ )
          {
            hScanSamSplit[SampleIterator] = new TH1D((name+samples[ivs]->GetSampleName(is)).c_str(), (name+samples[ivs]->GetSampleName(is)).c_str(), n_points, lower, upper);
            hScanSamSplit[SampleIterator]->SetTitle(std::string(std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());
            SampleIterator++;
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
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          samples[ivs]->reweight();
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
          int SampleIterator = 0;
          for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
          {
            for(_int_ is = 0; is < samples[ivs]->GetNsamples(); is++)
            {
              sampleSplitllh[SampleIterator] = samples[ivs]->getSampleLikelihood(is);
              SampleIterator++;
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
          int SampleIterator = 0;
          for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
          {
            for(_int_ is = 0; is < samples[ivs]->GetNsamples(); is++)
            {
              hScanSamSplit[is]->SetBinContent(j+1, 2*sampleSplitllh[is]);
              SampleIterator++;
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

      if(PlotAllNDsamplesLLH)
      {
        int SampleIterator = 0;
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          for(_int_ is = 0; is < samples[ivs]->GetNsamples(); is++)
          {
            SampleSplit_LLH[SampleIterator]->cd();
            hScanSamSplit[SampleIterator]->Write();
            delete hScanSamSplit[SampleIterator];
            SampleIterator++;
          }
        }
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
    int SampleIterator = 0;
    for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
    {
      for(_int_ is = 0; is < samples[ivs]->GetNsamples(); is++ )
      {
        SampleSplit_LLH[SampleIterator]->Write();
        delete SampleSplit_LLH[SampleIterator];
        SampleIterator++;
      }
    }
  }
}

// *************************
//LLH scan is good first estimate of step scale
void FitterBase::GetStepScaleBasedOnLLHScan() {
// *************************
  TDirectory *Sample_LLH = outputFile->Get<TDirectory>("Sample_LLH");
  MACH3LOG_INFO("Starting Get Step Scale Based On LLHScan");

  if(Sample_LLH->IsZombie())
  {
    MACH3LOG_WARN("Couldn't find Sample_LLH, it looks like LLH scan wasn't run, will do this now");
    RunLLHScan();
  }

  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    bool isxsec = (std::string((*it)->getName()) == "xsec_cov");

    const int npars = (*it)->GetNumParams();
    std::vector<double> StepScale(npars);
    for (int i = 0; i < npars; ++i)
    {
      std::string name = (*it)->GetParName(i);
      // For xsec we can get the actual name, hurray for being informative
      if (isxsec) name = (*it)->GetParFancyName(i);

      StepScale[i] = (*it)->GetIndivStepScale(i);
      TH1D* LLHScan = static_cast<TH1D*>(Sample_LLH->Get((name + "_sam").c_str()));
      if(LLHScan == nullptr)
      {
        MACH3LOG_WARN("Couldn't find LLH scan, for {}, skipping", name);
        continue;
      }
      const double LLH_val = std::max(LLHScan->GetBinContent(1), LLHScan->GetBinContent(LLHScan->GetNbinsX()));
      //If there is no sensitivity leave it
      if(LLH_val < 0.001) continue;

      // EM: assuming that the likelihood is gaussian, approximate sigma value is given by variation/sqrt(-2LLH)
      // can evaluate this at any point, simple to evaluate it in the first bin of the LLH scan
      // KS: We assume variation is 1 sigma, each dial has different scale so it becomes faff...
      const double Var = 1.;
      const double approxSigma = TMath::Abs(Var)/std::sqrt(LLH_val);

      // Based on Ewan comment I just took the 1sigma width from the LLH, assuming it was Gaussian, but then had to also scale by 2.38/sqrt(N_params)
      const double NewStepScale = approxSigma * 2.38/std::sqrt(npars);
      StepScale[i] = NewStepScale;
      MACH3LOG_DEBUG("Sigma: {}", approxSigma);
      MACH3LOG_DEBUG("optimal Step Size: {}", NewStepScale);
    }
    (*it)->setIndivStepScale(StepScale);
    (*it)->SaveUpdatedMatrixConfig();
  }
}

// *************************
// Run 2D LLH scan
void FitterBase::Run2DLLHScan() {
// *************************

  // Save the settings into the output file
  SaveSettings();

  MACH3LOG_INFO("Starting 2D LLH Scan");

  TDirectory *Sample_2DLLH = outputFile->mkdir("Sample_2DLLH");
  std::vector<std::string> SkipVector;
  if(fitMan->raw()["General"]["LLHScanSkipVector"])
  {
    SkipVector = fitMan->raw()["General"]["LLHScanSkipVector"].as<std::vector<std::string>>();
    MACH3LOG_INFO("Found skip vector with {} entries", SkipVector.size());
  }

  // Number of points we do for each LLH scan
  constexpr int n_points = 20;
  // We print 5 reweights
  constexpr int countwidth = double(n_points)/double(5);

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
    int npars = (*it)->GetNumParams();
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
        lower_x = prior_x - nSigma*std::sqrt(((*it)->getEigenValues())(i));
        upper_x = prior_x + nSigma*std::sqrt(((*it)->getEigenValues())(i));
        MACH3LOG_INFO("eval {} = {:.2f}", i, (*it)->getEigenValues()(i));
        MACH3LOG_INFO("prior {} = {:.2f}", i, prior_x);
        MACH3LOG_INFO("lower {} = {:.2f}", i, lower_x);
        MACH3LOG_INFO("upper {} = {:.2f}", i, upper_x);
        MACH3LOG_INFO("nSigma = {:.2f}", nSigma);
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

        skip = false;
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
          MACH3LOG_INFO("eval {} = {:.2f}", i, (*it)->getEigenValues()(j));
          MACH3LOG_INFO("prior {} = {:.2f}", i, prior_y);
          MACH3LOG_INFO("lower {} = {:.2f}", i, lower_y);
          MACH3LOG_INFO("upper {} = {:.2f}", i, upper_y);
          MACH3LOG_INFO("nSigma = {:.2f}", nSigma);
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

        auto hScanSam = std::make_unique<TH2D>((name_x + "_" + name_y + "_sam").c_str(), (name_x + "_" + name_y + "_sam").c_str(),
                                                n_points, lower_x, upper_x, n_points, lower_y, upper_y);
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
            //double *fake = 0;
            for(unsigned int ivs = 0; ivs < samples.size(); ivs++) {
              samples[ivs]->reweight();
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
}

// *************************
void FitterBase::RunSigmaVar() {
// *************************
  // Save the settings into the output file
  SaveSettings();

  MACH3LOG_INFO("Starting Sigma Variation");

  // Number of variations we want
  constexpr int numVar = 5;
  //-3 -1 0 +1 +3 sigma variation
  constexpr int sigmaArray[numVar] = {-3, -1, 0, 1, 3};

  outputFile->cd();

  //KS: this is only relevant if PlotByMode is turned on
  //Checking each mode is time consuming so we only consider one which are relevant for particular analysis
  constexpr int nRelevantModes = 2;
  MaCh3Modes_t RelevantModes[nRelevantModes] = {Modes->GetMode("CCQE"), Modes->GetMode("2p2h")};
  bool DoByMode = GetFromManager<int>(fitMan->raw()["General"]["DoByMode"], false);

  //KS: If true it will make additional plots with LLH sample contribution in each bin, should make it via config file...
  bool PlotLLHperBin = false;

  std::vector<std::string> SkipVector;
  if(fitMan->raw()["General"]["LLHScanSkipVector"])
  {
    SkipVector = fitMan->raw()["General"]["LLHScanSkipVector"].as<std::vector<std::string>>();
    MACH3LOG_INFO("Found skip vector with {} entries", SkipVector.size());
  }

  bool isxsec = false;
  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    TMatrixDSym *Cov = (*it)->getCovMatrix();

    if((*it)->IsPCA())
    {
      MACH3LOG_ERROR("Using PCAed matrix not implemented within sigma var code, I am sorry :(");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    if (std::string((*it)->getName()) == "xsec_cov")
    {
      isxsec = true;
    } else {
      isxsec = false;
    }

    // Loop over xsec parameters
    for (int i = 0; i < (*it)->GetNumParams(); ++i)
    {
      // Get the parameter name
      std::string name = (*it)->GetParName(i);
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

      outputFile->cd();
      TDirectory * dirArryDial = outputFile->mkdir(name.c_str());
      TDirectory **dirArrySample = new TDirectory*[TotalNSamples];

      int SampleIterator = 0;
      // Get each sample and how it's responded to our reweighted parameter
      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        for(_int_ k = 0; k < samples[ivs]->GetNsamples(); k++ )
        {
          std::string title = std::string(samples[ivs]->getPDF(k)->GetName());
          dirArryDial->cd();
          dirArrySample[SampleIterator] = dirArryDial->mkdir(title.c_str());
          SampleIterator++;
        }
      }

      // Get the initial value of ith parameter
      double init = (*it)->getParInit(i);

      TH1D ***sigmaArray_x = new TH1D**[numVar]();
      TH1D ***sigmaArray_y = new TH1D**[numVar]();
      TH1D ***sigmaArray_x_norm = new TH1D**[numVar]();
      TH1D ***sigmaArray_y_norm = new TH1D**[numVar]();

      // Set up for single mode
      TH1D ****sigmaArray_mode_x = nullptr;
      TH1D ****sigmaArray_mode_y = nullptr;
      if (DoByMode)
      {
        sigmaArray_mode_x = new TH1D***[numVar]();
        sigmaArray_mode_y = new TH1D***[numVar]();
      }

      MACH3LOG_INFO("{:<20}{}", name, "");

      // Make asymmetric errors just in case
      for (int j = 0; j < numVar; ++j)
      {
        // New value = prior + variation*1sigma uncertainty
        double paramVal = (*it)->getParInit(i)+sigmaArray[j]*std::sqrt((*Cov)(i,i));

        // Check the bounds on the parameter
        if (paramVal > (*it)->GetUpperBound(i)) {
          paramVal = (*it)->GetUpperBound(i);
        } else if (paramVal < (*it)->GetLowerBound(i)) {
          paramVal = (*it)->GetLowerBound(i);
        }

        // Set the parameter
        (*it)->setParProp(i, paramVal);
        // And reweight the sample
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++) {
          samples[ivs]->reweight();
        }

        sigmaArray_x[j] = new TH1D*[TotalNSamples]();
        sigmaArray_y[j] = new TH1D*[TotalNSamples]();
        sigmaArray_x_norm[j] = new TH1D*[TotalNSamples]();
        sigmaArray_y_norm[j] = new TH1D*[TotalNSamples]();

        if (DoByMode)
        {
          sigmaArray_mode_x[j] = new TH1D**[TotalNSamples]();
          sigmaArray_mode_y[j] = new TH1D**[TotalNSamples]();
        }

        SampleIterator = 0;
        // Get each sample and how it's responded to our reweighted parameter
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          for (_int_ k = 0; k < samples[ivs]->GetNsamples(); ++k)
          {
            // Make a string of the double
            std::ostringstream ss;
            ss << paramVal;
            std::string parVarTitle = std::string(name) + "_" + ss.str();

            // This is a TH2D
            TH2Poly* currSamp = static_cast<TH2Poly*>(samples[ivs]->getPDF(k)->Clone());
            currSamp->SetDirectory(0);
            // Set a descriptiv-ish title
            std::string title_long = std::string(currSamp->GetName())+"_"+parVarTitle;

            // Enable the mode histograms AFTER addSelection is called
            //Get the 1d binning we want. Let's just use SetupBinning to get this as it already exists
            std::vector<double> xbins;
            std::vector<double> ybins;
            samples[ivs]->SetupBinning(k, xbins, ybins);

            //KS:here we loop over all reaction modes defined in "RelevantModes[nRelevantModes]"
            if (DoByMode)
            {
              sigmaArray_mode_x[j][SampleIterator] = new TH1D*[nRelevantModes]();
              sigmaArray_mode_y[j][SampleIterator] = new TH1D*[nRelevantModes]();
              // Now get the TH2D mode variations
              TH2Poly** currSampMode = new TH2Poly*[nRelevantModes]();
              std::string mode_title_long;
              for(int ir = 0; ir < nRelevantModes; ir++)
              {
                currSampMode[ir] = static_cast<TH2Poly*>(samples[ivs]->getPDFMode(k, RelevantModes[ir])->Clone());
                currSampMode[ir]->SetDirectory(0);

                mode_title_long = title_long + "_" + Modes->GetMaCh3ModeName(RelevantModes[ir]);
                currSampMode[ir]->SetNameTitle(mode_title_long.c_str(), mode_title_long.c_str());
                dirArrySample[SampleIterator]->cd();
                currSampMode[ir]->Write();

                sigmaArray_mode_x[j][SampleIterator][ir] = PolyProjectionX(currSampMode[ir], (mode_title_long+"_xProj").c_str(), xbins);
                sigmaArray_mode_x[j][SampleIterator][ir]->SetDirectory(0);
                sigmaArray_mode_y[j][SampleIterator][ir] = PolyProjectionY(currSampMode[ir], (mode_title_long+"_yProj").c_str(), ybins);
                sigmaArray_mode_y[j][SampleIterator][ir]->SetDirectory(0);
                delete currSampMode[ir];
              }
              delete[] currSampMode;
            }

            //KS: This will give different results depending if data or Asimov, both have their uses
            if (PlotLLHperBin)
            {
              TH2Poly* currLLHSamp = static_cast<TH2Poly*>(samples[ivs]->getPDF(k)->Clone());
              currLLHSamp->SetDirectory(0);
              currLLHSamp->Reset("");
              currLLHSamp->Fill(0.0, 0.0, 0.0);

              TH2Poly* MCpdf = static_cast<TH2Poly*>(samples[ivs]->getPDF(k));
              TH2Poly* Datapdf = static_cast<TH2Poly*>(samples[ivs]->getData(k));
              TH2Poly* W2pdf = static_cast<TH2Poly*>(samples[ivs]->getW2(k));
              for(int bin = 1; bin < currLLHSamp->GetNumberOfBins()+1; bin++)
              {
                const double mc = MCpdf->GetBinContent(bin);
                const double dat = Datapdf->GetBinContent(bin);
                const double w2 = W2pdf->GetBinContent(bin);
                currLLHSamp->SetBinContent(bin, samples[ivs]->getTestStatLLH(dat, mc, w2));
              }
              currLLHSamp->SetNameTitle((title_long+"_LLH").c_str() ,(title_long+"_LLH").c_str());
              dirArrySample[SampleIterator]->cd();
              currLLHSamp->Write();
              delete currLLHSamp;
            }

            // Project down onto x axis
            sigmaArray_x[j][SampleIterator] = PolyProjectionX(currSamp, (title_long+"_xProj").c_str(), xbins);
            sigmaArray_x[j][SampleIterator]->SetDirectory(0);
            sigmaArray_x[j][SampleIterator]->GetXaxis()->SetTitle(currSamp->GetXaxis()->GetTitle());
            sigmaArray_y[j][SampleIterator] = PolyProjectionY(currSamp, (title_long+"_yProj").c_str(), ybins);
            sigmaArray_y[j][SampleIterator]->SetDirectory(0);
            sigmaArray_y[j][SampleIterator]->GetXaxis()->SetTitle(currSamp->GetYaxis()->GetTitle());

            sigmaArray_x_norm[j][SampleIterator] = static_cast<TH1D*>(sigmaArray_x[j][SampleIterator]->Clone());
            sigmaArray_x_norm[j][SampleIterator]->SetDirectory(0);
            sigmaArray_x_norm[j][SampleIterator]->Scale(1., "width");
            sigmaArray_y_norm[j][SampleIterator] = static_cast<TH1D*>(sigmaArray_y[j][SampleIterator]->Clone());
            sigmaArray_y_norm[j][SampleIterator]->SetDirectory(0);
            sigmaArray_y_norm[j][SampleIterator]->Scale(1., "width");

            currSamp->SetNameTitle(title_long.c_str(), title_long.c_str());
            dirArrySample[k]->cd();
            currSamp->Write();

            sigmaArray_x[j][k]->Write();
            sigmaArray_y[j][k]->Write();
            delete currSamp;
            SampleIterator++;
          }//End loop over samples
        }
      } // End looping over variation

      // Restore the parameter to prior value
      (*it)->setParProp(i, init);

      SampleIterator = 0;
      // Get each sample and how it's responded to our reweighted parameter
      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        for (_int_ k = 0; k < samples[ivs]->GetNsamples(); ++k)
        {
          std::string title = std::string(samples[ivs]->getPDF(k)->GetName()) + "_" + name;
          TGraphAsymmErrors *var_x = MakeAsymGraph(sigmaArray_x[1][SampleIterator], sigmaArray_x[2][SampleIterator], sigmaArray_x[3][SampleIterator], (title+"_X").c_str());
          TGraphAsymmErrors *var_y = MakeAsymGraph(sigmaArray_y[1][SampleIterator], sigmaArray_y[2][SampleIterator], sigmaArray_y[3][SampleIterator], (title+"_Y").c_str());

          TGraphAsymmErrors *var_x_norm = MakeAsymGraph(sigmaArray_x_norm[1][SampleIterator], sigmaArray_x_norm[2][SampleIterator], sigmaArray_x_norm[3][SampleIterator], (title+"_X_norm").c_str());
          TGraphAsymmErrors *var_y_norm = MakeAsymGraph(sigmaArray_y_norm[1][SampleIterator], sigmaArray_y_norm[2][SampleIterator], sigmaArray_y_norm[3][SampleIterator], (title+"_Y_norm").c_str());

          dirArrySample[SampleIterator]->cd();
          var_x->Write();
          var_y->Write();
          var_x_norm->Write();
          var_y_norm->Write();
          delete var_x;
          delete var_y;
          delete var_x_norm;
          delete var_y_norm;

          //KS: here we loop over all reaction modes defined in "RelevantModes[nRelevantModes]"
          if (DoByMode)
          {
            for(int ir = 0; ir < nRelevantModes;ir++)
            {
              TGraphAsymmErrors* var_mode_x = MakeAsymGraph(sigmaArray_mode_x[1][SampleIterator][ir], sigmaArray_mode_x[2][SampleIterator][ir], sigmaArray_mode_x[3][SampleIterator][ir], (title+"_"+Modes->GetMaCh3ModeName(RelevantModes[ir])+"_X").c_str());
              TGraphAsymmErrors* var_mode_y = MakeAsymGraph(sigmaArray_mode_y[1][SampleIterator][ir], sigmaArray_mode_y[2][SampleIterator][ir], sigmaArray_mode_y[3][SampleIterator][ir], (title+"_"+Modes->GetMaCh3ModeName(RelevantModes[ir])+"_Y").c_str());

              dirArrySample[SampleIterator]->cd();
              var_mode_x->Write();
              var_mode_y->Write();

              delete var_mode_x;
              delete var_mode_y;
            } // end for nRelevantModes
          } // end if mode

          SampleIterator++;
        }//End loop over samples(k)
      }//end looping over sample object
      SampleIterator = 0;
      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        for (_int_ k = 0; k < samples[ivs]->GetNsamples(); ++k)
        {
          dirArrySample[SampleIterator]->Close();
          delete dirArrySample[SampleIterator];
          SampleIterator++;
        }
      }
      delete[] dirArrySample;

      for (int j = 0; j < numVar; ++j)
      {
        SampleIterator = 0;
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          for (_int_ k = 0; k < samples[ivs]->GetNsamples(); ++k)
          {
            delete sigmaArray_x[j][SampleIterator];
            delete sigmaArray_y[j][SampleIterator];
            delete sigmaArray_x_norm[j][SampleIterator];
            delete sigmaArray_y_norm[j][SampleIterator];
            SampleIterator++;
          }
        }
        delete[] sigmaArray_x[j];
        delete[] sigmaArray_y[j];
        delete[] sigmaArray_x_norm[j];
        delete[] sigmaArray_y_norm[j];
      }
      delete[] sigmaArray_x;
      delete[] sigmaArray_y;
      delete[] sigmaArray_x_norm;
      delete[] sigmaArray_y_norm;

      dirArryDial->Close();
      delete dirArryDial;

      if (DoByMode)
      {
        for (int j = 0; j < numVar; ++j)
        {
          SampleIterator = 0;
          for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
          {
            for (_int_ k = 0; k < samples[ivs]->GetNsamples(); ++k)
            {
              for(int ir = 0; ir < nRelevantModes;ir++)
              {
                delete sigmaArray_mode_x[j][SampleIterator][ir];
                delete sigmaArray_mode_y[j][SampleIterator][ir];
              }// end for nRelevantModes
              delete[] sigmaArray_mode_x[j][SampleIterator];
              delete[] sigmaArray_mode_y[j][SampleIterator];
              SampleIterator++;
            }// end for samples
          }
        }
        delete[] sigmaArray_mode_x;
        delete[] sigmaArray_mode_y;
      }
    } // end looping over xsec parameters (i)
  } // end looping over covarianceBase objects
}
