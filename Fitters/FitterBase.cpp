#include "FitterBase.h"
#include "Samples/SampleHandlerFD.h"

_MaCh3_Safe_Include_Start_ //{
#include "TRandom.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
_MaCh3_Safe_Include_End_ //}

#pragma GCC diagnostic ignored "-Wuseless-cast"

// *************************
// Initialise the Manager and make it an object of FitterBase class
// Now we can dump Manager settings to the output file
FitterBase::FitterBase(Manager * const man) : fitMan(man) {
// *************************
  AlgorithmName = "";
  //Get mach3 modes from Manager
  random = std::make_unique<TRandom3>(Get<int>(fitMan->raw()["General"]["Seed"], __FILE__, __LINE__));

  // Counter of the accepted # of steps
  accCount = 0;
  step = 0;
  stepStart = 0;

  clock = std::make_unique<TStopwatch>();
  stepClock = std::make_unique<TStopwatch>();
  #ifdef DEBUG
  // Fit summary and debug info
  debug = GetFromManager<bool>(fitMan->raw()["General"]["Debug"], false, __FILE__ , __LINE__);
  #endif

  auto outfile = Get<std::string>(fitMan->raw()["General"]["OutputFile"], __FILE__ , __LINE__);
  // Save output every auto_save steps
  //you don't want this too often https://root.cern/root/html606/TTree_8cxx_source.html#l01229
  auto_save = Get<int>(fitMan->raw()["General"]["MCMC"]["AutoSave"], __FILE__ , __LINE__);

  // Set the output file
  outputFile = M3::Open(outfile, "RECREATE", __FILE__, __LINE__);
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
  SampleFolder = outputFile->mkdir("SampleFolder");
  outputFile->cd();

  #ifdef DEBUG
  // Prepare the output log file
  if (debug) debugFile.open((outfile+".log").c_str());
  #endif

  TotalNSamples = 0;
  fTestLikelihood = GetFromManager<bool>(fitMan->raw()["General"]["Fitter"]["FitTestLikelihood"], false, __FILE__ , __LINE__);
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

  if (std::getenv("MaCh3_ROOT") == nullptr) {
    MACH3LOG_ERROR("Need MaCh3_ROOT environment variable");
    MACH3LOG_ERROR("Please remember about source bin/setup.MaCh3.sh");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (std::getenv("MACH3") == nullptr) {
    MACH3LOG_ERROR("Need MACH3 environment variable");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  std::string header_path = std::string(std::getenv("MACH3"));
  header_path += "/version.h";
  FILE* file = fopen(header_path.c_str(), "r");
  //KS: It is better to use experiment specific header file. If given experiment didn't provide it we gonna use one given by Core MaCh3.
  if (!file) {
    header_path = std::string(std::getenv("MaCh3_ROOT"));
    header_path += "/version.h";
  } else {
    fclose(file);
  }

  // EM: embed the cmake generated version.h file
  TMacro versionHeader("version_header", "version_header");
  versionHeader.ReadFile(header_path.c_str());
  versionHeader.Write();

  if(GetName() == ""){
    MACH3LOG_ERROR("Name of currently used algorithm is {}", GetName());
    MACH3LOG_ERROR("Have you forgotten to modify AlgorithmName?");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
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
    MACH3LOG_INFO("{}: Cov name: {}, it has {} params", i, systematics[i]->GetName(), systematics[i]->GetNumParams());
  MACH3LOG_INFO("Number of SampleHandlers: {}", samples.size());
  for(unsigned int i = 0; i < samples.size(); ++i) {
    MACH3LOG_INFO("{}: SampleHandler name: {}, it has {} samples",i , samples[i]->GetName(), samples[i]->GetNsamples());
    for(int iSam = 0; iSam < samples[i]->GetNsamples(); ++iSam) {
      MACH3LOG_INFO("   {}: Sample name: {}, with {} osc channels",iSam , samples[i]->GetSampleTitle(iSam), samples[i]->GetNOscChannels(iSam));
    }
  }
  //TN: Have to close the folder in order to write it to disk before SaveOutput is called in the destructor
  CovFolder->Close();
  SampleFolder->Close();

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

    // Do we want to save proposal? This will break plotting scripts and is heave for disk space and step time. Only use when debugging
    bool SaveProposal = GetFromManager<bool>(fitMan->raw()["General"]["SaveProposal"], false, __FILE__ , __LINE__);

    if(SaveProposal) MACH3LOG_INFO("Will save in the chain proposal parameters and LogL");
    // Prepare the output trees
    for (ParameterHandlerBase *cov : systematics) {
      cov->SetBranches(*outTree, SaveProposal);
    }

    outTree->Branch("LogL", &logLCurr, "LogL/D");
    if(SaveProposal) outTree->Branch("LogLProp", &logLProp, "LogLProp/D");
    outTree->Branch("accProb", &accProb, "accProb/D");
    outTree->Branch("step", &step, "step/i");
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
    }

    for (size_t i = 0; i < systematics.size(); ++i) {
      std::stringstream oss, oss2;
      oss << "LogL_systematic_" << systematics[i]->GetName();
      oss2 << oss.str() << "/D";
      outTree->Branch(oss.str().c_str(), &syst_llh[i], oss2.str().c_str());
    }
  }
  else
  {
    outTree->Branch("LogL", &logLCurr, "LogL/D");
    outTree->Branch("accProb", &accProb, "accProb/D");
    outTree->Branch("step", &step, "step/i");
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
void FitterBase::SanitiseInputs() {
// *******************
  for (size_t i = 0; i < samples.size(); ++i) {
    samples[i]->CleanMemoryBeforeFit();
  }
}

// *******************
void FitterBase::SaveOutput() {
// *******************
  if(FileSaved) return;
  //Stop Clock
  clock->Stop();

  //KS: Some version of ROOT keep spamming about accessing already deleted object which is wrong and not helpful...
  int originalErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  outputFile->cd();
  outTree->Write();

  MACH3LOG_INFO("{} steps took {:.2e} seconds to complete. ({:.2e}s / step).", step - stepStart, clock->RealTime(), clock->RealTime() / static_cast<double>(step - stepStart));
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

  gErrorIgnoreLevel = originalErrorLevel;
}

// *************************
// Add SampleHandler object to the Markov Chain
void FitterBase::AddSampleHandler(SampleHandlerBase * const sample) {
// *************************
  // Check if any subsample name collides with already-registered subsamples
  for (const auto &s : samples) {
    for (int iExisting = 0; iExisting < s->GetNsamples(); ++iExisting) {
      for (int iNew = 0; iNew < sample->GetNsamples(); ++iNew) {
        if (s->GetSampleTitle(iExisting) == sample->GetSampleTitle(iNew)) {
          MACH3LOG_ERROR(
            "Duplicate sample title '{}' in handler {} detected: "
            "same title exist in handler ", sample->GetSampleTitle(iNew),
            sample->GetName(), s->GetName());
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
    }
  }

  for (const auto &s : samples) {
    if (s->GetName() == sample->GetName()) {
      MACH3LOG_WARN("SampleHandler with name '{}' already exists!", sample->GetName());
      MACH3LOG_WARN("Is it intended?");
      //throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }
  // Save additional info from samples
  SampleFolder->cd();

  sample->SaveAdditionalInfo(SampleFolder);
  TotalNSamples += sample->GetNsamples();
  MACH3LOG_INFO("Adding {} object, with {} samples", sample->GetName(), sample->GetNsamples());
  samples.push_back(sample);
  outputFile->cd();
}

// *************************
// Add flux systematics, cross-section systematics, ND systematics to the chain
void FitterBase::AddSystObj(ParameterHandlerBase * const cov) {
// *************************
  MACH3LOG_INFO("Adding systematic object {}, with {} params", cov->GetName(), cov->GetNumParams());
  // KS: Need to make sure we don't have params with same name, otherwise ROOT I/O and parts of MaCh3 will be terribly confused...
  for (size_t s = 0; s < systematics.size(); ++s)
  {
    for (int iPar = 0; iPar < systematics[s]->GetNumParams(); ++iPar)
    {
      for (int i = 0; i < cov->GetNumParams(); ++i)
      {
        if(systematics[s]->GetParName(iPar) == cov->GetParName(i)){
          MACH3LOG_ERROR("ParameterHandler {} has param '{}' which already exists in in {}, with name {}",
                         cov->GetName(), cov->GetParName(i), systematics[s]->GetName(), systematics[s]->GetParName(iPar));
          throw MaCh3Exception(__FILE__ , __LINE__ );
        }
        // Same for fancy name
        if(systematics[s]->GetParFancyName(iPar) == cov->GetParFancyName(i)){
          MACH3LOG_ERROR("ParameterHandler {} has param '{}' which already exists in {}, with name {}",
                         cov->GetName(), cov->GetParFancyName(i), systematics[s]->GetName(), systematics[s]->GetParFancyName(iPar));
          throw MaCh3Exception(__FILE__ , __LINE__ );
        }
      }
    }
  }

  systematics.push_back(cov);

  CovFolder->cd();
  std::vector<double> n_vec(cov->GetNumParams());
  for (int i = 0; i < cov->GetNumParams(); ++i) {
    n_vec[i] = cov->GetParInit(i);
  }
  cov->GetCovMatrix()->Write(cov->GetName().c_str());

  TH2D* CorrMatrix = cov->GetCorrelationMatrix();
  CorrMatrix->Write((cov->GetName() + std::string("_Corr")).c_str());
  delete CorrMatrix;

  // If we have yaml config file for covariance let's save it
  YAML::Node Config = cov->GetConfig();
  if(!Config.IsNull())
  {
    TMacro ConfigSave = YAMLtoTMacro(Config, (std::string("Config_") + cov->GetName()));
    ConfigSave.Write();
  }

  outputFile->cd();
}

// *******************
void FitterBase::StartFromPreviousFit(const std::string& FitName) {
// *******************
  MACH3LOG_INFO("Getting starting position from {}", FitName);
  TFile *infile = M3::Open(FitName, "READ", __FILE__, __LINE__);
  TTree *posts = infile->Get<TTree>("posteriors");
  unsigned int step_val = 0;
  double log_val = M3::_LARGE_LOGL_;
  posts->SetBranchAddress("step",&step_val);
  posts->SetBranchAddress("LogL",&log_val);

  for (size_t s = 0; s < systematics.size(); ++s)
  {
    TDirectory* CovarianceFolder = infile->Get<TDirectory>("CovarianceFolder");

    std::string ConfigName = "Config_" + systematics[s]->GetName();
    TMacro *ConfigCov = CovarianceFolder->Get<TMacro>(ConfigName.c_str());
    // KS: Not every covariance uses yaml, if it uses yaml make sure they are identical
    if (ConfigCov != nullptr) {
      // Config which was in MCMC from which we are starting
      YAML::Node CovSettings = TMacroToYAML(*ConfigCov);
      // Config from currently used cov object
      YAML::Node ConfigCurrent = systematics[s]->GetConfig();

      if (!compareYAMLNodes(CovSettings, ConfigCurrent))
      {
        MACH3LOG_ERROR("Yaml configs in previous chain (from path {}) and current one are different", FitName);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }

      delete ConfigCov;
    }

    CovarianceFolder->Close();
    delete CovarianceFolder;

    std::vector<double> branch_vals;
    std::vector<std::string> branch_name;
    systematics[s]->MatchMaCh3OutputBranches(posts, branch_vals, branch_name);
    posts->GetEntry(posts->GetEntries()-1);

    systematics[s]->SetParameters(branch_vals);
    systematics[s]->AcceptStep();

    MACH3LOG_INFO("Printing new starting values for: {}", systematics[s]->GetName());
    systematics[s]->PrintNominalCurrProp();

    // Resetting branch addressed to nullptr as we don't want to write into a delected vector out of scope...
    for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
      posts->SetBranchAddress(systematics[s]->GetParName(i).c_str(), nullptr);
    }
  }
  logLCurr = log_val;
  logLProp = log_val;
  delete posts;
  infile->Close();
  delete infile;

  for (size_t s = 0; s < systematics.size(); ++s) {
    if(systematics[s]->GetDoAdaption()){ //Use separate throw matrix for xsec
      systematics[s]->SetNumberOfSteps(step_val);
    }
  }
}

// *******************
// Process the MCMC output to get postfit etc
void FitterBase::ProcessMCMC() {
// *******************
  if (fitMan == nullptr) return;

  // Process the MCMC
  if (GetFromManager<bool>(fitMan->raw()["General"]["ProcessMCMC"], false, __FILE__ , __LINE__)){
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
      outputFile = new TFile(Get<std::string>(fitMan->raw()["General"]["OutputFile"], __FILE__, __LINE__).c_str(), "UPDATE");
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
  MACH3LOG_INFO("All tests will be performed with {} threads", M3::GetNThreads());

  // Reweight the MC
  for(unsigned int ivs = 0; ivs < samples.size(); ++ivs)
  {
    TStopwatch clockRace;
    clockRace.Start();
    for(int Lap = 0; Lap < NLaps; ++Lap) {
      samples[ivs]->Reweight();
    }
    clockRace.Stop();
    MACH3LOG_INFO("It took {:.4f} s to reweights {} times sample: {}", clockRace.RealTime(), NLaps, samples[ivs]->GetName());
    MACH3LOG_INFO("On average {:.6f}", clockRace.RealTime()/NLaps);
  }

  for(unsigned int ivs = 0; ivs < samples.size(); ++ivs)
  {
    TStopwatch clockRace;
    clockRace.Start();
    for(int Lap = 0; Lap < NLaps; ++Lap) {
      samples[ivs]->GetLikelihood();
    }
    clockRace.Stop();
    MACH3LOG_INFO("It took {:.4f} s to calculate  GetLikelihood {} times sample:  {}", clockRace.RealTime(), NLaps, samples[ivs]->GetName());
    MACH3LOG_INFO("On average {:.6f}", clockRace.RealTime()/NLaps);
  }
  // Get vector of proposed steps. If we want to run LLH scan or something else after we need to revert changes after proposing steps multiple times
  std::vector<std::vector<double>> StepsValuesBefore(systematics.size());
  for (size_t s = 0; s < systematics.size(); ++s) {
    StepsValuesBefore[s] = systematics[s]->GetProposed();
  }
  for (size_t s = 0; s < systematics.size(); ++s) {
    TStopwatch clockRace;
    clockRace.Start();
    for(int Lap = 0; Lap < NLaps; ++Lap) {
      systematics[s]->ProposeStep();
    }
    clockRace.Stop();
    MACH3LOG_INFO("It took {:.4f} s to propose step {} times cov:  {}",  clockRace.RealTime(), NLaps, systematics[s]->GetName());
    MACH3LOG_INFO("On average {:.6f}", clockRace.RealTime()/NLaps);
  }
  for (size_t s = 0; s < systematics.size(); ++s) {
    systematics[s]->SetParameters(StepsValuesBefore[s]);
  }

  for (size_t s = 0; s < systematics.size(); ++s) {
    TStopwatch clockRace;
    clockRace.Start();
    for(int Lap = 0; Lap < NLaps; ++Lap) {
      systematics[s]->GetLikelihood();
    }
    clockRace.Stop();
    MACH3LOG_INFO("It took {:.4f} s to calculate  get likelihood {} times cov:  {}",  clockRace.RealTime(), NLaps, systematics[s]->GetName());
    MACH3LOG_INFO("On average {:.6f}", clockRace.RealTime()/NLaps);
  }
  MACH3LOG_INFO("End of race");
}

// *************************
bool FitterBase::GetScanRange(std::map<std::string, std::vector<double>>& scanRanges) const {
// *************************
  bool isScanRanges = false;
  // YSP: Set up a mapping to store parameters with user-specified ranges, suggested by D. Barrow
  if(fitMan->raw()["LLHScan"]["ScanRanges"]){
    YAML::Node scanRangesList = fitMan->raw()["LLHScan"]["ScanRanges"];
    for (auto it = scanRangesList.begin(); it != scanRangesList.end(); ++it) {
      std::string itname = it->first.as<std::string>();
      std::vector<double> itrange = it->second.as<std::vector<double>>();
      // Set the mapping as param_name:param_range
      scanRanges[itname] = itrange;
    }
    isScanRanges = true;
  } else {
    MACH3LOG_INFO("There are no user-defined parameter ranges, so I'll use default param bounds for LLH Scans");
  }
  return isScanRanges;
}

// *************************
bool FitterBase::CheckSkipParameter(const std::vector<std::string>& SkipVector, const std::string& ParamName) const {
// *************************
  bool skip = false;
  for(unsigned int is = 0; is < SkipVector.size(); ++is)
  {
    if(ParamName.substr(0, SkipVector[is].length()) == SkipVector[is])
    {
      skip = true;
      break;
    }
  }
  return skip;
}


// *************************
void FitterBase::GetParameterScanRange(const ParameterHandlerBase* cov, const int i, double& CentralValue,
                                       double& lower, double& upper, const int n_points, const std::string& suffix) const {
// *************************
  // YSP: Set up a mapping to store parameters with user-specified ranges, suggested by D. Barrow
  std::map<std::string, std::vector<double>> scanRanges;
  const bool isScanRanges = GetScanRange(scanRanges);

  double nSigma = GetFromManager<int>(fitMan->raw()["LLHScan"]["LLHScanSigma"], 1., __FILE__, __LINE__);
  bool IsPCA = cov->IsPCA();

  // Get the parameter name
  std::string name = cov->GetParFancyName(i);
  if (IsPCA) name += "_PCA";

  // Get the parameter priors and bounds
  CentralValue = cov->GetParProp(i);
  if (IsPCA) CentralValue = cov->GetPCAHandler()->GetParPropPCA(i);

  double prior = cov->GetParInit(i);
  if (IsPCA) prior = cov->GetPCAHandler()->GetPreFitValuePCA(i);

  if (std::abs(CentralValue - prior) > 1e-10) {
    MACH3LOG_INFO("For {} scanning around value {} rather than prior {}", name, CentralValue, prior);
  }

  // Get the covariance matrix and do the +/- nSigma
  // Set lower and upper bounds relative the CentralValue
  // Set the parameter ranges between which LLH points are scanned
  lower = CentralValue - nSigma*cov->GetDiagonalError(i);
  upper = CentralValue + nSigma*cov->GetDiagonalError(i);
  // If PCA, transform these parameter values to the PCA basis
  if (IsPCA) {
    lower = CentralValue - nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(i));
    upper = CentralValue + nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(i));
    MACH3LOG_INFO("eval {} = {:.2f}", i, cov->GetPCAHandler()->GetEigenValues()(i));
    MACH3LOG_INFO("CV {} = {:.2f}", i, CentralValue);
    MACH3LOG_INFO("lower {} = {:.2f}", i, lower);
    MACH3LOG_INFO("upper {} = {:.2f}", i, upper);
    MACH3LOG_INFO("nSigma = {:.2f}", nSigma);
  }

  // Implementation suggested by D. Barrow
  // If param ranges are specified in scanRanges node, extract it from there
  if(isScanRanges){
    // Find matching entries through std::maps
    auto it = scanRanges.find(name);
    if (it != scanRanges.end() && it->second.size() == 2) { //Making sure the range is has only two entries
      lower = it->second[0];
      upper = it->second[1];
      MACH3LOG_INFO("Found matching param name for setting specified range for {}", name);
      MACH3LOG_INFO("Range for {} = [{:.2f}, {:.2f}]", name, lower, upper);
    }
  }

  // Cross-section and flux parameters have boundaries that we scan between, check that these are respected in setting lower and upper variables
  // This also applies for other parameters like osc, etc.
  lower = std::max(lower, cov->GetLowerBound(i));
  upper = std::min(upper, cov->GetUpperBound(i));
  MACH3LOG_INFO("Scanning {} {} with {} steps, from [{:.2f} , {:.2f}], CV = {:.2f}", suffix, name, n_points, lower, upper, CentralValue);
}

// *************************
// Run LLH scan
void FitterBase::RunLLHScan() {
// *************************
  // Save the settings into the output file
  SaveSettings();

  MACH3LOG_INFO("Starting {}", __func__);

  //KS: Turn it on if you want LLH scan for each ND sample separately, which increase time significantly but can be useful for validating new samples or dials.
  bool PlotLLHScanBySample = GetFromManager<bool>(fitMan->raw()["LLHScan"]["LLHScanBySample"], false, __FILE__ , __LINE__);
  auto SkipVector = GetFromManager<std::vector<std::string>>(fitMan->raw()["LLHScan"]["LLHScanSkipVector"], {}, __FILE__ , __LINE__);

  // Now finally get onto the LLH scan stuff
  // Very similar code to MCMC but never start MCMC; just scan over the parameter space
  std::vector<TDirectory *> Cov_LLH(systematics.size());
  for(unsigned int ivc = 0; ivc < systematics.size(); ++ivc )
  {
    std::string NameTemp = systematics[ivc]->GetName();
    NameTemp = NameTemp.substr(0, NameTemp.find("_cov")) + "_LLH";
    Cov_LLH[ivc] = outputFile->mkdir(NameTemp.c_str());
  }

  std::vector<TDirectory *> SampleClass_LLH(samples.size());
  for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
  {
    std::string NameTemp = samples[ivs]->GetName();
    SampleClass_LLH[ivs] = outputFile->mkdir(NameTemp.c_str());
  }

  TDirectory *Sample_LLH = outputFile->mkdir("Sample_LLH");
  TDirectory *Total_LLH = outputFile->mkdir("Total_LLH");

  std::vector<TDirectory *>SampleSplit_LLH;
  if(PlotLLHScanBySample)
  {
    SampleSplit_LLH.resize(TotalNSamples);
    int SampleIterator = 0;
    for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
    {
      for(int is = 0; is < samples[ivs]->GetNsamples(); ++is )
      {
        SampleSplit_LLH[SampleIterator] = outputFile->mkdir((samples[ivs]->GetSampleTitle(is)+ "_LLH").c_str());
        SampleIterator++;
      }
    }
  }
  // Number of points we do for each LLH scan
  const int n_points = GetFromManager<int>(fitMan->raw()["LLHScan"]["LLHScanPoints"], 100, __FILE__ , __LINE__);

  // We print 5 reweights
  const int countwidth = int(double(n_points)/double(5));
   
  // Loop over the covariance classes
  for (ParameterHandlerBase *cov : systematics)
  {
    // Scan over all the parameters
    // Get the number of parameters
    int npars = cov->GetNumParams();
    bool IsPCA = cov->IsPCA();
    if (IsPCA) npars = cov->GetNParameters();
    for (int i = 0; i < npars; ++i)
    {
      // Get the parameter name
      std::string name = cov->GetParFancyName(i);
      if (IsPCA) name += "_PCA";
      // KS: Check if we want to skip this parameter
      if(CheckSkipParameter(SkipVector, name)) continue;
      // Get the parameter central and bounds
      double CentralValue, lower, upper;
      GetParameterScanRange(cov, i, CentralValue, lower, upper, n_points);
      // Make the TH1D
      auto hScan = std::make_unique<TH1D>((name + "_full").c_str(), (name + "_full").c_str(), n_points, lower, upper);
      hScan->SetTitle((std::string("2LLH_full, ") + name + ";" + name + "; -2(ln L_{sample} + ln L_{xsec+flux} + ln L_{det})").c_str());
      hScan->SetDirectory(nullptr);

      auto hScanSam = std::make_unique<TH1D>((name + "_sam").c_str(), (name + "_sam").c_str(), n_points, lower, upper);
      hScanSam->SetTitle((std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());
      hScanSam->SetDirectory(nullptr);

      std::vector<std::unique_ptr<TH1D>> hScanSample(samples.size());
      std::vector<double> nSamLLH(samples.size());
      for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
      {
        std::string NameTemp = samples[ivs]->GetName();
        hScanSample[ivs] = std::make_unique<TH1D>((name+"_"+NameTemp).c_str(), (name+"_" + NameTemp).c_str(), n_points, lower, upper);
        hScanSample[ivs]->SetDirectory(nullptr);
        hScanSample[ivs]->SetTitle(("2LLH_" + NameTemp + ", " + name + ";" + name + "; -2(ln L_{" + NameTemp +"})").c_str());
        nSamLLH[ivs] = 0.;
      }

      std::vector<std::unique_ptr<TH1D>> hScanCov(systematics.size());
      std::vector<double> nCovLLH(systematics.size());
      for(unsigned int ivc = 0; ivc < systematics.size(); ++ivc )
      {
        std::string NameTemp = systematics[ivc]->GetName();
        NameTemp = NameTemp.substr(0, NameTemp.find("_cov"));
        hScanCov[ivc] = std::make_unique<TH1D>((name + "_" + NameTemp).c_str(), (name + "_" + NameTemp).c_str(), n_points, lower, upper);
        hScanCov[ivc]->SetDirectory(nullptr);
        hScanCov[ivc]->SetTitle(("2LLH_" + NameTemp + ", " + name + ";" + name + "; -2(ln L_{" + NameTemp +"})").c_str());
        nCovLLH[ivc] = 0.;
      }

      std::vector<TH1D*> hScanSamSplit;
      std::vector<double> sampleSplitllh;
      if(PlotLLHScanBySample)
      {
        int SampleIterator = 0;
        for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
        {
          hScanSamSplit.resize(TotalNSamples);
          sampleSplitllh.resize(TotalNSamples);
          for(int is = 0; is < samples[ivs]->GetNsamples(); ++is )
          {
            hScanSamSplit[SampleIterator] = new TH1D((name+samples[ivs]->GetSampleTitle(is)).c_str(), (name+samples[ivs]->GetSampleTitle(is)).c_str(), n_points, lower, upper);
            hScanSamSplit[SampleIterator]->SetTitle((std::string("2LLH_sam, ") + name + ";" + name + "; -2(ln L_{sample})").c_str());
            SampleIterator++;
          }
        }
      }

      // Scan over the parameter space
      for (int j = 0; j < n_points; ++j)
      {
        if (j % countwidth == 0)
          MaCh3Utils::PrintProgressBar(j, n_points);

        // For PCA we have to do it differently
        if (IsPCA) {
          cov->GetPCAHandler()->SetParPropPCA(i, hScan->GetBinCenter(j+1));
        } else {
          // Set the parameter
          cov->SetParProp(i, hScan->GetBinCenter(j+1));
        }

        // Reweight the MC
        for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
        {
          samples[ivs]->Reweight();
        }
        //Total LLH
        double totalllh = 0;

        // Get the -log L likelihoods
        double samplellh = 0;

        for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
        {
          nSamLLH[ivs] = samples[ivs]->GetLikelihood();
          samplellh += nSamLLH[ivs];
        }

        for(unsigned int ivc = 0; ivc < systematics.size(); ++ivc )
        {
          nCovLLH[ivc] = systematics[ivc]->GetLikelihood();
          totalllh += nCovLLH[ivc];
        }

        totalllh += samplellh;

        if(PlotLLHScanBySample)
        {
          int SampleIterator = 0;
          for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
          {
            for(int is = 0; is < samples[ivs]->GetNsamples(); ++is)
            {
              sampleSplitllh[SampleIterator] = samples[ivs]->GetSampleLikelihood(is);
              SampleIterator++;
            }
          }
        }

        for(unsigned int ivs = 0; ivs < samples.size(); ++ivs ) {
          hScanSample[ivs]->SetBinContent(j+1, 2*nSamLLH[ivs]);
        }
        for(unsigned int ivc = 0; ivc < systematics.size(); ++ivc ) {
          hScanCov[ivc]->SetBinContent(j+1, 2*nCovLLH[ivc]);
        }

        hScanSam->SetBinContent(j+1, 2*samplellh);
        hScan->SetBinContent(j+1, 2*totalllh);

        if(PlotLLHScanBySample)
        {
          int SampleIterator = 0;
          for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
          {
            for(int is = 0; is < samples[ivs]->GetNsamples(); ++is)
            {
              hScanSamSplit[SampleIterator]->SetBinContent(j+1, 2*sampleSplitllh[SampleIterator]);
              SampleIterator++;
            }
          }
        }
      }
      for(unsigned int ivc = 0; ivc < systematics.size(); ++ivc )
      {
        Cov_LLH[ivc]->cd();
        hScanCov[ivc]->Write();
      }

      for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
      {
        SampleClass_LLH[ivs]->cd();
        hScanSample[ivs]->Write();
      }
      Sample_LLH->cd();
      hScanSam->Write();
      Total_LLH->cd();
      hScan->Write();

      if(PlotLLHScanBySample)
      {
        int SampleIterator = 0;
        for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
        {
          for(int is = 0; is < samples[ivs]->GetNsamples(); ++is)
          {
            SampleSplit_LLH[SampleIterator]->cd();
            hScanSamSplit[SampleIterator]->Write();
            delete hScanSamSplit[SampleIterator];
            SampleIterator++;
          }
        }
      }

      // Reset the parameters to their CentralValue central values
      if (IsPCA) {
        cov->GetPCAHandler()->SetParPropPCA(i, CentralValue);
      } else {
        cov->SetParProp(i, CentralValue);
      }
    }//end loop over systematics
  }//end loop covariance classes

  for(unsigned int ivc = 0; ivc < systematics.size(); ++ivc )
  {
    Cov_LLH[ivc]->Write();
    delete Cov_LLH[ivc];
  }

  for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
  {
    SampleClass_LLH[ivs]->Write();
    delete SampleClass_LLH[ivs];
  }

  Sample_LLH->Write();
  delete Sample_LLH;

  Total_LLH->Write();
  delete Total_LLH;

  if(PlotLLHScanBySample)
  {
    int SampleIterator = 0;
    for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
    {
      for(int is = 0; is < samples[ivs]->GetNsamples(); ++is )
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

  for (ParameterHandlerBase *cov : systematics)
  {
    const int npars = cov->GetNumParams();
    std::vector<double> StepScale(npars);
    for (int i = 0; i < npars; ++i)
    {
      std::string name = cov->GetParFancyName(i);

      StepScale[i] = cov->GetIndivStepScale(i);
      TH1D* LLHScan = Sample_LLH->Get<TH1D>((name+"_sam").c_str());
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
    cov->SetIndivStepScale(StepScale);
    cov->SaveUpdatedMatrixConfig();
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
  auto SkipVector = GetFromManager<std::vector<std::string>>(fitMan->raw()["LLHScan"]["LLHScanSkipVector"], {}, __FILE__ , __LINE__);;

  // Number of points we do for each LLH scan
  const int n_points = GetFromManager<int>(fitMan->raw()["LLHScan"]["2DLLHScanPoints"], 20, __FILE__ , __LINE__);
  // We print 5 reweights
  const int countwidth = int(double(n_points)/double(5));

  // Loop over the covariance classes
  for (ParameterHandlerBase *cov : systematics)
  {
    // Scan over all the parameters
    // Get the number of parameters
    int npars = cov->GetNumParams();
    bool IsPCA = cov->IsPCA();
    if (IsPCA) npars = cov->GetNParameters();

    for (int i = 0; i < npars; ++i)
    {
      std::string name_x = cov->GetParFancyName(i);
      if (IsPCA) name_x += "_PCA";
      // Get the parameter central and bounds
      double central_x, lower_x, upper_x;
      GetParameterScanRange(cov, i, central_x, lower_x, upper_x, n_points, "X");

      // KS: Check if we want to skip this parameter
      if(CheckSkipParameter(SkipVector, name_x)) continue;

      for (int j = 0; j < i; ++j)
      {
        std::string name_y = cov->GetParFancyName(j);
        if (IsPCA) name_y += "_PCA";
        // KS: Check if we want to skip this parameter
        if(CheckSkipParameter(SkipVector, name_y)) continue;

        // Get the parameter central and bounds
        double central_y, lower_y, upper_y;
        GetParameterScanRange(cov, j, central_y, lower_y, upper_y, n_points, "Y");

        auto hScanSam = std::make_unique<TH2D>((name_x + "_" + name_y + "_sam").c_str(), (name_x + "_" + name_y + "_sam").c_str(),
                                                n_points, lower_x, upper_x, n_points, lower_y, upper_y);
        hScanSam->SetDirectory(nullptr);
        hScanSam->GetXaxis()->SetTitle(name_x.c_str());
        hScanSam->GetYaxis()->SetTitle(name_y.c_str());
        hScanSam->GetZaxis()->SetTitle("2LLH_sam");

        // Scan over the parameter space
        for (int x = 0; x < n_points; ++x)
        {
          if (x % countwidth == 0)
            MaCh3Utils::PrintProgressBar(x, n_points);

          for (int y = 0; y < n_points; ++y)
          {
            // For PCA we have to do it differently
            if (IsPCA) {
              cov->GetPCAHandler()->SetParPropPCA(i, hScanSam->GetXaxis()->GetBinCenter(x+1));
              cov->GetPCAHandler()->SetParPropPCA(j, hScanSam->GetYaxis()->GetBinCenter(y+1));
            } else {
              // Set the parameter
              cov->SetParProp(i, hScanSam->GetXaxis()->GetBinCenter(x+1));
              cov->SetParProp(j, hScanSam->GetYaxis()->GetBinCenter(y+1));
            }
            // Reweight the MC
            for(unsigned int ivs = 0; ivs < samples.size(); ++ivs) {
              samples[ivs]->Reweight();
            }

            // Get the -log L likelihoods
            double samplellh = 0;
            for(unsigned int ivs = 0; ivs < samples.size(); ++ivs) {
              samplellh += samples[ivs]->GetLikelihood();
            }
            hScanSam->SetBinContent(x+1, y+1, 2*samplellh);
          }// end loop over y points
        } // end loop over x points

        Sample_2DLLH->cd();
        hScanSam->Write();
        // Reset the parameters to their central central values
        if (IsPCA) {
          cov->GetPCAHandler()->SetParPropPCA(i, central_x);
          cov->GetPCAHandler()->SetParPropPCA(j, central_y);
        } else {
          cov->SetParProp(i, central_x);
          cov->SetParProp(j, central_y);
        }
      } //end loop over systematics y
    }//end loop over systematics X
  }//end loop covariance classes
  Sample_2DLLH->Write();
  delete Sample_2DLLH;
}

// *************************
// Run a general multi-dimensional LLH scan
void FitterBase::RunLLHMap() {
// *************************
  // Save the settings into the output file
  SaveSettings();

  MACH3LOG_INFO("Starting {}", __func__);

  //KS: Turn it on if you want LLH scan for each ND sample separately, which increase time significantly but can be useful for validating new samples or dials.
  bool PlotLLHScanBySample = GetFromManager<bool>(fitMan->raw()["LLHScan"]["LLHScanBySample"], false, __FILE__ , __LINE__);
  auto ParamsOfInterest = GetFromManager<std::vector<std::string>>(fitMan->raw()["LLHScan"]["LLHParameters"], {}, __FILE__, __LINE__);

  if(ParamsOfInterest.empty()) {
    MACH3LOG_WARN("There were no LLH parameters of interest specified to run the LLHMap! LLHMap will not run at all ...");
    return;
  }

  // Parameters IDs within the covariance objects
  // ParamsCovIDs = {Name, CovObj, IDinCovObj}
  std::vector<std::tuple<std::string, ParameterHandlerBase*, int>> ParamsCovIDs;
  for(auto& p : ParamsOfInterest) {
    bool found = false;
    for(auto cov : systematics) {
      for(int c = 0; c < cov->GetNumParams(); ++c) {
        if(cov->GetParName(c) == p || cov->GetParFancyName(c) == p) {
          bool add = true;
          for(auto& pc : ParamsCovIDs) {
            if(std::get<1>(pc) == cov && std::get<2>(pc) == c)
            {
              MACH3LOG_WARN("Parameter {} as {}({}) listed multiple times for LLHMap, omitting and using only once!", p, cov->GetName(), c);
              add = false;
              break;
            }
          }

          if(add)
            ParamsCovIDs.push_back(std::make_tuple(p, cov, c));

          found = true;
          break;
        }
      }
      if(found)
        break;
    }
    if(found)
      MACH3LOG_INFO("Parameter {} found in {} at an index {}.", p, std::get<1>(ParamsCovIDs.back())->GetName(), std::get<2>(ParamsCovIDs.back()));
    else
      MACH3LOG_WARN("Parameter {} not found in any of the systematic covariance objects. Will not scan over this one!", p);
  }

  // ParamsRanges["parameter"] = {nPoints, {low, high}}
  std::map<std::string, std::pair<int, std::pair<double, double>>> ParamsRanges;

  MACH3LOG_INFO("======================================================================================");
  MACH3LOG_INFO("Performing a general multi-dimensional LogL map scan over following parameters ranges:");
  MACH3LOG_INFO("======================================================================================");
  unsigned long TotalPoints = 1;

  double nSigma = GetFromManager<int>(fitMan->raw()["LLHScan"]["LLHScanSigma"], 1., __FILE__, __LINE__);

  // TN: Setting up the scan ranges might look like a re-implementation of the
  // FitterBase::GetScanRange, but I guess the use-case here is a bit different.
  // Anyway, just in case, we can discuss and rewrite to everyone's liking!
  for(auto& p : ParamsCovIDs) {
    // Auxiliary vars to help readability
    std::string name = std::get<0>(p);
    int i = std::get<2>(p);
    ParameterHandlerBase* cov = std::get<1>(p);

    ParamsRanges[name].first = GetFromManager<int>(fitMan->raw()["LLHScan"]["LLHScanPoints"], 20, __FILE__, __LINE__);
    if(CheckNodeExists(fitMan->raw(),"LLHScan","ScanPoints"))
      ParamsRanges[name].first = GetFromManager<int>(fitMan->raw()["LLHScan"]["ScanPoints"][name], ParamsRanges[name].first, __FILE__, __LINE__);

    // Get the parameter priors and bounds
    double CentralValue = cov->GetParProp(i);

    bool IsPCA = cov->IsPCA();
    if (IsPCA)
      CentralValue = cov->GetPCAHandler()->GetParPropPCA(i);

    double prior = cov->GetParInit(i);
    if (IsPCA)
      prior = cov->GetPCAHandler()->GetPreFitValuePCA(i);

    if (std::abs(CentralValue - prior) > 1e-10) {
      MACH3LOG_INFO("For {} scanning around value {} rather than prior {}", name, CentralValue, prior);
    }
    // Get the covariance matrix and do the +/- nSigma
    // Set lower and upper bounds relative the CentralValue
    // Set the parameter ranges between which LLH points are scanned
    double lower = CentralValue - nSigma*cov->GetDiagonalError(i);
    double upper = CentralValue + nSigma*cov->GetDiagonalError(i);
    // If PCA, transform these parameter values to the PCA basis
    if (IsPCA) {
      lower = CentralValue - nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(i));
      upper = CentralValue + nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(i));
      MACH3LOG_INFO("eval {} = {:.2f}", i, cov->GetPCAHandler()->GetEigenValues()(i));
      MACH3LOG_INFO("CV {} = {:.2f}", i, CentralValue);
      MACH3LOG_INFO("lower {} = {:.2f}", i, lower);
      MACH3LOG_INFO("upper {} = {:.2f}", i, upper);
      MACH3LOG_INFO("nSigma = {:.2f}", nSigma);
    }

    ParamsRanges[name].second = {lower,upper};

    if(CheckNodeExists(fitMan->raw(),"LLHScan","ScanRanges"))
      ParamsRanges[name].second = GetFromManager<std::pair<double,double>>(fitMan->raw()["LLHScan"]["ScanRanges"][name], ParamsRanges[name].second, __FILE__, __LINE__);

    MACH3LOG_INFO("{} from {:.4f} to {:.4f} with a {:.5f} step ({} points total)",
                  name, ParamsRanges[name].second.first, ParamsRanges[name].second.second,
                  (ParamsRanges[name].second.second - ParamsRanges[name].second.first)/(ParamsRanges[name].first - 1.),
                  ParamsRanges[name].first);

    TotalPoints *= ParamsRanges[name].first;
  }

  // TN: Waiting for C++ 20 std::format() function
  MACH3LOG_INFO("In total, looping over {} points, from {} parameters. Estimates for run time:", TotalPoints, ParamsCovIDs.size());
  MACH3LOG_INFO("   1 s per point = {} hours", double(TotalPoints)/3600.);
  MACH3LOG_INFO(" 0.1 s per point = {} hours", double(TotalPoints)/36000.);
  MACH3LOG_INFO("0.01 s per point = {} hours", double(TotalPoints)/360000.);
  MACH3LOG_INFO("==================================================================================");

  const int countwidth = int(double(TotalPoints)/double(20));

  // Tree to store LogL values
  auto LLHMap = new TTree("llhmap", "LLH Map");

  std::vector<double> CovLogL(systematics.size());
  for(unsigned int ivc = 0; ivc < systematics.size(); ++ivc )
  {
    std::string NameTemp = systematics[ivc]->GetName();
    NameTemp = NameTemp.substr(0, NameTemp.find("_cov")) + "_LLH";
    LLHMap->Branch(NameTemp.c_str(), &CovLogL[ivc]);
  }

  std::vector<double> SampleClassLogL(samples.size());
  for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
  {
    std::string NameTemp = samples[ivs]->GetName()+"_LLH";
    LLHMap->Branch(NameTemp.c_str(), &SampleClassLogL[ivs]);
  }

  double SampleLogL, TotalLogL;
  LLHMap->Branch("Sample_LLH", &SampleLogL);
  LLHMap->Branch("Total_LLH", &TotalLogL);

  std::vector<double>SampleSplitLogL;
  if(PlotLLHScanBySample)
  {
    SampleSplitLogL.resize(TotalNSamples);
    int SampleIterator = 0;
    for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
    {
      for(int is = 0; is < samples[ivs]->GetNsamples(); ++is )
      {
        std::string NameTemp =samples[ivs]->GetSampleTitle(is)+"_LLH";
        LLHMap->Branch(NameTemp.c_str(), &SampleSplitLogL[SampleIterator]);
        SampleIterator++;
      }
    }
  }

  std::vector<double> ParamsValues(ParamsCovIDs.size());
  for(unsigned int i=0; i < ParamsCovIDs.size(); ++i)
    LLHMap->Branch(std::get<0>(ParamsCovIDs[i]).c_str(), &ParamsValues[i]);

  // Setting up the scan
  // Starting at index {0,0,0,...}
  std::vector<unsigned long> idx(ParamsCovIDs.size(), 0);

  // loop over scanned points sp
  for(unsigned long sp = 0; sp < TotalPoints; ++sp)
  {
    // At each point need to find the indeces and test values to calculate LogL
    for(unsigned int n = 0; n < ParamsCovIDs.size(); ++n)
    {
      // Auxiliaries
      std::string name = std::get<0>(ParamsCovIDs[n]);
      int points  = ParamsRanges[name].first;
      double low  = ParamsRanges[name].second.first;
      double high = ParamsRanges[name].second.second;

      // Find the n-th index of the sp-th scanned point
      unsigned long dev = 1;
      for(unsigned int m = 0; m <= n; ++m)
        dev *= ParamsRanges[std::get<0>(ParamsCovIDs[m])].first;

      idx[n] = sp % dev;
      if (n > 0)
        idx[n] = idx[n] / ( dev / points );

      // Parameter test value = low + ( high - low ) * idx / ( #points - 1 )
      ParamsValues[n] = low + (high-low) * double(idx[n])/double(points-1);

      // Now set the covariance objects
      // Auxiliary
      ParameterHandlerBase* cov = std::get<1>(ParamsCovIDs[n]);
      int i = std::get<2>(ParamsCovIDs[n]);

      if(cov->IsPCA())
        cov->GetPCAHandler()->SetParPropPCA(i, ParamsValues[n]);
      else
        cov->SetParProp(i, ParamsValues[n]);
    }

    // Reweight samples and calculate LogL
    TotalLogL = .0;
    SampleLogL = .0;

    for(unsigned int ivs = 0; ivs < samples.size(); ++ivs)
      samples[ivs]->Reweight();

    for(unsigned int ivs = 0; ivs < samples.size(); ++ivs)
    {
      SampleClassLogL[ivs] = 2.*samples[ivs]->GetLikelihood();
      SampleLogL += SampleClassLogL[ivs];
    }
    TotalLogL += SampleLogL;

    // CovObjs LogL
    for(unsigned int ivc = 0; ivc < systematics.size(); ++ivc )
    {
      CovLogL[ivc] = 2.*systematics[ivc]->GetLikelihood();
      TotalLogL += CovLogL[ivc];
    }

    if(PlotLLHScanBySample)
    {
      int SampleIterator = 0;
      for(unsigned int ivs = 0; ivs < samples.size(); ++ivs )
      {
        for(int is = 0; is < samples[ivs]->GetNsamples(); ++is)
        {
          SampleSplitLogL[SampleIterator] = 2.*samples[ivs]->GetSampleLikelihood(is);
          SampleIterator++;
        }
      }
    }

    LLHMap->Fill();

    if (sp % countwidth == 0)
      MaCh3Utils::PrintProgressBar(sp, TotalPoints);
  }

  outputFile->cd();
  LLHMap->Write();
}

// *************************
// For comparison with P-Theta we usually have to apply different parameter values then usual 1, 3 sigma
void FitterBase::CustomRange(const std::string& ParName, const double sigma, double& ParamShiftValue) const {
// *************************
  if(!fitMan->raw()["SigmaVar"]["CustomRange"]) return;

  auto Config = fitMan->raw()["SigmaVar"]["CustomRange"];

  const auto sigmaStr = std::to_string(static_cast<int>(std::round(sigma)));

  if (Config[ParName] && Config[ParName][sigmaStr]) {
    ParamShiftValue = Config[ParName][sigmaStr].as<double>();
    MACH3LOG_INFO("  ::: setting custom range from config ::: {} -> {}", ParName, ParamShiftValue);
  }
}

// *************************
/// Helper to write histograms
void WriteHistograms(TH1 *hist, const std::string& baseName) {
// *************************
  if (!hist) return;
  hist->SetTitle(baseName.c_str());
  // Get the class name of the histogram
  TString className = hist->ClassName();

  // Set the appropriate axis title based on the histogram type
  if (className.Contains("TH1")) {
    hist->GetYaxis()->SetTitle("Events");
  } else if (className.Contains("TH2")) {
    hist->GetZaxis()->SetTitle("Events");
  }
  hist->SetDirectory(nullptr);
  hist->Write(baseName.c_str());
}

// *************************
/// Generic histogram writer - should make main code more palatable
void WriteHistogramsByMode(SampleHandlerBase *sample,
                           const std::string& suffix,
                           const bool by_mode,
                           const bool by_channel,
                           const std::vector<TDirectory*>& SampleDir) {
// *************************
  MaCh3Modes *modes = sample->GetMaCh3Modes();
  for (int iSample = 0; iSample < sample->GetNsamples(); ++iSample) {
    SampleDir[iSample]->cd();
    const std::string sampleName = sample->GetSampleTitle(iSample);
    for(int iDim1 = 0; iDim1 < sample->GetNDim(iSample); iDim1++) {
      std::string ProjectionName = sample->GetKinVarName(iSample, iDim1);
      std::string ProjectionSuffix = "_1DProj" + std::to_string(iDim1);

      // Probably a better way of handling this logic
      if (by_mode) {
        for (int iMode = 0; iMode < modes->GetNModes(); ++iMode) {
          auto modeHist = sample->Get1DVarHistByModeAndChannel(iSample, ProjectionName, iMode);
          WriteHistograms(modeHist, sampleName + "_" + modes->GetMaCh3ModeName(iMode) + ProjectionSuffix + suffix);
          delete modeHist;
        }
      }

      if (by_channel) {
        for (int iChan = 0; iChan < sample->GetNOscChannels(iSample); ++iChan) {
          auto chanHist = sample->Get1DVarHistByModeAndChannel(iSample, ProjectionName, -1, iChan); // -1 skips over mode plotting
          WriteHistograms(chanHist, sampleName + "_" + sample->GetFlavourName(iSample, iChan) + ProjectionSuffix + suffix);
          delete chanHist;
        }
      }

      if (by_mode && by_channel) {
        for (int iMode = 0; iMode < modes->GetNModes(); ++iMode) {
          for (int iChan = 0; iChan < sample->GetNOscChannels(iSample); ++iChan) {
            auto hist = sample->Get1DVarHistByModeAndChannel(iSample, ProjectionName, iMode, iChan);
            WriteHistograms(hist, sampleName + "_" + modes->GetMaCh3ModeName(iMode) + "_" + sample->GetFlavourName(iSample, iChan) + ProjectionSuffix + suffix);
            delete hist;
          }
        }
      }

      if (!by_mode && !by_channel) {
        auto hist = sample->Get1DVarHist(iSample, ProjectionName);
        WriteHistograms(hist, sampleName + ProjectionSuffix + suffix);
        delete hist;
        // Only for 2D and Beyond
        for (int iDim2 = iDim1 + 1; iDim2 < sample->GetNDim(iSample); ++iDim2) {
          // Get the names for the two dimensions
          std::string XVarName = sample->GetKinVarName(iSample, iDim1);
          std::string YVarName = sample->GetKinVarName(iSample, iDim2);

          // Get the 2D histogram for this pair
          auto hist2D = sample->Get2DVarHist(iSample, XVarName, YVarName);

          // Write the histogram
          std::string suffix2D = "_2DProj_" + std::to_string(iDim1) + "_vs_" + std::to_string(iDim2) + suffix;
          WriteHistograms(hist2D, sampleName + suffix2D);

          // Clean up
          delete hist2D;
        }
      }
    }
  }
}

// *************************
void FitterBase::RunSigmaVar() {
// *************************
  // Save the settings into the output file
  SaveSettings();

  bool plot_by_mode = GetFromManager<bool>(fitMan->raw()["SigmaVar"]["PlotByMode"], false);
  bool plot_by_channel = GetFromManager<bool>(fitMan->raw()["SigmaVar"]["PlotByChannel"], false);
  auto SkipVector = GetFromManager<std::vector<std::string>>(fitMan->raw()["SigmaVar"]["SkipVector"], {}, __FILE__ , __LINE__);

  if (plot_by_mode) MACH3LOG_INFO("Plotting by sample and mode");
  if (plot_by_channel) MACH3LOG_INFO("Plotting by sample and channel");
  if (!plot_by_mode && !plot_by_channel) MACH3LOG_INFO("Plotting by sample only");
  if (plot_by_mode && plot_by_channel) MACH3LOG_INFO("Plotting by sample, mode and channel");

  auto SigmaArray = GetFromManager<std::vector<double>>(fitMan->raw()["SigmaVar"]["SigmaArray"], {-3, -1, 0, 1, 3}, __FILE__ , __LINE__);
  if (std::find(SigmaArray.begin(), SigmaArray.end(), 0.0) == SigmaArray.end()) {
    MACH3LOG_ERROR(":: SigmaArray does not contain 0! Current contents: {} ::", fmt::join(SigmaArray, ", "));
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TDirectory* SigmaDir = outputFile->mkdir("SigmaVar");
  outputFile->cd();

  for (size_t s = 0; s < systematics.size(); ++s)
  {
    for(int i = 0; i < systematics[s]->GetNumParams(); i++)
    {
      std::string ParName = systematics[s]->GetParFancyName(i);
      // KS: Check if we want to skip this parameter
      if(CheckSkipParameter(SkipVector, ParName)) continue;

      MACH3LOG_INFO(":: Param {} ::", systematics[s]->GetParFancyName(i));

      TDirectory* ParamDir = SigmaDir->mkdir(ParName.c_str());
      ParamDir->cd();

      const double ParamCentralValue = systematics[s]->GetParProp(i);
      const double Prior = systematics[s]->GetParInit(i);
      const double ParamLower = systematics[s]->GetLowerBound(i);
      const double ParamUpper = systematics[s]->GetUpperBound(i);

      if (std::abs(ParamCentralValue - Prior) > 1e-10) {
        MACH3LOG_INFO("For {} scanning around value {} rather than prior {}", ParName, ParamCentralValue, Prior);
      }

      for(unsigned int iSample = 0; iSample < samples.size(); ++iSample)
      {
        auto* MaCh3Sample = samples[iSample];
        std::vector<TDirectory*> SampleDir(MaCh3Sample->GetNsamples());
        for (int SampleIndex = 0; SampleIndex < MaCh3Sample->GetNsamples(); ++SampleIndex) {
          SampleDir[SampleIndex] = ParamDir->mkdir(MaCh3Sample->GetSampleTitle(SampleIndex).c_str());
        }

        for (size_t j = 0; j < SigmaArray.size(); ++j) {
          double sigma = SigmaArray[j];

          double ParamShiftValue = ParamCentralValue + sigma * std::sqrt((*systematics[s]->GetCovMatrix())(i,i));
          ParamShiftValue = std::max(std::min(ParamShiftValue, ParamUpper), ParamLower);

          /// Apply custom range to make easier comparison with p-theta
          CustomRange(ParName, sigma, ParamShiftValue);

          MACH3LOG_INFO("  - set to {:<5.2f} ({:<2} sigma shift)", ParamShiftValue, sigma);
          systematics[s]->SetParProp(i, ParamShiftValue);

          std::ostringstream valStream;
          valStream << std::fixed << std::setprecision(2) << ParamShiftValue;
          std::string valueStr = valStream.str();

          std::ostringstream sigmaStream;
          sigmaStream << std::fixed << std::setprecision(2) << std::abs(sigma);
          std::string sigmaStr = sigmaStream.str();

          std::string suffix;
          if (sigma == 0) {
            suffix = "_" + ParName + "_nom_val_" + valueStr;
          } else {
            std::string sign = (sigma > 0) ? "p" : "n";
            suffix = "_" + ParName + "_sig_" + sign + sigmaStr + "_val_" + valueStr;
          }

          systematics[s]->SetParProp(i, ParamShiftValue);
          MaCh3Sample->Reweight();

          WriteHistogramsByMode(MaCh3Sample, suffix, plot_by_mode, plot_by_channel, SampleDir);
        }
        for (int subSampleIndex = 0; subSampleIndex < MaCh3Sample->GetNsamples(); ++subSampleIndex) {
          SampleDir[subSampleIndex]->Close();
          delete SampleDir[subSampleIndex];
        }
        ParamDir->cd();
      }

      systematics[s]->SetParProp(i, ParamCentralValue);
      MACH3LOG_INFO("  - set back to CV {:<5.2f}", ParamCentralValue);
      MACH3LOG_INFO("");
      ParamDir->Close();
      delete ParamDir;
      SigmaDir->cd();
    } // end loop over params
  } // end loop over systemics
  SigmaDir->Close();
  delete SigmaDir;

  outputFile->cd();
}
