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
// Initialise the manager and make it an object of FitterBase class
// Now we can dump manager settings to the output file
FitterBase::FitterBase(manager * const man) : fitMan(man) {
// *************************
  //Get mach3 modes from manager
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
  for(unsigned int i = 0; i < samples.size(); ++i)
    MACH3LOG_INFO("{}: SampleHandler name: {}, it has {} samples, {} OscChannels",i , samples[i]->GetTitle(), samples[i]->GetNsamples(), samples[i]->GetNOscChannels());

  //TN: Have to close the folder in order to write it to disk before SaveOutput is called in the destructor
  CovFolder->Close();
  
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
    for (ParameterHandlerBase *cov : systematics) {
      cov->SetBranches(*outTree, SaveProposal);
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

  outputFile->cd();
  outTree->Write();

  MACH3LOG_INFO("{} steps took {:.2f} seconds to complete. ({:.2f}s / step).", step - stepStart, clock->RealTime(), clock->RealTime() / static_cast<double>(step - stepStart));
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
// Add SampleHandler object to the Markov Chain
void FitterBase::AddSampleHandler(SampleHandlerBase * const sample) {
// *************************
  //Check if the sample has a unique name
  for (const auto &s : samples) {
    if (s->GetTitle() == sample->GetTitle()) {
      MACH3LOG_ERROR("SampleHandler with name '{}' already exists!", sample->GetTitle());
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }

  TotalNSamples += sample->GetNsamples();
  MACH3LOG_INFO("Adding {} object, with {} samples", sample->GetTitle(), sample->GetNsamples());
  samples.push_back(sample);
}

// *************************
// Add flux systematics, cross-section systematics, ND systematics to the chain
void FitterBase::AddSystObj(ParameterHandlerBase * const cov) {
// *************************
  MACH3LOG_INFO("Adding systematic object {}, with {} params", cov->GetName(), cov->GetNumParams());
  systematics.push_back(cov);

  CovFolder->cd();
  std::vector<double> n_vec(cov->GetNumParams());
  for (int i = 0; i < cov->GetNumParams(); ++i)
    n_vec[i] = cov->GetParInit(i);

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

  TFile *infile = new TFile(FitName.c_str(), "READ");
  TTree *posts = infile->Get<TTree>("posteriors");
  int step_val = 0;
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
        MACH3LOG_ERROR("Yaml configs in previous chain and current one are different", FitName);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }

      delete ConfigCov;
    }

    CovarianceFolder->Close();
    delete CovarianceFolder;

    std::vector<double> branch_vals(systematics[s]->GetNumParams(), M3::_BAD_DOUBLE_);
    for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
      posts->SetBranchAddress(systematics[s]->GetParName(i).c_str(), &branch_vals[i]);
    }
    posts->GetEntry(posts->GetEntries()-1);

    for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
      if(branch_vals[i] == M3::_BAD_DOUBLE_)
      {
        MACH3LOG_ERROR("Parameter {} is unvitalised with value {}", i, branch_vals[i]);
        MACH3LOG_ERROR("Please check more precisely chain you passed {}", FitName);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
    }
    systematics[s]->SetParameters(branch_vals);
    systematics[s]->AcceptStep();

    MACH3LOG_INFO("Printing new starting values for: {}", systematics[s]->GetName());
    systematics[s]->PrintNominalCurrProp();

    // Resetting branch adressed to nullptr as we don't want to write into a delected vector out of scope...
    for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
      posts->SetBranchAddress(systematics[s]->GetParName(i).c_str(), nullptr);
    }
  }
  logLCurr = log_val;

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
      outputFile = new TFile(Get<std::string>(fitMan->raw()["General"]["Output"]["Filename"], __FILE__, __LINE__).c_str(), "UPDATE");
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
    MACH3LOG_INFO("It took {:.4f} s to reweights {} times sample: {}", clockRace.RealTime(), NLaps, samples[ivs]->GetTitle());
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
    MACH3LOG_INFO("It took {:.4f} s to calculate  GetLikelihood {} times sample:  {}", clockRace.RealTime(), NLaps, samples[ivs]->GetTitle());
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
bool FitterBase::GetScaneRange(std::map<std::string, std::vector<double>>& scanRanges) {
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
// Run LLH scan
void FitterBase::RunLLHScan() {
// *************************
  // Save the settings into the output file
  SaveSettings();

  MACH3LOG_INFO("Starting LLH Scan");

  //KS: Turn it on if you want LLH scan for each ND sample separately, which increase time significantly but can be useful for validating new samples or dials.
  bool PlotLLHScanBySample = GetFromManager<bool>(fitMan->raw()["LLHScan"]["LLHScanBySample"], false, __FILE__ , __LINE__);;
  auto SkipVector = GetFromManager<std::vector<std::string>>(fitMan->raw()["LLHScan"]["LLHScanSkipVector"], {}, __FILE__ , __LINE__);;

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
    std::string NameTemp = samples[ivs]->GetTitle();
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
        SampleSplit_LLH[SampleIterator] = outputFile->mkdir((samples[ivs]->GetSampleName(is)+ "_LLH").c_str());
        SampleIterator++;
      }
    }
  }
  // Number of points we do for each LLH scan
  const int n_points = GetFromManager<int>(fitMan->raw()["LLHScan"]["LLHScanPoints"], 100, __FILE__ , __LINE__);
  double nSigma = GetFromManager<int>(fitMan->raw()["LLHScan"]["LLHScanSigma"], 1., __FILE__, __LINE__);

  // We print 5 reweights
  const int countwidth = int(double(n_points)/double(5));

  // YSP: Set up a mapping to store parameters with user-specified ranges, suggested by D. Barrow
  std::map<std::string, std::vector<double>> scanRanges;
  const bool isScanRanges = GetScaneRange(scanRanges);
   
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

      // Get the parameter priors and bounds
      double prior = cov->GetParInit(i);
      if (IsPCA) prior = cov->GetPCAHandler()->GetParCurrPCA(i);

      // Get the covariance matrix and do the +/- nSigma
      // Set lower and upper bounds relative the prior
      // Set the parameter ranges between which LLH points are scanned
      double lower = prior - nSigma*cov->GetDiagonalError(i);
      double upper = prior + nSigma*cov->GetDiagonalError(i);
      // If PCA, transform these parameter values to the PCA basis
      if (IsPCA) {
        lower = prior - nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(i));
        upper = prior + nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(i));
        MACH3LOG_INFO("eval {} = {:.2f}", i, cov->GetPCAHandler()->GetEigenValues()(i));
        MACH3LOG_INFO("prior {} = {:.2f}", i, prior);
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
      MACH3LOG_INFO("Scanning {} with {} steps, from [{:.2f} , {:.2f}], prior = {:.2f}", name, n_points, lower, upper, prior);

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
        std::string NameTemp = samples[ivs]->GetTitle();
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
            hScanSamSplit[SampleIterator] = new TH1D((name+samples[ivs]->GetSampleName(is)).c_str(), (name+samples[ivs]->GetSampleName(is)).c_str(), n_points, lower, upper);
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

      // Reset the parameters to their prior central values
      if (IsPCA) {
        cov->GetPCAHandler()->SetParPropPCA(i, prior);
      } else {
        cov->SetParProp(i, prior);
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

  std::map<std::string, std::vector<double>> scanRanges;
  const bool isScanRanges = GetScaneRange(scanRanges);

  const double nSigma = GetFromManager<int>(fitMan->raw()["LLHScan"]["LLHScanSigma"], 1., __FILE__, __LINE__);

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

      // Get the parameter priors and bounds
      double prior_x = cov->GetParInit(i);
      if (IsPCA) prior_x = cov->GetPCAHandler()->GetParCurrPCA(i);

      // Get the covariance matrix and do the +/- nSigma
      // Set lower and upper bounds relative the prior
      double lower_x = prior_x - nSigma*cov->GetDiagonalError(i);
      double upper_x = prior_x + nSigma*cov->GetDiagonalError(i);
      // If PCA, transform these parameter values to the PCA basis
      if (IsPCA) {
        lower_x = prior_x - nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(i));
        upper_x = prior_x + nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(i));
        MACH3LOG_INFO("eval {} = {:.2f}", i, cov->GetPCAHandler()->GetEigenValues()(i));
        MACH3LOG_INFO("prior {} = {:.2f}", i, prior_x);
        MACH3LOG_INFO("lower {} = {:.2f}", i, lower_x);
        MACH3LOG_INFO("upper {} = {:.2f}", i, upper_x);
        MACH3LOG_INFO("nSigma = {:.2f}", nSigma);
      }
      // If param ranges are specified in scanRanges node, extract it from there
      if(isScanRanges){
        // Find matching entries through std::maps
        auto it = scanRanges.find(name_x);
        if (it != scanRanges.end() && it->second.size() == 2) { //Making sure the range is has only two entries
          lower_x = it->second[0];
          upper_x = it->second[1];
          MACH3LOG_INFO("Found matching param name for setting specified range for {}", name_x);
          MACH3LOG_INFO("Range for {} = [{:.2f}, {:.2f}]", name_x, lower_x, upper_x);
        }
      }

      // Cross-section and flux parameters have boundaries that we scan between, check that these are respected in setting lower and upper variables
      lower_x = std::max(lower_x, cov->GetLowerBound(i));
      upper_x = std::min(upper_x, cov->GetUpperBound(i));
      // KS: Check if we want to skip this parameter
      if(CheckSkipParameter(SkipVector, name_x)) continue;

      for (int j = 0; j < i; ++j)
      {
        std::string name_y = cov->GetParFancyName(j);
        if (IsPCA) name_y += "_PCA";
        // KS: Check if we want to skip this parameter
        if(CheckSkipParameter(SkipVector, name_y)) continue;

        // Get the parameter priors and bounds
        double prior_y = cov->GetParInit(j);
        if (IsPCA) prior_y = cov->GetPCAHandler()->GetParCurrPCA(j);

        // Set lower and upper bounds relative the prior
        double lower_y = prior_y - nSigma*cov->GetDiagonalError(j);
        double upper_y = prior_y + nSigma*cov->GetDiagonalError(j);
        // If PCA, transform these parameter values to the PCA basis
        if (IsPCA) {
          lower_y = prior_y - nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(j));
          upper_y = prior_y + nSigma*std::sqrt((cov->GetPCAHandler()->GetEigenValues())(j));
          MACH3LOG_INFO("eval {} = {:.2f}", i, cov->GetPCAHandler()->GetEigenValues()(j));
          MACH3LOG_INFO("prior {} = {:.2f}", i, prior_y);
          MACH3LOG_INFO("lower {} = {:.2f}", i, lower_y);
          MACH3LOG_INFO("upper {} = {:.2f}", i, upper_y);
          MACH3LOG_INFO("nSigma = {:.2f}", nSigma);
        }
        // If param ranges are specified in scanRanges node, extract it from there
        if(isScanRanges){
          // Find matching entries through std::maps
          auto it = scanRanges.find(name_y);
          if (it != scanRanges.end() && it->second.size() == 2) { //Making sure the range is has only two entries
            lower_y = it->second[0];
            upper_y = it->second[1];
            MACH3LOG_INFO("Found matching param name for setting specified range for {}", name_y);
            MACH3LOG_INFO("Range for {} = [{:.2f}, {:.2f}]", name_y, lower_y, upper_y);
          }
        }

        // Cross-section and flux parameters have boundaries that we scan between, check that these are respected in setting lower and upper variables
        lower_y = std::max(lower_y, cov->GetLowerBound(j));
        upper_y = std::min(upper_y, cov->GetUpperBound(j));
        MACH3LOG_INFO("Scanning X {} with {} steps, from {:.2f} - {:.2f}, prior = {}", name_x, n_points, lower_x, upper_x, prior_x);
        MACH3LOG_INFO("Scanning Y {} with {} steps, from {:.2f} - {:.2f}, prior = {}", name_y, n_points, lower_y, upper_y, prior_y);

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
        // Reset the parameters to their prior central values
        if (IsPCA) {
          cov->GetPCAHandler()->SetParPropPCA(i, prior_x);
          cov->GetPCAHandler()->SetParPropPCA(j, prior_y);
        } else {
          cov->SetParProp(i, prior_x);
          cov->SetParProp(j, prior_y);
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

  bool DoByMode = GetFromManager<int>(fitMan->raw()["General"]["DoByMode"], false, __FILE__ , __LINE__);

  //KS: If true it will make additional plots with LLH sample contribution in each bin, should make it via config file...
  bool PlotLLHperBin = false;
  auto SkipVector = GetFromManager<std::vector<std::string>>(fitMan->raw()["LLHScan"]["LLHScanSkipVector"], {}, __FILE__ , __LINE__);;

  for (ParameterHandlerBase *cov : systematics)
  {
    TMatrixDSym *Cov = cov->GetCovMatrix();

    if(cov->IsPCA())
    {
      MACH3LOG_ERROR("Using PCAed matrix not implemented within sigma var code, I am sorry :(");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    // Loop over xsec parameters
    for (int i = 0; i < cov->GetNumParams(); ++i)
    {
      // Get the parameter name
      std::string name = cov->GetParFancyName(i);
      // KS: Check if we want to skip this parameter
      if(CheckSkipParameter(SkipVector, name)) continue;

      outputFile->cd();
      TDirectory* dirArryDial = outputFile->mkdir(name.c_str());
      std::vector<TDirectory*> dirArrySample(TotalNSamples);

      int SampleIterator = 0;
      // Get each sample and how it's responded to our reweighted parameter
      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        for(int k = 0; k < samples[ivs]->GetNsamples(); k++ )
        {
          std::string title = std::string(samples[ivs]->GetPDF(k)->GetName());
          dirArryDial->cd();
          dirArrySample[SampleIterator] = dirArryDial->mkdir(title.c_str());
          SampleIterator++;
        }
      }

      // Get the initial value of ith parameter
      double init = cov->GetParInit(i);

      std::vector<std::vector<std::unique_ptr<TH1D>>> sigmaArray_x(numVar);
      std::vector<std::vector<std::unique_ptr<TH1D>>> sigmaArray_y(numVar);
      std::vector<std::vector<std::unique_ptr<TH1D>>> sigmaArray_x_norm(numVar);
      std::vector<std::vector<std::unique_ptr<TH1D>>> sigmaArray_y_norm(numVar);

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
        double paramVal = cov->GetParInit(i)+sigmaArray[j]*std::sqrt((*Cov)(i,i));

        // Check the bounds on the parameter
        paramVal = std::max(cov->GetLowerBound(i), std::min(paramVal, cov->GetUpperBound(i)));

        // Set the parameter
        cov->SetParProp(i, paramVal);
        // And reweight the sample
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++) {
          samples[ivs]->Reweight();
        }

        sigmaArray_x[j].resize(TotalNSamples);
        sigmaArray_y[j].resize(TotalNSamples);
        sigmaArray_x_norm[j].resize(TotalNSamples);
        sigmaArray_y_norm[j].resize(TotalNSamples);

        if (DoByMode)
        {
          sigmaArray_mode_x[j] = new TH1D**[TotalNSamples]();
          sigmaArray_mode_y[j] = new TH1D**[TotalNSamples]();
        }

        SampleIterator = 0;
        // Get each sample and how it's responded to our reweighted parameter
        for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
        {
          for (int k = 0; k < samples[ivs]->GetNsamples(); ++k)
          {
            // Make a string of the double
            std::ostringstream ss;
            ss << paramVal;
            std::string parVarTitle = name + "_" + ss.str();

            auto currSamp = M3::Clone<TH2Poly>(static_cast<TH2Poly*>(samples[ivs]->GetPDF(k)));
            // Set a descriptiv-ish title
            std::string title_long = std::string(currSamp->GetName())+"_"+parVarTitle;

            // Enable the mode histograms AFTER addSelection is called
            //Get the 1d binning we want. Let's just use SetupBinning to get this as it already exists
            std::vector<double> xbins;
            std::vector<double> ybins;
            samples[ivs]->SetupBinning(M3::int_t(k), xbins, ybins);

            //KS:here we loop over all reaction modes defined in "RelevantModes[nRelevantModes]"
            if (DoByMode)
            {
              //KS: this is only relevant if PlotByMode is turned on
              //Checking each mode is time consuming so we only consider one which are relevant for particular analysis
              MaCh3Modes_t RelevantModes[nRelevantModes] = {samples[ivs]->GetMaCh3Modes()->GetMode("CCQE"), samples[ivs]->GetMaCh3Modes()->GetMode("2p2h")};

              sigmaArray_mode_x[j][SampleIterator] = new TH1D*[nRelevantModes]();
              sigmaArray_mode_y[j][SampleIterator] = new TH1D*[nRelevantModes]();
              // Now get the TH2D mode variations
              std::string mode_title_long;

              for(int ir = 0; ir < nRelevantModes; ir++)
              {
                auto currSampMode = M3::Clone<TH2Poly>(static_cast<TH2Poly*>(samples[ivs]->GetPDFMode(k, RelevantModes[ir])));

                mode_title_long = title_long + "_" + samples[ivs]->GetMaCh3Modes()->GetMaCh3ModeName(RelevantModes[ir]);
                currSampMode->SetNameTitle(mode_title_long.c_str(), mode_title_long.c_str());
                dirArrySample[SampleIterator]->cd();
                currSampMode->Write();

                sigmaArray_mode_x[j][SampleIterator][ir] = PolyProjectionX(currSampMode.get(), (mode_title_long+"_xProj").c_str(), xbins);
                sigmaArray_mode_x[j][SampleIterator][ir]->SetDirectory(nullptr);
                sigmaArray_mode_y[j][SampleIterator][ir] = PolyProjectionY(currSampMode.get(), (mode_title_long+"_yProj").c_str(), ybins);
                sigmaArray_mode_y[j][SampleIterator][ir]->SetDirectory(nullptr);
              }
            }

            //KS: This will give different results depending if data or Asimov, both have their uses
            if (PlotLLHperBin)
            {
              auto currLLHSamp = M3::Clone<TH2Poly>(static_cast<TH2Poly*>(samples[ivs]->GetPDF(k)));
              currLLHSamp->Reset("");
              currLLHSamp->Fill(0.0, 0.0, 0.0);

              TH2Poly* MCpdf = static_cast<TH2Poly*>(samples[ivs]->GetPDF(k));
              TH2Poly* Datapdf = static_cast<TH2Poly*>(samples[ivs]->GetData(k));
              TH2Poly* W2pdf = samples[ivs]->GetW2(k);

              for(int bin = 1; bin < currLLHSamp->GetNumberOfBins()+1; bin++)
              {
                const double mc = MCpdf->GetBinContent(bin);
                const double dat = Datapdf->GetBinContent(bin);
                const double w2 = W2pdf->GetBinContent(bin);
                currLLHSamp->SetBinContent(bin, samples[ivs]->GetTestStatLLH(dat, mc, w2));
              }
              currLLHSamp->SetNameTitle((title_long+"_LLH").c_str() ,(title_long+"_LLH").c_str());
              dirArrySample[SampleIterator]->cd();
              currLLHSamp->Write();
            }

            // Project down onto x axis
            sigmaArray_x[j][SampleIterator] = std::unique_ptr<TH1D>(PolyProjectionX(currSamp.get(), (title_long+"_xProj").c_str(), xbins));
            sigmaArray_x[j][SampleIterator]->SetDirectory(nullptr);
            sigmaArray_x[j][SampleIterator]->GetXaxis()->SetTitle(currSamp->GetXaxis()->GetTitle());
            sigmaArray_y[j][SampleIterator] = std::unique_ptr<TH1D>(PolyProjectionY(currSamp.get(), (title_long+"_yProj").c_str(), ybins));
            sigmaArray_y[j][SampleIterator]->SetDirectory(nullptr);
            sigmaArray_y[j][SampleIterator]->GetXaxis()->SetTitle(currSamp->GetYaxis()->GetTitle());

            sigmaArray_x_norm[j][SampleIterator] = M3::Clone<TH1D>(sigmaArray_x[j][SampleIterator].get());
            sigmaArray_x_norm[j][SampleIterator]->Scale(1., "width");
            sigmaArray_y_norm[j][SampleIterator] = M3::Clone<TH1D>(sigmaArray_y[j][SampleIterator].get());
            sigmaArray_y_norm[j][SampleIterator]->SetDirectory(nullptr);
            sigmaArray_y_norm[j][SampleIterator]->Scale(1., "width");

            currSamp->SetNameTitle(title_long.c_str(), title_long.c_str());
            dirArrySample[k]->cd();
            currSamp->Write();

            sigmaArray_x[j][k]->Write();
            sigmaArray_y[j][k]->Write();
            SampleIterator++;
          }//End loop over samples
        }
      } // End looping over variation

      // Restore the parameter to prior value
      cov->SetParProp(i, init);

      SampleIterator = 0;
      // Get each sample and how it's responded to our reweighted parameter
      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        for (int k = 0; k < samples[ivs]->GetNsamples(); ++k)
        {
          std::string title = std::string(samples[ivs]->GetPDF(k)->GetName()) + "_" + name;
          auto var_x = MakeAsymGraph(sigmaArray_x[1][SampleIterator].get(), sigmaArray_x[2][SampleIterator].get(), sigmaArray_x[3][SampleIterator].get(), (title+"_X").c_str());
          auto var_y = MakeAsymGraph(sigmaArray_y[1][SampleIterator].get(), sigmaArray_y[2][SampleIterator].get(), sigmaArray_y[3][SampleIterator].get(), (title+"_Y").c_str());

          auto var_x_norm = MakeAsymGraph(sigmaArray_x_norm[1][SampleIterator].get(), sigmaArray_x_norm[2][SampleIterator].get(), sigmaArray_x_norm[3][SampleIterator].get(), (title+"_X_norm").c_str());
          auto var_y_norm = MakeAsymGraph(sigmaArray_y_norm[1][SampleIterator].get(), sigmaArray_y_norm[2][SampleIterator].get(), sigmaArray_y_norm[3][SampleIterator].get(), (title+"_Y_norm").c_str());

          dirArrySample[SampleIterator]->cd();
          var_x->Write();
          var_y->Write();
          var_x_norm->Write();
          var_y_norm->Write();

          //KS: here we loop over all reaction modes defined in "RelevantModes[nRelevantModes]"
          if (DoByMode)
          {
            //KS: this is only relevant if PlotByMode is turned on
            //Checking each mode is time consuming so we only consider one which are relevant for particular analysis
            MaCh3Modes_t RelevantModes[nRelevantModes] = {samples[ivs]->GetMaCh3Modes()->GetMode("CCQE"), samples[ivs]->GetMaCh3Modes()->GetMode("2p2h")};

            for(int ir = 0; ir < nRelevantModes;ir++)
            {
              auto var_mode_x = MakeAsymGraph(sigmaArray_mode_x[1][SampleIterator][ir], sigmaArray_mode_x[2][SampleIterator][ir], sigmaArray_mode_x[3][SampleIterator][ir], (title+"_"+samples[ivs]->GetMaCh3Modes()->GetMaCh3ModeName(RelevantModes[ir])+"_X").c_str());
              auto var_mode_y = MakeAsymGraph(sigmaArray_mode_y[1][SampleIterator][ir], sigmaArray_mode_y[2][SampleIterator][ir], sigmaArray_mode_y[3][SampleIterator][ir], (title+"_"+samples[ivs]->GetMaCh3Modes()->GetMaCh3ModeName(RelevantModes[ir])+"_Y").c_str());

              dirArrySample[SampleIterator]->cd();
              var_mode_x->Write();
              var_mode_y->Write();
            } // end for nRelevantModes
          } // end if mode

          SampleIterator++;
        }//End loop over samples(k)
      }//end looping over sample object
      SampleIterator = 0;
      for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
      {
        for (int k = 0; k < samples[ivs]->GetNsamples(); ++k)
        {
          dirArrySample[SampleIterator]->Close();
          delete dirArrySample[SampleIterator];
          SampleIterator++;
        }
      }

      dirArryDial->Close();
      delete dirArryDial;

      if (DoByMode)
      {
        for (int j = 0; j < numVar; ++j)
        {
          SampleIterator = 0;
          for(unsigned int ivs = 0; ivs < samples.size(); ivs++ )
          {
            for (int k = 0; k < samples[ivs]->GetNsamples(); ++k)
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


// *************************
// For comparison with P-Theta we usually have to apply different parameter values then usual 1, 3 sigma
void FitterBase::CustomRange(const std::string& ParName, const double sigma, double& ParamShiftValue) {
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
  hist->GetYaxis()->SetTitle("Events");
  hist->SetDirectory(nullptr);
  hist->Write(baseName.c_str());
}

// *************************
/// Generic histogram writer - should make main code more palatable
void WriteHistogramsByMode(SampleHandlerFD *sample, const std::string& suffix, const bool by_mode, const bool by_channel) {
// *************************
  std::string sampleName = sample->GetTitle();
  MaCh3Modes *modes = sample->GetMaCh3Modes();

  // Probably a better way of handling this logic
  if (by_mode) {
    for (int iMode = 0; iMode < modes->GetNModes(); ++iMode) {
      auto modeHist = sample->Get1DVarHistByModeAndChannel(sample->GetXBinVarName(), iMode);
      WriteHistograms(modeHist, sampleName + "_" + modes->GetMaCh3ModeName(iMode) + suffix);
      delete modeHist;
    }
  }

  if (by_channel) {
    for (int iChan = 0; iChan < sample->GetNOscChannels(); ++iChan) {
      auto chanHist = sample->Get1DVarHistByModeAndChannel(sample->GetXBinVarName(), -1, iChan); // -1 skips over mode plotting
      WriteHistograms(chanHist, sampleName + "_" + sample->GetFlavourName(iChan) + suffix);
      delete chanHist;
    }
  }

  if (by_mode && by_channel) {
    for (int iMode = 0; iMode < modes->GetNModes(); ++iMode) {
      for (int iChan = 0; iChan < sample->GetNOscChannels(); ++iChan) {
        auto hist = sample->Get1DVarHistByModeAndChannel(sample->GetXBinVarName(), iMode, iChan);
        WriteHistograms(hist, sampleName + "_" + modes->GetMaCh3ModeName(iMode) + "_" + sample->GetFlavourName(iChan) + suffix);
        delete hist;
      }
    }
  }

  if (!by_mode && !by_channel) {
    auto hist = sample->Get1DVarHistByModeAndChannel(sample->GetXBinVarName());
    WriteHistograms(hist, sampleName + suffix);
    delete hist;
  }
}

// *************************
void FitterBase::RunSigmaVarFD() {
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
      double ParamNomValue = systematics[s]->GetParProp(i);
      double ParamLower = systematics[s]->GetLowerBound(i);
      double ParamUpper = systematics[s]->GetUpperBound(i);

      for(unsigned int iSample = 0; iSample < samples.size(); ++iSample)
      {
        if(samples[iSample]->GetNsamples() > 1){
          MACH3LOG_ERROR(":: Sample has more than one sample {} ::", samples[iSample]->GetNsamples());
          throw MaCh3Exception(__FILE__ , __LINE__ );
        }

        auto* MaCh3Sample = dynamic_cast<SampleHandlerFD*>(samples[iSample]);
        if (!MaCh3Sample) {
          MACH3LOG_ERROR(":: Sample {} do not inherit from  SampleHandlerFD this is not implemented::", samples[i]->GetTitle());
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        TDirectory* SampleDir = ParamDir->mkdir(MaCh3Sample->GetTitle().c_str());
        SampleDir->cd();

        for (size_t j = 0; j < SigmaArray.size(); ++j) {
          double sigma = SigmaArray[j];

          double ParamShiftValue = ParamNomValue + sigma * std::sqrt((*systematics[s]->GetCovMatrix())(i,i));
          ParamShiftValue = std::max(std::min(ParamShiftValue, ParamUpper), ParamLower);

          /// Apply custom range to make easier comparison with p-theta
          CustomRange(ParName, sigma, ParamShiftValue);

          MACH3LOG_INFO(
            "  - set to {:<5.2f} ({:<2} sigma shift)",
                        ParamShiftValue,
                        sigma
          );

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

          WriteHistogramsByMode(MaCh3Sample, suffix, plot_by_mode, plot_by_channel);
        }
        SampleDir->Close();
        delete SampleDir;
        ParamDir->cd();
      }

      systematics[s]->SetParProp(i, ParamNomValue);
      MACH3LOG_INFO("  - set back to nominal value {:<5.2f}", ParamNomValue);
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
