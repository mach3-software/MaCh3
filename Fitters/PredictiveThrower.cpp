#include "PredictiveThrower.h"
#include "Samples/SampleHandlerFD.h"
#include "Parameters/ParameterHandlerGeneric.h"
#include "TH3.h"

// *************************
PredictiveThrower::PredictiveThrower(Manager *man) : FitterBase(man) {
// *************************
  AlgorithmName = "PredictiveThrower";
  if(!CheckNodeExists(fitMan->raw(), "Predictive")) {
   MACH3LOG_ERROR("Predictive is missing in your main yaml config");
   throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  StandardFluctuation = GetFromManager<bool>(fitMan->raw()["Predictive"]["StandardFluctuation"], true, __FILE__, __LINE__ );

  if(StandardFluctuation) MACH3LOG_INFO("Using standard method of statistical fluctuation");
  else MACH3LOG_INFO("Using alternative method of statistical fluctuation, which is much slower");

  ModelSystematic = nullptr;
  // Use the full likelihood for the Prior/Posterior predictive pvalue
  FullLLH = GetFromManager<bool>(fitMan->raw()["Predictive"]["FullLLH"], false, __FILE__, __LINE__ );
  NModelParams = 0;

  Is_PriorPredictive = Get<bool>(fitMan->raw()["Predictive"]["PriorPredictive"], __FILE__, __LINE__);
  Ntoys = Get<int>(fitMan->raw()["Predictive"]["Ntoy"], __FILE__, __LINE__);

  ReweightWeight.resize(Ntoys);
  PenaltyTerm.resize(Ntoys);
}

// *************************
// Destructor:
PredictiveThrower::~PredictiveThrower() {
// *************************

}

// *************************
void PredictiveThrower::SetParamters() {
// *************************
  // WARNING This should be removed in the future
  auto DoNotThrowLegacyCov = GetFromManager<std::vector<std::string>>(fitMan->raw()["Predictive"]["DoNotThrowLegacyCov"], {}, __FILE__, __LINE__);
  /// Have ability to not throw legacy matrices
  for (size_t i = 0; i < DoNotThrowLegacyCov.size(); ++i) {
    for (size_t s = 0; s < systematics.size(); ++s) {
      if (systematics[s]->GetName() == DoNotThrowLegacyCov[i]) {
        systematics[s]->SetParameters();
        break;
      }
    }
  }

  // Set groups to prefit values if they were set to not be varies
  if(ModelSystematic && ParameterGroupsNotVaried.size() > 0) {
    ModelSystematic->SetGroupOnlyParameters(ParameterGroupsNotVaried);
  }

  /// Alternatively vary only selected params
  if (ModelSystematic && !ParameterOnlyToVary.empty()) {
    for (int i = 0; i < ModelSystematic->GetNumParams(); ++i) {
      // KS: If parameter is in map then we are skipping this, otherwise for params that we don't want to vary we simply set it to prior
      if (ParameterOnlyToVary.find(i) == ParameterOnlyToVary.end()) {
        ModelSystematic->SetParProp(i, ModelSystematic->GetParInit(i));
      }
    }
  }
}

// *************************
void PredictiveThrower::SetupSampleInformation() {
// *************************
  TotalNumberOfSamples = 0;
  for (size_t iPDF = 0; iPDF < samples.size(); iPDF++) {
    TotalNumberOfSamples += samples[iPDF]->GetNsamples();
  }

  MC_Hist_Toy.resize(TotalNumberOfSamples);
  W2_Hist_Toy.resize(TotalNumberOfSamples);
  Data_Hist.resize(TotalNumberOfSamples);
  MC_Nom_Hist.resize(TotalNumberOfSamples);
  W2_Nom_Hist.resize(TotalNumberOfSamples);

  SampleInfo.resize(TotalNumberOfSamples+1);


  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    MC_Hist_Toy[sample].resize(Ntoys);
    W2_Hist_Toy[sample].resize(Ntoys);
  }
  int counter = 0;
  for (size_t iPDF = 0; iPDF < samples.size(); iPDF++) {
    for (int SampleIndex = 0; SampleIndex < samples[iPDF]->GetNsamples(); ++SampleIndex) {
      SampleInfo[counter].Name = samples[iPDF]->GetSampleTitle(SampleIndex);
      SampleInfo[counter].LocalId = SampleIndex;
      SampleInfo[counter].SamHandler = samples[iPDF];
      SampleInfo[counter].Dimenstion = SampleInfo[counter].SamHandler->GetNDim(SampleIndex);
      counter++;
    }
  }
  SampleInfo[TotalNumberOfSamples].Name= "Total";
}

// *************************
// Produce MaCh3 toys:
void PredictiveThrower::SetupToyGeneration() {
// *************************
  int counter = 0;
  for (size_t s = 0; s < systematics.size(); ++s) {
    auto* MaCh3Params = dynamic_cast<ParameterHandlerGeneric*>(systematics[s]);
    if(MaCh3Params) {
      ModelSystematic = MaCh3Params;
      counter++;
    }
  }

  SetupSampleInformation();

  if(Is_PriorPredictive) {
    MACH3LOG_INFO("You've chosen to run Prior Predictive Distribution");
  } else {
    auto PosteriorFileName = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);
    //KS: We use MCMCProcessor to get names of covariances that were actually used to produce given chain
    MCMCProcessor Processor(PosteriorFileName);
    Processor.Initialise();

    // For throwing FD predictions from ND-only chain we have to allow having different yaml configs
    auto AllowDifferentConfigs = GetFromManager<bool>(fitMan->raw()["Predictive"]["AllowDifferentConfigs"], false, __FILE__, __LINE__);

    ///Let's ask the manager what are the file with covariance matrix
    YAML::Node ConfigInChain = Processor.GetCovConfig(kXSecPar);
    if(ModelSystematic) {
      YAML::Node ConfigNow = ModelSystematic->GetConfig();
      if (!compareYAMLNodes(ConfigNow, ConfigInChain))
      {
        if(AllowDifferentConfigs){
          MACH3LOG_WARN("Yaml configs used for your ParameterHandler for chain you want sample from ({}) and one currently initialised are different", PosteriorFileName);
        } else {
          MACH3LOG_ERROR("Yaml configs used for your ParameterHandler for chain you want sample from ({}) and one currently initialised are different", PosteriorFileName);
          throw MaCh3Exception(__FILE__ , __LINE__ );
        }
      }
    }
  }
  if(counter > 1) {
    MACH3LOG_ERROR("Found {} ParmaterHandler inheriting from ParameterHandlerGeneric, I can accept at most 1", counter);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (size_t s = 0; s < systematics.size(); ++s) {
    NModelParams += systematics[s]->GetNumParams();
  }

  if (ModelSystematic) {
    auto ThrowParamGroupOnly = GetFromManager<std::vector<std::string>>(fitMan->raw()["Predictive"]["ThrowParamGroupOnly"], {}, __FILE__, __LINE__);
    auto UniqueParamGroup = ModelSystematic->GetUniqueParameterGroups();
    auto ParameterOnlyToVaryString = GetFromManager<std::vector<std::string>>(fitMan->raw()["Predictive"]["ThrowSinlgeParams"], {}, __FILE__, __LINE__);

    if (!ThrowParamGroupOnly.empty() && !ParameterOnlyToVaryString.empty()) {
      MACH3LOG_ERROR("Can't use ThrowParamGroupOnly and ThrowSinlgeParams at the same time");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    if (!ParameterOnlyToVaryString.empty()) {
      MACH3LOG_INFO("I will throw only: {}", fmt::join(ParameterOnlyToVaryString, ", "));
      std::vector<int> ParameterVary(ParameterOnlyToVaryString.size());

      for (size_t i = 0; i < ParameterOnlyToVaryString.size(); ++i) {
        ParameterVary[i] = ModelSystematic->GetParIndex(ParameterOnlyToVaryString[i]);
        if (ParameterVary[i] == M3::_BAD_INT_) {
          MACH3LOG_ERROR("Can't proceed if param {} is missing", ParameterOnlyToVaryString[i]);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
      ParameterOnlyToVary = std::unordered_set<int>(ParameterVary.begin(), ParameterVary.end());
    } else {
      MACH3LOG_INFO("I have following parameter groups: {}", fmt::join(UniqueParamGroup, ", "));
      if (ThrowParamGroupOnly.empty()) {
        MACH3LOG_INFO("I will vary all");
      } else {
        std::unordered_set<std::string> throwOnlySet(ThrowParamGroupOnly.begin(), ThrowParamGroupOnly.end());
        ParameterGroupsNotVaried.clear();

        for (const auto& group : UniqueParamGroup) {
          if (throwOnlySet.find(group) == throwOnlySet.end()) {
            ParameterGroupsNotVaried.push_back(group);
          }
        }

        MACH3LOG_INFO("I will vary: {}", fmt::join(ThrowParamGroupOnly, ", "));
        MACH3LOG_INFO("Exclude: {}", fmt::join(ParameterGroupsNotVaried, ", "));
      }
    }
  }
}

// *************************
// Try loading toys
bool PredictiveThrower::LoadToys() {
// *************************
  auto PosteriorFileName = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);
  // Open the ROOT file
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  TFile* file = TFile::Open(PosteriorFileName.c_str(), "READ");

  gErrorIgnoreLevel = originalErrorWarning;
  TDirectory* ToyDir = nullptr;
  if (!file || file->IsZombie()) {
    return false;
  } else {
    // Check for the "toys" directory
    if ((ToyDir = file->GetDirectory("Toys"))) {
      MACH3LOG_INFO("Found toys in Posterior file will attempt toy reading");
    } else {
      file->Close();
      delete file;
      return false;
    }
  }

  // Finally get the TTree branch with the penalty vectors for each of the toy throws
  TTree* PenaltyTree = static_cast<TTree*>(file->Get("ToySummary"));
  if (!PenaltyTree) {
    MACH3LOG_WARN("ToySummary TTree not found in file.");
    file->Close();
    delete file;
    return false;
  }

  Ntoys = static_cast<int>(PenaltyTree->GetEntries());
  int ConfigNtoys = Get<int>(fitMan->raw()["Predictive"]["Ntoy"], __FILE__, __LINE__);;
  if (Ntoys != ConfigNtoys) {
    MACH3LOG_WARN("Found different number of toys in saved file than asked to run!");
    MACH3LOG_INFO("Will read _ALL_ toys in the file");
    MACH3LOG_INFO("Ntoys in file: {}", Ntoys);
    MACH3LOG_INFO("Ntoys specified: {}", ConfigNtoys);
  }

  PenaltyTerm.resize(Ntoys);
  ReweightWeight.resize(Ntoys);

  double Penalty = 0, Weight = 1;
  PenaltyTree->SetBranchAddress("Penalty", &Penalty);
  PenaltyTree->SetBranchAddress("Weight", &Weight);
  PenaltyTree->SetBranchAddress("NModelParams", &NModelParams);

  for (int i = 0; i < Ntoys; ++i) {
    PenaltyTree->GetEntry(i);
    if (FullLLH) {
      PenaltyTerm[i] = Penalty;
    } else {
      PenaltyTerm[i] = 0.0;
    }

    ReweightWeight[i] = Weight;
  }
  // Resize all vectors and get sample names
  SetupSampleInformation();

  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    TH1* DataHist1D = static_cast<TH1*>(ToyDir->Get((SampleInfo[sample].Name + "_data").c_str()));
    Data_Hist[sample] = M3::Clone(DataHist1D);

    TH1* MCHist1D = static_cast<TH1*>(ToyDir->Get((SampleInfo[sample].Name + "_mc").c_str()));
    MC_Nom_Hist[sample] = M3::Clone(MCHist1D);

    TH1* W2Hist1D = static_cast<TH1*>(ToyDir->Get((SampleInfo[sample].Name + "_w2").c_str()));
    W2_Nom_Hist[sample] = M3::Clone(W2Hist1D);
  }


  for (int iToy = 0; iToy < Ntoys; ++iToy)
  {
    if (iToy % 100 == 0) MACH3LOG_INFO("   Loaded toy {}", iToy);

    for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
      TH1* MCHist1D = static_cast<TH1*>(ToyDir->Get((SampleInfo[sample].Name + "_mc_" + std::to_string(iToy)).c_str()));
      TH1* W2Hist1D = static_cast<TH1*>(ToyDir->Get((SampleInfo[sample].Name + "_w2_" + std::to_string(iToy)).c_str()));

      MC_Hist_Toy[sample][iToy] = M3::Clone(MCHist1D);
      W2_Hist_Toy[sample][iToy] = M3::Clone(W2Hist1D);
    }
  }

  file->Close();
  delete file;
  return true;
}

// *************************
std::vector<std::string> PredictiveThrower::GetStoredFancyName(ParameterHandlerBase* Systematics) const {
// *************************
  TDirectory * ogdir = gDirectory;

  std::vector<std::string> FancyNames;
  std::string Name = std::string("Config_") + Systematics->GetName();
  auto PosteriorFileName = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);

  TFile* file = TFile::Open(PosteriorFileName.c_str(), "READ");
  TDirectory* CovarianceFolder = file->GetDirectory("CovarianceFolder");

  TMacro* FoundMacro = static_cast<TMacro*>(CovarianceFolder->Get(Name.c_str()));
  if(FoundMacro == nullptr) {
    file->Close();
    delete file;
    if(ogdir){ ogdir->cd(); }

    return FancyNames;
  }
  MACH3LOG_DEBUG("Found config for {}", Name);
  YAML::Node Settings = TMacroToYAML(*FoundMacro);

  int params = int(Settings["Systematics"].size());
  FancyNames.resize(params);
  int iPar = 0;
  for (auto const &param : Settings["Systematics"]) {
    FancyNames[iPar] = Get<std::string>(param["Systematic"]["Names"]["FancyName"], __FILE__ , __LINE__);
    iPar++;
  }
  file->Close();
  delete file;
  if(ogdir){ ogdir->cd(); }
  return FancyNames;
}


// *************************
void PredictiveThrower::WriteToy(TDirectory* ToyDirectory,
                                 TDirectory* Toy_1DDirectory,
                                 TDirectory* Toy_2DDirectory,
                                 const int iToy) {
// *************************
  int SampleCounter = 0;
  for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
  {
    auto* SampleHandler = samples[iPDF];
    for (int iSample = 0; iSample < SampleHandler->GetNsamples(); ++iSample)
    {
      ToyDirectory->cd();

      auto SampleName = SampleHandler->GetSampleTitle(iSample);
      TH1* MCHist = SampleHandler->GetMCHist(iSample);
      MC_Hist_Toy[SampleCounter][iToy] = M3::Clone(MCHist, SampleName + "_mc_" + std::to_string(iToy));
      MC_Hist_Toy[SampleCounter][iToy]->Write();

      TH1* W2Hist = SampleHandler->GetW2Hist(iSample);
      W2_Hist_Toy[SampleCounter][iToy] = M3::Clone(W2Hist, SampleName + "_w2_" + std::to_string(iToy));
      W2_Hist_Toy[SampleCounter][iToy]->Write();

      // now get 1D projection for every dimension
      Toy_1DDirectory->cd();
      for(int iDim = 0; iDim < SampleHandler->GetNDim(iSample); iDim++) {
        std::string ProjectionName = SampleHandler->GetKinVarName(iSample, iDim);
        std::string ProjectionSuffix = "_1DProj" + std::to_string(iDim) + "_" + std::to_string(iToy);

        auto hist = SampleHandler->Get1DVarHist(iSample, ProjectionName);
        hist->SetTitle((SampleName + ProjectionSuffix).c_str());
        hist->SetName((SampleName + ProjectionSuffix).c_str());
        hist->Write();
        delete hist;
      }

      Toy_2DDirectory->cd();
      // now get 2D projection for every combination
      for(int iDim1 = 0; iDim1 < SampleHandler->GetNDim(iSample); iDim1++) {
        for (int iDim2 = iDim1 + 1; iDim2 < SampleHandler->GetNDim(iSample); ++iDim2) {
          // Get the names for the two dimensions
          std::string XVarName = SampleHandler->GetKinVarName(iSample, iDim1);
          std::string YVarName = SampleHandler->GetKinVarName(iSample, iDim2);

          // Get the 2D histogram for this pair
          auto hist2D = SampleHandler->Get2DVarHist(iSample, XVarName, YVarName);

          // Write the histogram
          std::string suffix2D = "_2DProj_" + std::to_string(iDim1) + "_vs_" + std::to_string(iDim2) + "_" + std::to_string(iToy);
          hist2D->SetTitle((SampleName + suffix2D).c_str());
          hist2D->SetName((SampleName + suffix2D).c_str());
          hist2D->Write();
          delete hist2D;
        }
      }
      SampleCounter++;
    }
  }
}

// *************************
// Produce MaCh3 toys:
void PredictiveThrower::ProduceToys() {
// *************************
  /// If we found toys then skip process of making new toys
  if(LoadToys()) return;

  /// Setup useful information for toy generation
  SetupToyGeneration();

  auto PosteriorFileName = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);

  MACH3LOG_INFO("Starting {}", __func__);

  outputFile->cd();
  double Penalty = 0, Weight = 1;
  int Draw = 0;

  TTree *ToyTree = new TTree("ToySummary", "ToySummary");
  ToyTree->Branch("Penalty", &Penalty, "Penalty/D");
  ToyTree->Branch("Weight", &Weight, "Weight/D");
  ToyTree->Branch("Draw", &Draw, "Draw/I");
  ToyTree->Branch("NModelParams", &NModelParams, "NModelParams/I");

  // KS: define branches so we can keep track of what params we are throwing
  std::vector<double> ParamValues(NModelParams);
  std::vector<const double*> ParampPointers(NModelParams);
  int ParamCounter = 0;
  for (size_t iSys = 0; iSys < systematics.size(); iSys++)
  {
    for (int iPar = 0; iPar < systematics[iSys]->GetNumParams(); iPar++)
    {
      ParampPointers[ParamCounter] = systematics[iSys]->RetPointer(iPar);
      std::string Name = systematics[iSys]->GetParFancyName(iPar);
      //CW: Also strip out - signs because it messes up TBranches
      while (Name.find("-") != std::string::npos) {
        Name.replace(Name.find("-"), 1, std::string("_"));
      }
      ToyTree->Branch(Name.c_str(), &ParamValues[ParamCounter], (Name + "/D").c_str());
      ParamCounter++;
    }
  }
  TDirectory* ToyDirectory = outputFile->mkdir("Toys");
  ToyDirectory->cd();
  int SampleCounter = 0;
  for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
  {
    auto* MaCh3Sample = dynamic_cast<SampleHandlerFD*>(samples[iPDF]);
    for (int SampleIndex = 0; SampleIndex < MaCh3Sample->GetNsamples(); ++SampleIndex)
    {
      // Get nominal spectra and event rates
      TH1* DataHist = MaCh3Sample->GetDataHist(SampleIndex);
      Data_Hist[SampleCounter] = M3::Clone(DataHist, MaCh3Sample->GetSampleTitle(SampleIndex) + "_data");
      Data_Hist[SampleCounter]->Write((MaCh3Sample->GetSampleTitle(SampleIndex) + "_data").c_str());

      TH1* MCHist = MaCh3Sample->GetMCHist(SampleIndex);
      MC_Nom_Hist[SampleCounter] = M3::Clone(MCHist, MaCh3Sample->GetSampleTitle(SampleIndex) + "_mc");
      MC_Nom_Hist[SampleCounter]->Write((MaCh3Sample->GetSampleTitle(SampleIndex) + "_mc").c_str());

      TH1* W2Hist = MaCh3Sample->GetW2Hist(SampleIndex);
      W2_Nom_Hist[SampleCounter] = M3::Clone(W2Hist, MaCh3Sample->GetSampleTitle(SampleIndex) + "_w2");
      W2_Nom_Hist[SampleCounter]->Write((MaCh3Sample->GetSampleTitle(SampleIndex) + "_w2").c_str());
      SampleCounter++;
    }
  }

  TDirectory* Toy_1DDirectory = outputFile->mkdir("Toys_1DHistVar");
  TDirectory* Toy_2DDirectory = outputFile->mkdir("Toys_2DHistVar");
  /// this store value of parameters sampled from a chain
  std::vector<std::vector<double>> branch_vals(systematics.size());
  std::vector<std::vector<std::string>> branch_name(systematics.size());

  TChain* PosteriorFile = nullptr;
  unsigned int burn_in = 0;
  unsigned int maxNsteps = 0;
  unsigned int Step = 0;
  if(!Is_PriorPredictive)
  {
    PosteriorFile = new TChain("posteriors");
    PosteriorFile->Add(PosteriorFileName.c_str());

    PosteriorFile->SetBranchAddress("step", &Step);
    if (PosteriorFile->GetBranch("Weight")) {
      PosteriorFile->SetBranchStatus("Weight", true);
      PosteriorFile->SetBranchAddress("Weight", &Weight);
    } else {
      MACH3LOG_WARN("Not applying reweighting weight");
      Weight = 1.0;
    }

    for (size_t s = 0; s < systematics.size(); ++s) {
      auto fancy_names = GetStoredFancyName(systematics[s]);
      systematics[s]->MatchMaCh3OutputBranches(PosteriorFile, branch_vals[s], branch_name[s], fancy_names);
    }

    //Get the burn-in from the config
    burn_in = Get<unsigned int>(fitMan->raw()["Predictive"]["BurnInSteps"], __FILE__, __LINE__);

    //DL: Adding sanity check for chains shorter than burn in
    maxNsteps = static_cast<unsigned int>(PosteriorFile->GetMaximum("step"));
    if(burn_in >= maxNsteps)
    {
      MACH3LOG_ERROR("You are running on a chain shorter than burn in cut");
      MACH3LOG_ERROR("Maximal value of nSteps: {}, burn in cut {}", maxNsteps, burn_in);
      MACH3LOG_ERROR("You will run into infinite loop");
      MACH3LOG_ERROR("You can make new chain or modify burn in cut");
      throw MaCh3Exception(__FILE__,__LINE__);
    }
  }

  TStopwatch TempClock;
  TempClock.Start();
  for(int i = 0; i < Ntoys; i++)
  {
    if(Ntoys >= 10 && i % (Ntoys/10) == 0) {
      MaCh3Utils::PrintProgressBar(i, Ntoys);
    }

    if(!Is_PriorPredictive){
      int entry = 0;
      Step = 0;

      //YSP: Ensures you get an entry from the mcmc even when burn_in is set to zero (Although not advised :p ).
      //Take 200k burn in steps, WP: Eb C in 1st peaky
      // If we have combined chains by hadd need to check the step in the chain
      // Note, entry is not necessarily same as step due to merged ROOT files, so can't choose entry in the range BurnIn - nEntries :(
      while(Step < burn_in){
        entry = random->Integer(static_cast<unsigned int>(PosteriorFile->GetEntries()));
        PosteriorFile->GetEntry(entry);
      }
      Draw = entry;
    }
    for (size_t s = 0; s < systematics.size(); ++s)
    {
      //KS: Below line can help you get prior predictive distributions which are helpful for getting pre and post ND fit spectra
      //YSP: If not set in the config, the code runs SK Posterior Predictive distributions by default. If true, then the code runs SK prior predictive.
      if(Is_PriorPredictive) {
        systematics[s]->ThrowParameters();
      } else {
        systematics[s]->SetParameters(branch_vals[s]);
      }
    }

    // This set some params to prior value this way you can evaluate errors from subset of errors
    SetParamters();

    Penalty = 0;
    if(FullLLH) {
      for (size_t s = 0; s < systematics.size(); ++s) {
        //KS: do times 2 because banff reports chi2
        Penalty = 2.0 * systematics[s]->GetLikelihood();
      }
    }

    PenaltyTerm[i] = Penalty;
    ReweightWeight[i] = Weight;

    for (size_t iPDF = 0; iPDF < samples.size(); iPDF++) {
      samples[iPDF]->Reweight();
    }
    // Save histograms to file
    WriteToy(ToyDirectory, Toy_1DDirectory, Toy_2DDirectory, i);

    // Fill parameter value so we know throw values
    for (size_t iPar = 0; iPar < ParamValues.size(); iPar++) {
      ParamValues[iPar] = *ParampPointers[iPar];
    }

    ToyTree->Fill();
  }//end of toys loop
  TempClock.Stop();

  if(PosteriorFile) delete PosteriorFile;
  ToyDirectory->Close(); delete ToyDirectory;
  Toy_1DDirectory->Close(); delete Toy_1DDirectory;
  Toy_2DDirectory->Close(); delete Toy_2DDirectory;

  outputFile->cd();
  ToyTree->Write(); delete ToyTree;

  MACH3LOG_INFO("{} took {:.2f}s to finish for {} toys", __func__, TempClock.RealTime(), Ntoys);
}

// *************************
void PredictiveThrower::Study1DProjections(const std::vector<TDirectory*>& SampleDirectories) const {
// *************************
  MACH3LOG_INFO("Starting {}", __func__);

  TDirectory * ogdir = gDirectory;
  auto PosteriorFileName = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);
  // Open the ROOT file
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  TFile* file = TFile::Open(PosteriorFileName.c_str(), "READ");

  gErrorIgnoreLevel = originalErrorWarning;
  TDirectory* ToyDir = file->GetDirectory("Toys_1DHistVar");
  // If toys not amiable in posterior file this means they must be in output file
  if(ToyDir == nullptr) {
    ToyDir = outputFile->GetDirectory("Toys_1DHistVar");
  }
  // [sample], [toy], [dim]
  std::vector<std::vector<std::vector<std::unique_ptr<TH1D>>>> ProjectionToys(TotalNumberOfSamples);
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    ProjectionToys[sample].resize(Ntoys);
    const int nDims = SampleInfo[sample].Dimenstion;
    for (int iToy = 0; iToy < Ntoys; ++iToy) {
      ProjectionToys[sample][iToy].resize(nDims);
    }
  }

  for (int iToy = 0; iToy < Ntoys; ++iToy) {
    if (iToy % 100 == 0) MACH3LOG_INFO("   Loaded Projection toys {}", iToy);
    for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
      const int nDims = SampleInfo[sample].Dimenstion;
      for(int iDim = 0; iDim < nDims; iDim ++){
        std::string ProjectionSuffix = "_1DProj" + std::to_string(iDim) + "_" + std::to_string(iToy);
        TH1D* MCHist1D = static_cast<TH1D*>(ToyDir->Get((SampleInfo[sample].Name + ProjectionSuffix).c_str()));
        ProjectionToys[sample][iToy][iDim] = M3::Clone(MCHist1D);
      }
    } // end loop over samples
  } // end loop over toys
  file->Close(); delete file;
  if(ogdir){ ogdir->cd(); }

  ProduceSpectra(ProjectionToys, SampleDirectories, "mc");
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    const int nDims = SampleInfo[sample].Dimenstion;
    // KS: We only care about doing projections for 2D, for 1D we have well 1D, for beyond 2D we have flattened TH1D
    if(nDims == 2){
      auto hist = Data_Hist[sample].get();
      SampleDirectories[sample]->cd();

      std::string nameX = "Data_" + SampleInfo[sample].Name + "_Dim0";
      std::string nameY = "Data_" + SampleInfo[sample].Name + "_Dim1";

      if(std::string(hist->ClassName()) == "TH2Poly") {
        TAxis* xax = ProjectionToys[sample][0][0]->GetXaxis();
        TAxis* yax = ProjectionToys[sample][0][1]->GetXaxis();

        std::vector<double> XBinning(xax->GetNbins()+1);
        std::vector<double> YBinning(yax->GetNbins()+1);

        for(int i=0;i<=xax->GetNbins();++i)
          XBinning[i] = xax->GetBinLowEdge(i+1);

        for(int i=0;i<=yax->GetNbins();++i)
          YBinning[i] = yax->GetBinLowEdge(i+1);

        TH1D* ProjectionX = PolyProjectionX(static_cast<TH2Poly*>(hist), nameX.c_str(), XBinning, false);
        TH1D* ProjectionY = PolyProjectionY(static_cast<TH2Poly*>(hist), nameY.c_str(), YBinning, false);

        ProjectionX->SetDirectory(nullptr);
        ProjectionY->SetDirectory(nullptr);

        ProjectionX->Write(nameX.c_str());
        ProjectionY->Write(nameY.c_str());

        delete ProjectionX;
        delete ProjectionY;
      } else { //TH2D
        TH1D* ProjectionX = static_cast<TH2D*>(hist)->ProjectionX(nameX.c_str());
        TH1D* ProjectionY = static_cast<TH2D*>(hist)->ProjectionY(nameY.c_str());

        ProjectionX->SetDirectory(nullptr);
        ProjectionY->SetDirectory(nullptr);

        ProjectionX->Write(nameX.c_str());
        ProjectionY->Write(nameY.c_str());
        delete ProjectionX;
        delete ProjectionY;
      }
    }
  }
}

// *************************
void PredictiveThrower::ProduceSpectra(const std::vector<std::vector<std::vector<std::unique_ptr<TH1D>>>>& Toys,
                                       const std::vector<TDirectory*>& SampleDirectories,
                                       const std::string suffix) const {
// *************************
  std::vector<std::vector<double>> MaxValue(TotalNumberOfSamples);

  // 1. Create Max value
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    const int nDims = SampleInfo[sample].Dimenstion;
    MaxValue[sample].assign(nDims, 0);
  }

  // 2. Find maximum entries over all toys
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    for (int toy = 0; toy < Ntoys; ++toy) {
      const int nDims = SampleInfo[sample].Dimenstion;
      for (int dim = 0; dim < nDims; dim++) {
        double max_val = Toys[sample][toy][dim]->GetMaximum();
        MaxValue[sample][dim] = std::max(MaxValue[sample][dim], max_val);
      }
    }
  }

  // 3. Make actual spectra histogram (this is because making ROOT histograms is not save)
  // And we now have actual max values
  std::vector<std::vector<std::unique_ptr<TH2D>>> Spectra(TotalNumberOfSamples);
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    const int nDims = SampleInfo[sample].Dimenstion;
    Spectra[sample].resize(nDims);
    for (int dim = 0; dim < nDims; dim++) {
      // Get MC histogram x-axis binning
      TH1D* refHist = Toys[sample][0][dim].get();

      const int n_bins_x = refHist->GetNbinsX();
      std::vector<double> x_bin_edges(n_bins_x + 1);
      for (int b = 0; b < n_bins_x; ++b) {
        x_bin_edges[b] = refHist->GetXaxis()->GetBinLowEdge(b + 1);
      }
      x_bin_edges[n_bins_x] = refHist->GetXaxis()->GetBinUpEdge(n_bins_x);

      constexpr int n_bins_y = 400;
      constexpr double y_min = 0.0;
      const double y_max = MaxValue[sample][dim] * 1.05;

      // Create TH2D with variable binning on x axis
      Spectra[sample][dim] = std::make_unique<TH2D>(
        (SampleInfo[sample].Name + "_" + suffix + "_dim" + std::to_string(dim)).c_str(),   // name
        (SampleInfo[sample].Name + "_" + suffix + "_dim" + std::to_string(dim)).c_str(),   // title
        n_bins_x, x_bin_edges.data(),                   // x axis bins
        n_bins_y, y_min, y_max                          // y axis bins
      );

      Spectra[sample][dim]->GetXaxis()->SetTitle(refHist->GetXaxis()->GetTitle());
      Spectra[sample][dim]->GetYaxis()->SetTitle("Events");

      Spectra[sample][dim]->SetDirectory(nullptr);
      Spectra[sample][dim]->Sumw2(true);
    }
  }

  // 4. now we can actually fill our projections
  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2)
  #endif
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    for (int toy = 0; toy < Ntoys; ++toy) {
      const int nDims = SampleInfo[sample].Dimenstion;
      for (int dim = 0; dim < nDims; dim++) {
        FastViolinFill(Spectra[sample][dim].get(), Toys[sample][toy][dim].get());
      }
    }
  }

  // 5. Save histograms which is not thread save
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    SampleDirectories[sample]->cd();
    const int nDims = SampleInfo[sample].Dimenstion;
    for (long unsigned int dim = 0; dim < Spectra[sample].size(); dim++) {
      Spectra[sample][dim]->Write();
      // For case of 2D make additional histograms
      if(nDims == 2) {
        const std::string name = SampleInfo[sample].Name + "_" + suffix+ "_PostPred_dim" + std::to_string(dim);
        auto Summary = MakeSummaryFromSpectra(Spectra[sample][dim].get(), name);
        Summary->Write();
      }
    }
  }
}

// *************************
std::string PredictiveThrower::GetBinName(TH1* hist,
                                          const bool uniform,
                                          const int Dim,
                                          const std::vector<int>& bins) const {
// *************************
  std::string BinName = "";
  if(Dim == 1) { // True 1D distribution using TH1D
    const int b = bins[0];
    const TAxis* ax = hist->GetXaxis();
    const double low = ax->GetBinLowEdge(b);
    const double up  = ax->GetBinUpEdge(b);

    BinName = fmt::format("Dim0 ({:g}, {:g})", low, up);
  } else if (Dim == 2) { // True 2D dsitrubitons
    if(uniform == true) { //using TH2D
      const int bx = bins[0];
      const int by = bins[1];
      const TAxis* ax = hist->GetXaxis();
      const TAxis* ay = hist->GetYaxis();
      BinName = fmt::format("Dim0 ({:g}, {:g}), ", ax->GetBinLowEdge(bx), ax->GetBinUpEdge(bx));
      BinName += fmt::format("Dim1 ({:g}, {:g})", ay->GetBinLowEdge(by), ay->GetBinUpEdge(by));
    } else { // using TH2Poly
      TH2PolyBin* bin = static_cast<TH2PolyBin*>(static_cast<TH2Poly*>(hist)->GetBins()->At(bins[0]-1));
      // Just make a little fancy name
      BinName += fmt::format("Dim{} ({:g}, {:g})", 0, bin->GetXMin(), bin->GetXMax());
      BinName += fmt::format("Dim{} ({:g}, {:g})", 1, bin->GetYMin(), bin->GetYMax());
    }
  } else { // N-dimensional distribution using flatten TH1D
    BinName = hist->GetXaxis()->GetBinLabel(bins[0]);
  }
  return BinName;
}

// *************************
std::vector<std::unique_ptr<TH1D>> PredictiveThrower::PerBinHistogram(TH1* hist,
                                                                      const int SampleId,
                                                                      const int Dim,
                                                                      const std::string& suffix) const {
// *************************
  std::vector<std::unique_ptr<TH1D>> PosteriorHistVec;
  constexpr int nBins = 100;
  const std::string Sample_Name = SampleInfo[SampleId].Name;
  if (Dim == 2) {
      if(std::string(hist->ClassName()) == "TH2Poly") {
        for (int i = 1; i <= static_cast<TH2Poly*>(hist)->GetNumberOfBins(); ++i) {
          std::string ProjName = fmt::format("{} {} Bin: {}",
                                             Sample_Name, suffix,
                                             GetBinName(hist, false, Dim, {i}));
          //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
          // https://root.cern.ch/doc/master/classTH1.html#auto-bin
          auto PosteriorHist = std::make_unique<TH1D>(ProjName.c_str(), ProjName.c_str(), nBins, 1, -1);
          PosteriorHist->SetDirectory(nullptr);
          PosteriorHist->GetXaxis()->SetTitle("Events");
          PosteriorHistVec.push_back(std::move(PosteriorHist));
        } //end loop over bin
      } else {
        int nbinsx = hist->GetNbinsX();
        int nbinsy = hist->GetNbinsY();
        for (int iy = 1; iy <= nbinsy; ++iy) {
          for (int ix = 1; ix <= nbinsx; ++ix) {
            std::string ProjName = fmt::format("{} {} Bin: {}",
                                              Sample_Name, suffix,
                                              GetBinName(hist, true, Dim, {ix,iy}));
            //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
            // https://root.cern.ch/doc/master/classTH1.html#auto-bin
            auto PosteriorHist = std::make_unique<TH1D>(ProjName.c_str(), ProjName.c_str(), nBins, 1, -1);
            PosteriorHist->SetDirectory(nullptr);
            PosteriorHist->GetXaxis()->SetTitle("Events");
            PosteriorHistVec.push_back(std::move(PosteriorHist));
          }
        }
      }
  } else {
    int nbinsx = hist->GetNbinsX();
    PosteriorHistVec.reserve(nbinsx);
    for (int i = 1; i <= nbinsx; ++i) {
      std::string ProjName = fmt::format("{} {} Bin: {}",
                                          Sample_Name, suffix,
                                          GetBinName(hist, true, Dim, {i}));
      //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
      // https://root.cern.ch/doc/master/classTH1.html#auto-bin
      auto PosteriorHist = std::make_unique<TH1D>(ProjName.c_str(), ProjName.c_str(), nBins, 1, -1);
      PosteriorHist->SetDirectory(nullptr);
      PosteriorHist->GetXaxis()->SetTitle("Events");
      PosteriorHistVec.push_back(std::move(PosteriorHist));
    }
  }
  return PosteriorHistVec;
}

// *************************
std::vector<std::unique_ptr<TH1>> PredictiveThrower::MakePredictive(const std::vector<std::vector<std::unique_ptr<TH1>>>& Toys,
                                                                    const std::vector<TDirectory*>& Directory,
                                                                    const std::string& suffix,
                                                                    const bool DebugHistograms) {
// *************************
  std::vector<std::unique_ptr<TH1>> PostPred(TotalNumberOfSamples);
  std::vector<std::vector<std::unique_ptr<TH1D>>> Posterior_hist(TotalNumberOfSamples);
  // 1.initialisation
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    const int nDims = SampleInfo[sample].Dimenstion;
    const std::string Sample_Name = SampleInfo[sample].Name;
    Posterior_hist[sample] = PerBinHistogram(Toys[sample][0].get(), sample, nDims, suffix);
    auto PredictiveHist = M3::Clone(Toys[sample][0].get());
    // Clear the bin contents
    PredictiveHist->Reset();
    PredictiveHist->SetName((Sample_Name + "_" + suffix + "_PostPred").c_str());
    PredictiveHist->SetTitle((Sample_Name + "_" + suffix + "_PostPred").c_str());
    PredictiveHist->SetDirectory(nullptr);
    PostPred[sample] = std::move(PredictiveHist);
  }

  /// 2. Fill histograms, thread safe as all histograms are allocated before and we loop over samples
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    const int nDims = SampleInfo[sample].Dimenstion;
    auto& hist = Toys[sample][0];
    for (size_t iToy = 0; iToy < Toys[sample].size(); ++iToy) {
      if(nDims == 2) {
        if(std::string(hist->ClassName()) == "TH2Poly") {
          for (int i = 1; i <= static_cast<TH2Poly*>(hist.get())->GetNumberOfBins(); ++i) {
            double content = Toys[sample][iToy]->GetBinContent(i);
            Posterior_hist[sample][i-1]->Fill(content, ReweightWeight[iToy]);
          }
        } else {
          int nbinsx = hist->GetNbinsX();
          int nbinsy = hist->GetNbinsY();
          for (int iy = 1; iy <= nbinsy; ++iy) {
            for (int ix = 1; ix <= nbinsx; ++ix) {
              int Bin = (iy-1) * nbinsx + (ix-1);
              double content = Toys[sample][iToy]->GetBinContent(ix, iy);
              Posterior_hist[sample][Bin];
              Posterior_hist[sample][Bin]->Fill(content, ReweightWeight[iToy]);
            } // end loop over X bins
          } // end loop over Y bins
        }
      } else {
        int nbinsx = hist->GetNbinsX();
          for (int i = 1; i <= nbinsx; ++i) {
            double content = Toys[sample][iToy]->GetBinContent(i);
            Posterior_hist[sample][i-1]->Fill(content, ReweightWeight[iToy]);
          } // end loop over bins
        } // end if over dimensions
    } // end loop over toys
  } // end loop over samples

  // 3.save
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    const int nDims = SampleInfo[sample].Dimenstion;
    auto& hist = Toys[sample][0];
    Directory[sample]->cd();
    if(nDims == 2) {
      if(std::string(hist->ClassName()) == "TH2Poly") {
        for (int i = 1; i <= static_cast<TH2Poly*>(hist.get())->GetNumberOfBins(); ++i) {
          PostPred[sample]->SetBinContent(i, Posterior_hist[sample][i-1]->GetMean());
          // KS: If ROOT below 6.18 one need -1 only for error due to stupid bug...
          PostPred[sample]->SetBinError(i, Posterior_hist[sample][i-1]->GetRMS());
          if (DebugHistograms) Posterior_hist[sample][i-1]->Write();
        } // end loop over poly bins
      } else {
        int nbinsx = hist->GetNbinsX();
        int nbinsy = hist->GetNbinsY();
        for (int iy = 1; iy <= nbinsy; ++iy) {
          for (int ix = 1; ix <= nbinsx; ++ix) {
            int Bin = (iy-1) * nbinsx + (ix-1);
            if (DebugHistograms) Posterior_hist[sample][Bin]->Write();
            PostPred[sample]->SetBinContent(ix, iy, Posterior_hist[sample][Bin]->GetMean());
            PostPred[sample]->SetBinError(ix, iy, Posterior_hist[sample][Bin]->GetRMS());
          } // end loop over x
        } // end loop over y
      }
    } else {
      int nbinsx = hist->GetNbinsX();
      for (int i = 1; i <= nbinsx; ++i) {
        PostPred[sample]->SetBinContent(i, Posterior_hist[sample][i-1]->GetMean());
        PostPred[sample]->SetBinError(i, Posterior_hist[sample][i-1]->GetRMS());
        if (DebugHistograms) Posterior_hist[sample][i-1]->Write();
      }
    }
    PostPred[sample]->Write();
  } // end loop over samples
  return PostPred;
}

// *************************
// Perform predictive analysis
void PredictiveThrower::RunPredictiveAnalysis() {
// *************************
  MACH3LOG_INFO("Starting {}", __func__);
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  TStopwatch TempClock;
  TempClock.Start();

  auto DebugHistograms = GetFromManager<bool>(fitMan->raw()["Predictive"]["DebugHistograms"], false, __FILE__, __LINE__);
  auto StudyBeta = GetFromManager<bool>(fitMan->raw()["Predictive"]["StudyBetaParameters"], true, __FILE__, __LINE__);

  TDirectory* PredictiveDir = outputFile->mkdir("Predictive");
  std::vector<TDirectory*> SampleDirectories;
  SampleDirectories.resize(TotalNumberOfSamples+1);

  // open directory for every sample
  for (int sample = 0; sample < TotalNumberOfSamples+1; ++sample) {
    SampleDirectories[sample] = PredictiveDir->mkdir(SampleInfo[sample].Name.c_str());
  }

  // Produce Violin style spectra
  Study1DProjections(SampleDirectories);
  // Produce posterior predictive distribution
  auto PostPred_mc = MakePredictive(MC_Hist_Toy, SampleDirectories, "mc", DebugHistograms);
  // Calculate Posterior Predictive $p$-value
  PosteriorPredictivepValue(PostPred_mc, SampleDirectories);

  // Close directories
  for (int sample = 0; sample < TotalNumberOfSamples+1; ++sample) {
    SampleDirectories[sample]->Close();
    delete SampleDirectories[sample];
  }
  // Perform beta analysis for mc statical uncertainty
  if(StudyBeta) StudyBetaParameters(PredictiveDir);

  PredictiveDir->Close();
  delete PredictiveDir;

  outputFile->cd();

  TempClock.Stop();
  MACH3LOG_INFO("{} took {:.2f}s to finish for {} toys", __func__, TempClock.RealTime(), Ntoys);
}

// *************************
double PredictiveThrower::GetLLH(const std::unique_ptr<TH1>& DatHist,
                                 const std::unique_ptr<TH1>& MCHist,
                                 const std::unique_ptr<TH1>& W2Hist,
                                 const SampleHandlerBase* SampleHandler) {
// *************************
  double llh = 0.0;
  for (int i = 1; i <= DatHist->GetXaxis()->GetNbins(); ++i)
  {
    const double data = DatHist->GetBinContent(i);
    const double mc = MCHist->GetBinContent(i);
    const double w2 = W2Hist->GetBinContent(i);
    llh += SampleHandler->GetTestStatLLH(data, mc, w2);
  }
  //KS: do times 2 because banff reports chi2
  return 2*llh;
}

// ****************
//KS: We have two methods how to apply statistical fluctuation standard is faster hence is default
void PredictiveThrower::MakeFluctuatedHistogram(TH1* FluctHist, TH1* Hist, const int nDims) {
// ****************
  if(StandardFluctuation){
    if (nDims == 2) {
      if(std::string(Hist->ClassName()) == "TH2Poly") {
        MakeFluctuatedHistogramStandard(static_cast<TH2Poly*>(FluctHist), static_cast<TH2Poly*>(Hist), random.get());
      } else {
        MakeFluctuatedHistogramStandard(static_cast<TH2D*>(FluctHist), static_cast<TH2D*>(Hist), random.get());
      }
    } else {
      MakeFluctuatedHistogramStandard(static_cast<TH1D*>(FluctHist), static_cast<TH1D*>(Hist), random.get());
    }
  } else {
    if (nDims == 2) {
      if(std::string(Hist->ClassName()) == "TH2Poly") {
        MakeFluctuatedHistogramAlternative(static_cast<TH2Poly*>(FluctHist), static_cast<TH2Poly*>(Hist), random.get());
      } else {
        MakeFluctuatedHistogramAlternative(static_cast<TH2D*>(FluctHist), static_cast<TH2D*>(Hist), random.get());
      }
    } else {
        MakeFluctuatedHistogramAlternative(static_cast<TH1D*>(FluctHist), static_cast<TH1D*>(Hist), random.get());
    }
  }
}

// *************************
void PredictiveThrower::PosteriorPredictivepValue(const std::vector<std::unique_ptr<TH1>>& PostPred_mc,
                                                  const std::vector<TDirectory*>& SampleDir) {
// *************************
  // [Toys][Sample]
  std::vector<std::vector<double>> chi2_dat_vec(Ntoys);
  std::vector<std::vector<double>> chi2_mc_vec(Ntoys);
  std::vector<std::vector<double>> chi2_pred_vec(Ntoys);

  for(int iToy = 0; iToy < Ntoys; iToy++) {
    chi2_dat_vec[iToy].resize(TotalNumberOfSamples+1, 0);
    chi2_mc_vec[iToy].resize(TotalNumberOfSamples+1, 0);
    chi2_pred_vec[iToy].resize(TotalNumberOfSamples+1, 0);

    chi2_dat_vec[iToy].back() = PenaltyTerm[iToy];
    chi2_mc_vec[iToy].back() = PenaltyTerm[iToy];
    chi2_pred_vec[iToy].back() = PenaltyTerm[iToy];

    /// TODO This can be multithreaded but be careful for Clone!!!
    for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
      const int nDims = SampleInfo[iSample].Dimenstion;

      auto DrawFluctHist = M3::Clone(MC_Hist_Toy[iSample][iToy].get());
      auto PredFluctHist = M3::Clone(PostPred_mc[iSample].get());

      MakeFluctuatedHistogram(DrawFluctHist.get(), MC_Hist_Toy[iSample][iToy].get(), nDims);
      MakeFluctuatedHistogram(PredFluctHist.get(), PostPred_mc[iSample].get(), nDims);

      // Okay now we can do our chi2 calculation for our sample
      chi2_dat_vec[iToy][iSample]  = GetLLH(Data_Hist[iSample], MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], SampleInfo[iSample].SamHandler);
      chi2_mc_vec[iToy][iSample]   = GetLLH(DrawFluctHist, MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], SampleInfo[iSample].SamHandler);
      chi2_pred_vec[iToy][iSample] = GetLLH(PredFluctHist, MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], SampleInfo[iSample].SamHandler);

      chi2_dat_vec[iToy].back()  += chi2_dat_vec[iToy][iSample];
      chi2_mc_vec[iToy].back()   += chi2_mc_vec[iToy][iSample];
      chi2_pred_vec[iToy].back() += chi2_pred_vec[iToy][iSample];
    }
  }

  MakeChi2Plots(chi2_mc_vec,   "-2LLH (Draw Fluc, Draw)", chi2_dat_vec, "-2LLH (Data, Draw)", SampleDir, "_drawfluc_draw");
  MakeChi2Plots(chi2_pred_vec, "-2LLH (Pred Fluc, Draw)", chi2_dat_vec, "-2LLH (Data, Draw)", SampleDir, "_predfluc_draw");
}

// *************************
void PredictiveThrower::MakeChi2Plots(const std::vector<std::vector<double>>& Chi2_x,
                                      const std::string& Chi2_x_title,
                                      const std::vector<std::vector<double>>& Chi2_y,
                                      const std::string& Chi2_y_title,
                                      const std::vector<TDirectory*>& SampleDir,
                                      const std::string Title) {
// *************************
  for (int iSample = 0; iSample < TotalNumberOfSamples+1; ++iSample) {
    SampleDir[iSample]->cd();

    // Transpose to extract chi2 values for a given sample across all toys
    std::vector<double> chi2_y_sample(Ntoys);
    std::vector<double> chi2_x_per_sample(Ntoys);

    for (int iToy = 0; iToy < Ntoys; ++iToy) {
      chi2_y_sample[iToy] = Chi2_y[iToy][iSample];
      chi2_x_per_sample[iToy]  = Chi2_x[iToy][iSample];
    }

    const double min_val = std::min(*std::min_element(chi2_y_sample.begin(), chi2_y_sample.end()),
                              *std::min_element(chi2_x_per_sample.begin(), chi2_x_per_sample.end()));
    const double max_val = std::max(*std::max_element(chi2_y_sample.begin(), chi2_y_sample.end()),
                              *std::max_element(chi2_x_per_sample.begin(), chi2_x_per_sample.end()));

    auto chi2_hist = std::make_unique<TH2D>((SampleInfo[iSample].Name+ Title).c_str(),
                                            (SampleInfo[iSample].Name+ Title).c_str(),
                                            100, min_val, max_val, 100, min_val, max_val);
    chi2_hist->SetDirectory(nullptr);
    chi2_hist->GetXaxis()->SetTitle(Chi2_x_title.c_str());
    chi2_hist->GetYaxis()->SetTitle(Chi2_y_title.c_str());

    for (int iToy = 0; iToy < Ntoys; ++iToy) {
      chi2_hist->Fill(chi2_x_per_sample[iToy], chi2_y_sample[iToy]);
    }

    Get2DBayesianpValue(chi2_hist.get());
    chi2_hist->Write();
  }
}

// *************************
// Study Beta Parameters
void PredictiveThrower::StudyBetaParameters(TDirectory* PredictiveDir) {
// *************************
  bool StudyBeta = GetFromManager<bool>(fitMan->raw()["Predictive"]["StudyBetaParameters"], true, __FILE__, __LINE__ );
  if (StudyBeta == false) return;

  MACH3LOG_INFO("Starting {}", __func__);
  TDirectory* BetaDir = PredictiveDir->mkdir("BetaParameters");
  std::vector<std::vector<std::unique_ptr<TH1D>>> BetaHist(TotalNumberOfSamples);
  std::vector<TDirectory *> DirBeta(TotalNumberOfSamples);
  // initialise directory for each sample
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    BetaDir->cd();
    DirBeta[sample] = BetaDir->mkdir(SampleInfo[sample].Name.c_str());
  }

  /// 1. Initialise Beta histogram
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    const int nDims = SampleInfo[iSample].Dimenstion;
    // Use any histogram that defines the binning structure
    TH1* RefHist = Data_Hist[iSample].get();
    BetaHist[iSample] = PerBinHistogram(RefHist, iSample, nDims, "Beta_Parameter");
    // Change x-axis title
    for (size_t i = 0; i < BetaHist[iSample].size(); ++i) {
      BetaHist[iSample][i]->GetXaxis()->SetTitle("beta parameter");
    }
  }

  /// 2. Fill histograms, thread safe as all histograms are allocated before and we loop over samples
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    const int nDims = SampleInfo[iSample].Dimenstion;
    const auto likelihood = SampleInfo[iSample].SamHandler->GetTestStatistic();
    for (int iToy = 0; iToy < Ntoys; ++iToy) {
      if (nDims == 2) {
        if(std::string(Data_Hist[iSample]->ClassName()) == "TH2Poly") {
          for (int i = 1; i <= static_cast<TH2Poly*>(Data_Hist[iSample].get())->GetNumberOfBins(); ++i) {
            const double Data = Data_Hist[iSample]->GetBinContent(i);
            const double MC   = MC_Hist_Toy[iSample][iToy]->GetBinContent(i);
            const double w2   = W2_Hist_Toy[iSample][iToy]->GetBinContent(i);

            const double BetaParam = GetBetaParameter(Data, MC, w2, likelihood);
            BetaHist[iSample][i-1]->Fill(BetaParam, ReweightWeight[iToy]);
          } // end loop over poly bins
        } else {
          const int nX = Data_Hist[iSample]->GetNbinsX();
          const int nY = Data_Hist[iSample]->GetNbinsY();
          for (int iy = 1; iy <= nY; ++iy) {
            for (int ix = 1; ix <= nX; ++ix) {
              const int FlatBin = (iy-1) * nX + (ix-1);

              const double Data = Data_Hist[iSample]->GetBinContent(ix, iy);
              const double MC   = MC_Hist_Toy[iSample][iToy]->GetBinContent(ix, iy);
              const double w2   = W2_Hist_Toy[iSample][iToy]->GetBinContent(ix, iy);

              const double BetaParam = GetBetaParameter(Data, MC, w2, likelihood);
              BetaHist[iSample][FlatBin]->Fill(BetaParam, ReweightWeight[iToy]);
            }
          } // end loop over x
        } // end loop over y
      } else {
        int nbinsx = Data_Hist[iSample]->GetNbinsX();
        for (int ix = 1; ix <= nbinsx; ++ix) {
          /// ROOT enumerates from 1 while MaCh3 from 0
          const double Data = Data_Hist[iSample]->GetBinContent(ix);
          const double MC = MC_Hist_Toy[iSample][iToy]->GetBinContent(ix);
          const double w2 = W2_Hist_Toy[iSample][iToy]->GetBinContent(ix);

          const double BetaParam = GetBetaParameter(Data, MC, w2, likelihood);
          BetaHist[iSample][ix-1]->Fill(BetaParam, ReweightWeight[iToy]);
        } // end loop over bins
      }
    } // end loop over toys
  } // end loop over samples

  /// 3. Write to file
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    for (size_t iBin = 0; iBin < BetaHist[iSample].size(); iBin++) {
      DirBeta[iSample]->cd();
      BetaHist[iSample][iBin]->Write();
    }
    DirBeta[iSample]->Close();
    delete DirBeta[iSample];
  }
  BetaDir->Close();
  delete BetaDir;

  PredictiveDir->cd();
}
