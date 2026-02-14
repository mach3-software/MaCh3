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
  for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
  {
    auto* MaCh3Sample = dynamic_cast<SampleHandlerFD*>(samples[iPDF]);
    if (!MaCh3Sample) {
      MACH3LOG_ERROR("Sample {} do not inherit from SampleHandlerFD this is not implemented", samples[iPDF]->GetName());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
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
      SampleInfo[counter].SamHandler = dynamic_cast<SampleHandlerFD*>(samples[iPDF]);
      SampleInfo[counter].Binning = SampleInfo[counter].SamHandler->GetBinningHandler();
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

    SampleCounter = 0;
    for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
    {
      auto* MaCh3Sample = dynamic_cast<SampleHandlerFD*>(samples[iPDF]);
      for (int SampleIndex = 0; SampleIndex < MaCh3Sample->GetNsamples(); ++SampleIndex)
      {
        TH1* MCHist = MaCh3Sample->GetMCHist(SampleIndex);
        MC_Hist_Toy[SampleCounter][i] = M3::Clone(MCHist, MaCh3Sample->GetSampleTitle(SampleIndex) + "_mc_" + std::to_string(i));
        MC_Hist_Toy[SampleCounter][i]->Write();

        TH1* W2Hist = MaCh3Sample->GetW2Hist(SampleIndex);
        W2_Hist_Toy[SampleCounter][i] = M3::Clone(W2Hist, MaCh3Sample->GetSampleTitle(SampleIndex) + "_w2_" + std::to_string(i));
        W2_Hist_Toy[SampleCounter][i]->Write();
        SampleCounter++;
      }
    }

    // Fill parameter value so we know throw values
    for (size_t iPar = 0; iPar < ParamValues.size(); iPar++) {
      ParamValues[iPar] = *ParampPointers[iPar];
    }

    ToyTree->Fill();
  }//end of toys loop
  TempClock.Stop();

  if(PosteriorFile) delete PosteriorFile;
  ToyDirectory->Close();
  delete ToyDirectory;

  outputFile->cd();
  ToyTree->Write();
  delete ToyTree;

  MACH3LOG_INFO("{} took {:.2f}s to finish for {} toys", __func__, TempClock.RealTime(), Ntoys);
}

// *************************
std::vector<std::vector<std::unique_ptr<TH2D>>> PredictiveThrower::ProduceSpectra(
                                                    const std::vector<std::vector<std::unique_ptr<TH1>>>& Toys,
                                                    const std::vector<TDirectory*>& SampleDirectories,
                                                    const std::string suffix) {
// *************************
  std::vector<std::vector<double>> MaxValue(TotalNumberOfSamples);
  //Projections[sample][toy][dim]
  std::vector<std::vector<std::vector<std::unique_ptr<TH1>>>> Projections(TotalNumberOfSamples);

  // 1. Create 1D projections, this is not thread safe, also we will use it a lot later on
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    SampleDirectories[sample]->cd();

    const int nDims = SampleInfo[sample].Dimenstion;
    MaxValue[sample].assign(nDims, 0);
    Projections[sample].resize(Ntoys);
    for (int toy = 0; toy < Ntoys; ++toy) {
      if (nDims == 1) {
        // For 1D histograms, no projection needed.
        // Leave Projections[sample][toy][0] == nullptr
      } else if (nDims == 2) {
        Projections[sample][toy].resize(nDims);
        TH2* h2 = dynamic_cast<TH2*>(Toys[sample][toy].get());
        if (!h2) {
          MACH3LOG_ERROR("Histogram is not TH2 for nDims==2");
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        auto px = h2->ProjectionX();
        px->SetDirectory(nullptr);
        Projections[sample][toy][0] = M3::Clone(px);

        auto py = h2->ProjectionY();
        py->SetDirectory(nullptr);
        Projections[sample][toy][1] = M3::Clone(py);

        delete px;
        delete py;
      } else {
        MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, nDims);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }

  // 2. Find maximum entries over all toys
  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2)
  #endif
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    for (int toy = 0; toy < Ntoys; ++toy) {
      const int nDims = Toys[sample][0]->GetDimension();
      for (int dim = 0; dim < nDims; dim++) {
        double max_val = 0;
        if (nDims == 1) {
          max_val = Toys[sample][toy]->GetMaximum();
        }
        else if (nDims == 2) {
          if (dim == 0) {
            max_val = Projections[sample][toy][0]->GetMaximum();
          } else {
            max_val = Projections[sample][toy][1]->GetMaximum();
          }
        } else {
          MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, nDims);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
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
      TH1D* refHist = nullptr;
      if (nDims == 1) {
        refHist = static_cast<TH1D*>(Toys[sample][0].get());
      }
      else if (nDims == 2) {
        if (dim == 0) {
          refHist = static_cast<TH1D*>(Projections[sample][0][0].get());
        } else {
          refHist = static_cast<TH1D*>(Projections[sample][0][1].get());
        }
      }
      else {
        MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, nDims);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      if (!refHist) {
        MACH3LOG_ERROR("Failed to cast to {} dimensions in {}!", nDims, __func__);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      const int n_bins_x = refHist->GetNbinsX();
      std::vector<double> x_bin_edges(n_bins_x + 1);
      for (int b = 0; b <= n_bins_x; ++b) {
        x_bin_edges[b] = refHist->GetXaxis()->GetBinLowEdge(b + 1); // ROOT bins start at 1
      }
      // Last edge is upper edge of last bin:
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
  #pragma omp parallel for
  #endif
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    const int nDims = SampleInfo[sample].Dimenstion;
    for (int dim = 0; dim < nDims; dim++) {
      for (int toy = 0; toy < Ntoys; ++toy) {
        if (nDims == 1) {
          FastViolinFill(Spectra[sample][dim].get(), static_cast<TH1D*>(Toys[sample][toy].get()));
        }
        else if (nDims == 2) {
          FastViolinFill(Spectra[sample][dim].get(), static_cast<TH1D*>(Projections[sample][toy][dim].get()));
        }
        else {
          MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, nDims);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
    }
  }

  // 5. Save histograms which is not thread save
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    SampleDirectories[sample]->cd();
    for (long unsigned int dim = 0; dim < Spectra[sample].size(); dim++) {
      Spectra[sample][dim]->Write();
    }
  }

  return Spectra;
}


// *************************
std::vector<std::unique_ptr<TH1>> PredictiveThrower::MakePredictive(const std::vector<std::vector<std::unique_ptr<TH1>>>& Toys,
                                                                    const std::vector<TDirectory*>& Directory,
                                                                    const std::string& suffix,
                                                                    const bool DebugHistograms) {
// *************************
  /// @todo KS: If we break up into 1.initialisation->2.writing->3.saving into separate loops we can multithread
  std::vector<std::unique_ptr<TH1>> PostPred_mc(TotalNumberOfSamples);

  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    Directory[sample]->cd();
    const int nDims = SampleInfo[sample].Dimenstion;
    const std::string Sample_Name = SampleInfo[sample].Name;
    if (nDims == 1) {
      int nbinsx = Toys[sample][0]->GetNbinsX();
      const double* xbins = Toys[sample][0]->GetXaxis()->GetXbins()->GetArray();
      auto PredictiveHist = std::make_unique<TH1D>(
        (Sample_Name + "_" + suffix + "_PostPred").c_str(),
        (Sample_Name + "_" + suffix + "_PostPred").c_str(),
         nbinsx, xbins
      );
      PredictiveHist->GetXaxis()->SetTitle(Toys[sample][0]->GetXaxis()->GetTitle());
      PredictiveHist->GetYaxis()->SetTitle("Events");
      PredictiveHist->SetDirectory(nullptr);

      for (int i = 1; i <= nbinsx; ++i) {
        std::string ProjName = fmt::format("{} {} Bin: {}",
                                           SampleInfo[sample].Name, suffix,
                                           SampleInfo[sample].Binning->GetBinName(SampleInfo[sample].LocalId, i-1));
        //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
        // https://root.cern.ch/doc/master/classTH1.html#auto-bin
        auto PosteriorHist = std::make_unique<TH1D>(ProjName.c_str(), ProjName.c_str(), 100, 1, -1);
        PosteriorHist->SetDirectory(nullptr);
        PosteriorHist->GetXaxis()->SetTitle("Events");

        for (size_t iToy = 0; iToy < Toys[sample].size(); ++iToy) {
          double content = Toys[sample][iToy]->GetBinContent(i);
          PosteriorHist->Fill(content, ReweightWeight[iToy]);
        }

        if (DebugHistograms) PosteriorHist->Write();

        PredictiveHist->SetBinContent(i, PosteriorHist->GetMean());
        PredictiveHist->SetBinError(i, PosteriorHist->GetRMS());
      }
      PredictiveHist->Write();
      PostPred_mc[sample] = std::move(PredictiveHist);
    }
    else if (nDims == 2) {
      int nbinsx = Toys[sample][0]->GetNbinsX();
      int nbinsy = Toys[sample][0]->GetNbinsY();
      const double* xbins = Toys[sample][0]->GetXaxis()->GetXbins()->GetArray();
      const double* ybins = Toys[sample][0]->GetYaxis()->GetXbins()->GetArray();

      // Create 2D predictive histogram with same binning as Toys[sample][0]
      auto PredictiveHist = std::make_unique<TH2D>(
        (Sample_Name + "_" + suffix + "_PostPred").c_str(),
        (Sample_Name + "_" + suffix + "_PostPred").c_str(),
        nbinsx, xbins,
        nbinsy, ybins
      );
      PredictiveHist->GetXaxis()->SetTitle(Toys[sample][0]->GetXaxis()->GetTitle());
      PredictiveHist->GetYaxis()->SetTitle(Toys[sample][0]->GetYaxis()->GetTitle());
      PredictiveHist->SetDirectory(nullptr);

      for (int iy = 1; iy <= nbinsy; ++iy) {
        for (int ix = 1; ix <= nbinsx; ++ix) {
          // MaCh3 vs ROOT conventions, above loop should be over MaCh3 bins
          const int MaCh3Bin = SampleInfo[sample].Binning->GetBinSafe(SampleInfo[sample].LocalId, {ix-1, iy-1});
          std::string ProjName = fmt::format("{} {} Bin: {}",
                                      SampleInfo[sample].Name, suffix,
                                      SampleInfo[sample].Binning->GetBinName(SampleInfo[sample].LocalId, MaCh3Bin));
          //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
          // https://root.cern.ch/doc/master/classTH1.html#auto-bin
          auto PosteriorHist = std::make_unique<TH1D>(ProjName.c_str(), ProjName.c_str(), 100, 1, -1);
          PosteriorHist->SetDirectory(nullptr);
          PosteriorHist->GetXaxis()->SetTitle("Events");
          int bin = Toys[sample][0]->GetBin(ix, iy);
          for (size_t iToy = 0; iToy < Toys[sample].size(); ++iToy) {
            double content = Toys[sample][iToy]->GetBinContent(bin);
            PosteriorHist->Fill(content, ReweightWeight[iToy]);
          }

          if (DebugHistograms) PosteriorHist->Write();

          PredictiveHist->SetBinContent(ix, iy, PosteriorHist->GetMean());
          PredictiveHist->SetBinError(ix, iy, PosteriorHist->GetRMS());
        }
      }
      PredictiveHist->Write();
      PostPred_mc[sample] = std::move(PredictiveHist);
    }
    else {
      MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, nDims);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  return PostPred_mc;
}

// *************************
// Perform predictive analysis
void PredictiveThrower::RunPredictiveAnalysis() {
// *************************
  MACH3LOG_INFO("Starting {}", __func__);
  TStopwatch TempClock;
  TempClock.Start();

  auto DebugHistograms = GetFromManager<bool>(fitMan->raw()["Predictive"]["DebugHistograms"], false, __FILE__, __LINE__);

  TDirectory* PredictiveDir = outputFile->mkdir("Predictive");
  std::vector<TDirectory*> SampleDirectories;
  SampleDirectories.resize(TotalNumberOfSamples+1);

  // open directory for every sample
  for (int sample = 0; sample < TotalNumberOfSamples+1; ++sample) {
    SampleDirectories[sample] = PredictiveDir->mkdir(SampleInfo[sample].Name.c_str());
  }

  // Produce Violin style spectra
  auto Spectra_mc = ProduceSpectra(MC_Hist_Toy, SampleDirectories, "mc");
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
  StudyBetaParameters(PredictiveDir);

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

// *************************
void PredictiveThrower::PosteriorPredictivepValue(const std::vector<std::unique_ptr<TH1>>& PostPred_mc,
                                                  const std::vector<TDirectory*>& SampleDir) {
// *************************
  //(void) PostPred_w2;
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

    /// TODO This can be multithreaded
    for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
      const int nDims = SampleInfo[iSample].Dimenstion;

      auto DrawFluctHist = M3::Clone(MC_Hist_Toy[iSample][iToy].get());
      auto PredFluctHist = M3::Clone(PostPred_mc[iSample].get());
      if (nDims == 1) {
        MakeFluctuatedHistogramAlternative(static_cast<TH1D*>(DrawFluctHist.get()), static_cast<TH1D*>(MC_Hist_Toy[iSample][iToy].get()), random.get());
        MakeFluctuatedHistogramAlternative(static_cast<TH1D*>(PredFluctHist.get()), static_cast<TH1D*>(PostPred_mc[iSample].get()), random.get());
      }
      else if (nDims == 2) {
        MakeFluctuatedHistogramAlternative(static_cast<TH2D*>(DrawFluctHist.get()), static_cast<TH2D*>(MC_Hist_Toy[iSample][iToy].get()), random.get());
        MakeFluctuatedHistogramAlternative(static_cast<TH2D*>(PredFluctHist.get()), static_cast<TH2D*>(PostPred_mc[iSample].get()), random.get());
      }
      else {
        MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, nDims);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

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
    BetaHist[iSample].resize(SampleInfo[iSample].Binning->GetNBins(SampleInfo[iSample].LocalId));

    for (int iBin = 0; iBin < SampleInfo[iSample].Binning->GetNBins(SampleInfo[iSample].LocalId); iBin++ ) {
      std::string title = fmt::format("#beta param, {} Bin: {}",
                                      SampleInfo[iSample].Name,
                                      SampleInfo[iSample].Binning->GetBinName(SampleInfo[iSample].LocalId, iBin));
      //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
      // https://root.cern.ch/doc/master/classTH1.html#auto-bin
      BetaHist[iSample][iBin] = std::make_unique<TH1D>(title.c_str(), title.c_str(), 100, 1, -1);
      BetaHist[iSample][iBin]->SetDirectory(nullptr);
      BetaHist[iSample][iBin]->GetXaxis()->SetTitle("beta parameter");
    }
  }

  /// 2. Fill histograms, thread safe as all histograms are allocated before and we loop over samples
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    const int nDims = SampleInfo[iSample].Dimenstion;
    if (nDims == 1) {
      int nbinsx = Data_Hist[iSample]->GetNbinsX();
      for (int iBinX = 0; iBinX < nbinsx; ++iBinX) {
        const int RootBin = iBinX+1;
        for (int iToy = 0; iToy < Ntoys; ++iToy) {
          /// ROOT enumerates from 1 while MaCh3 from 0
          const double Data = Data_Hist[iSample]->GetBinContent(RootBin);
          const double MC = MC_Hist_Toy[iSample][iToy]->GetBinContent(RootBin);
          const double w2 = W2_Hist_Toy[iSample][iToy]->GetBinContent(RootBin);
          const auto likelihood = SampleInfo[iSample].SamHandler->GetTestStatistic();

          const double BetaParam = GetBetaParameter(Data, MC, w2, likelihood);
          BetaHist[iSample][iBinX]->Fill(BetaParam, ReweightWeight[iToy]);
        }
      }
    } else if (nDims == 2) {
      const int nX = Data_Hist[iSample]->GetNbinsX();
      const int nY = Data_Hist[iSample]->GetNbinsY();
      for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
          const int RootBin = Data_Hist[iSample]->GetBin(x+1, y+1);
          const int MaCh3Bin = SampleInfo[iSample].Binning->GetBinSafe(SampleInfo[iSample].LocalId, {x, y});

          for (int iToy = 0; iToy < Ntoys; ++iToy) {
            const double Data = Data_Hist[iSample]->GetBinContent(RootBin);
            const double MC   = MC_Hist_Toy[iSample][iToy]->GetBinContent(RootBin);
            const double w2   = W2_Hist_Toy[iSample][iToy]->GetBinContent(RootBin);
            const auto likelihood = SampleInfo[iSample].SamHandler->GetTestStatistic();

            const double BetaParam = GetBetaParameter(Data, MC, w2, likelihood);
            BetaHist[iSample][MaCh3Bin]->Fill(BetaParam, ReweightWeight[iToy]);
          }
        }
      }
    } else {
      MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, nDims);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  /// 3. Write to file
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    for (int iBin = 0; iBin < SampleInfo[iSample].Binning->GetNBins(SampleInfo[iSample].LocalId); iBin++ ) {
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
