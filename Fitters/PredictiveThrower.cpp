#include "PredictiveThrower.h"
#include "Parameters/ParameterHandlerGeneric.h"
#include "TH3.h"

//this file is choc full of usage of a root interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wuseless-cast"

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
void PredictiveThrower::SetParamters(std::vector<std::string>& ParameterGroupsNotVaried,
                                     std::unordered_set<int>& ParameterOnlyToVary) {
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
    TotalNumberOfSamples += samples[iPDF]->GetNSamples();
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
    for (int SampleIndex = 0; SampleIndex < samples[iPDF]->GetNSamples(); ++SampleIndex) {
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
void PredictiveThrower::SetupToyGeneration(std::vector<std::string>& ParameterGroupsNotVaried,
                                           std::unordered_set<int>& ParameterOnlyToVary,
                                           std::vector<const M3::float_t*>& BoundValuePointer,
                                           std::vector<std::pair<double, double>>& ParamBounds) {
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

  auto paramNode = fitMan->raw()["Predictive"]["ParameterBounds"];
  for (const auto& p : paramNode) {
    // Extract name
    std::string name = p[0].as<std::string>();

    // Extract bounds: min and max
    double minVal = p[1][0].as<double>();
    double maxVal = p[1][1].as<double>();
    ParamBounds.emplace_back(minVal, maxVal);

    for (size_t s = 0; s < systematics.size(); ++s) {
      for(int iPar = 0; iPar < systematics[s]->GetNParameters(); iPar++){
        if(systematics[s]->GetParFancyName(iPar) == name){
          BoundValuePointer.push_back(systematics[s]->RetPointer(iPar));
          break;
        }
      }
    }
    if(ParamBounds.size() != BoundValuePointer.size()){
      MACH3LOG_ERROR("Ddin't find paramter {}", name);
      throw MaCh3Exception(__FILE__,__LINE__);
    }
    MACH3LOG_INFO("Parameter: {} with : [{}, {}]", name, minVal, maxVal);
  }
  if(Is_PriorPredictive && ParamBounds.size() > 0) {
    MACH3LOG_ERROR("Additional bounds not supported by prior predictive right now");
    throw MaCh3Exception(__FILE__,__LINE__);
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
    for (int iSample = 0; iSample < SampleHandler->GetNSamples(); ++iSample)
    {
      ToyDirectory->cd();

      auto SampleName = SampleHandler->GetSampleTitle(iSample);
      const TH1* MCHist = SampleHandler->GetMCHist(iSample);
      MC_Hist_Toy[SampleCounter][iToy] = M3::Clone(MCHist, SampleName + "_mc_" + std::to_string(iToy));
      MC_Hist_Toy[SampleCounter][iToy]->Write();

      const TH1* W2Hist = SampleHandler->GetW2Hist(iSample);
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
        }
      }
      SampleCounter++;
    }
  }
}

// *************************
void PredictiveThrower::WriteByModeToys(TDirectory* ByModeDirectory,
                                        const int iToy) {
// *************************
  int SampleCounter = 0;
  for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
  {
    auto* SampleHandler = samples[iPDF];
    auto* modes = SampleHandler->GetMaCh3Modes();
    for (int iSample = 0; iSample < SampleHandler->GetNSamples(); ++iSample)
    {
      ByModeDirectory->cd();

      auto SampleName = SampleHandler->GetSampleTitle(iSample);
      for (int iMode = 0; iMode < modes->GetNModes()+1; ++iMode) {
        auto ModeName = modes->GetMaCh3ModeName(iMode);
        for(int iDim = 0; iDim < SampleHandler->GetNDim(iSample); iDim++) {
          std::string ProjectionName = SampleHandler->GetKinVarName(iSample, iDim);
          std::string PlotSuffix = "_1DProj" + std::to_string(iDim) + "_" + ModeName + "_" + std::to_string(iToy);

          auto hist = SampleHandler->Get1DVarHistByModeAndChannel(iSample, ProjectionName, iMode);
          hist->SetTitle((SampleName + PlotSuffix).c_str());
          hist->SetName((SampleName + PlotSuffix).c_str());
          hist->Write();
        } // end loop over dimension
      } // end loop over mode
      SampleCounter++;
    } // end loop over sample
  } // end loop over sample handler objects
}

// *************************
bool CheckBounds(const std::vector<const M3::float_t*>& BoundValuePointer,
                 const std::vector<std::pair<double,double>>& ParamBounds) {
// *************************
  for (size_t i = 0; i < BoundValuePointer.size(); ++i) {
    const double val = *(BoundValuePointer[i]);
    const double minVal = ParamBounds[i].first;
    const double maxVal = ParamBounds[i].second;

    if (val < minVal || val > maxVal)
      return false; // out of bounds
  }
  return true; // all values are within bounds
}

// *************************
// Produce MaCh3 toys:
void PredictiveThrower::ProduceToys() {
// *************************
  // If we found toys then skip process of making new toys
  if(LoadToys()) return;

  /// KS: Names of parameter groups that will not be varied
  std::vector<std::string> ParameterGroupsNotVaried;
  /// KS: Index of parameters that will be varied
  std::unordered_set<int> ParameterOnlyToVary;
  // For study where one would like to apply bounds
  std::vector<const M3::float_t*> BoundValuePointer;
  std::vector<std::pair<double, double>> ParamBounds;

  // Setup useful information for toy generation
  SetupToyGeneration(ParameterGroupsNotVaried, ParameterOnlyToVary,
                     BoundValuePointer, ParamBounds);

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
  std::vector<const M3::float_t*> ParampPointers(NModelParams);
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
    auto* MaCh3Sample = samples[iPDF];
    for (int SampleIndex = 0; SampleIndex < MaCh3Sample->GetNSamples(); ++SampleIndex)
    {
      // Get nominal spectra and event rates
      const TH1* DataHist = MaCh3Sample->GetDataHist(SampleIndex);
      Data_Hist[SampleCounter] = M3::Clone(DataHist, MaCh3Sample->GetSampleTitle(SampleIndex) + "_data");
      Data_Hist[SampleCounter]->Write((MaCh3Sample->GetSampleTitle(SampleIndex) + "_data").c_str());

      const TH1* MCHist = MaCh3Sample->GetMCHist(SampleIndex);
      MC_Nom_Hist[SampleCounter] = M3::Clone(MCHist, MaCh3Sample->GetSampleTitle(SampleIndex) + "_mc");
      MC_Nom_Hist[SampleCounter]->Write((MaCh3Sample->GetSampleTitle(SampleIndex) + "_mc").c_str());

      const TH1* W2Hist = MaCh3Sample->GetW2Hist(SampleIndex);
      W2_Nom_Hist[SampleCounter] = M3::Clone(W2Hist, MaCh3Sample->GetSampleTitle(SampleIndex) + "_w2");
      W2_Nom_Hist[SampleCounter]->Write((MaCh3Sample->GetSampleTitle(SampleIndex) + "_w2").c_str());
      SampleCounter++;
    }
  }

  TDirectory* Toy_1DDirectory = outputFile->mkdir("Toys_1DHistVar");
  TDirectory* Toy_2DDirectory = outputFile->mkdir("Toys_2DHistVar");
  auto doByMode = GetFromManager<bool>(fitMan->raw()["Predictive"]["ByMode"], false, __FILE__, __LINE__);
  TDirectory* ByModeDirectory = nullptr;
  if(doByMode) ByModeDirectory = outputFile->mkdir("Toys_ByMode");

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
      M3::Utils::PrintProgressBar(i, Ntoys);
    }
    if(!Is_PriorPredictive){
      int entry = 0;
      Step = 0;
      // KS This allow to set additional bounds like mass ordering
      bool WithinBounds = false;
      //YSP: Ensures you get an entry from the mcmc even when burn_in is set to zero (Although not advised :p ).
      //Take 200k burn in steps, WP: Eb C in 1st peaky
      // If we have combined chains by hadd need to check the step in the chain
      // Note, entry is not necessarily same as step due to merged ROOT files, so can't choose entry in the range BurnIn - nEntries :(
      while(Step < burn_in || !WithinBounds) {
        entry = random->Integer(static_cast<unsigned int>(PosteriorFile->GetEntries()));
        PosteriorFile->GetEntry(entry);
        // KS: This might be bit hacky... but BoundValuePointer refer to values in ParameterHandler
        // so we need to update them
        if(BoundValuePointer.size() > 0) {
          for (size_t s = 0; s < systematics.size(); ++s) {
            systematics[s]->SetParameters(branch_vals[s]);
          }
        }
        WithinBounds = CheckBounds(BoundValuePointer, ParamBounds);
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
    SetParamters(ParameterGroupsNotVaried, ParameterOnlyToVary);

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
    if(doByMode) WriteByModeToys(ByModeDirectory, i);

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
  if(doByMode){
    ByModeDirectory->Close();
    delete ByModeDirectory;
  }

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
void PredictiveThrower::StudyByMode1DProjections(const std::vector<TDirectory*>& SampleDirectories) const {
// *************************
  MACH3LOG_INFO("Starting {}", __func__);

  TDirectory * ogdir = gDirectory;
  auto PosteriorFileName = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);
  // Open the ROOT file
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  TFile* file = TFile::Open(PosteriorFileName.c_str(), "READ");

  gErrorIgnoreLevel = originalErrorWarning;
  TDirectory* ToyDir = file->GetDirectory("Toys_ByMode");
  // If toys not amiable in posterior file this means they must be in output file
  if(ToyDir == nullptr) {
    ToyDir = outputFile->GetDirectory("Toys_ByMode");
  }
  /// @todo KS: Here we assume each sample has same modes, this is because ProduceSpectra function,
  /// expects vector [sample], [toy], [dim], so we make ProjectionToys with [mode], [sample], [toy], [dim]
  /// so we can reuse this functionality
  auto* mode = SampleInfo[0].SamHandler->GetMaCh3Modes();
  auto NModes = mode->GetNModes()+1;
  // [mode], [sample], [toy], [dim]
  std::vector<std::vector<std::vector<std::vector<std::unique_ptr<TH1D>>>>> ProjectionToys(NModes);
  for(int iMode = 0; iMode < NModes; iMode++) {
    ProjectionToys[iMode].resize(TotalNumberOfSamples);
    for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
      ProjectionToys[iMode][sample].resize(Ntoys);
      const int nDims = SampleInfo[sample].Dimenstion;
      for (int iToy = 0; iToy < Ntoys; ++iToy) {
        ProjectionToys[iMode][sample][iToy].resize(nDims);
      }
    }
  }

  for (int iToy = 0; iToy < Ntoys; ++iToy) {
    if (iToy % 100 == 0) MACH3LOG_INFO("   Loaded Projection toys {}", iToy);
    for(int iMode = 0; iMode < NModes; iMode++) {
      auto ModeName = mode->GetMaCh3ModeName(iMode);
      for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
        const int nDims = SampleInfo[sample].Dimenstion;
        for(int iDim = 0; iDim < nDims; iDim ++) {
          std::string ProjectionSuffix = "_1DProj" + std::to_string(iDim) + "_" + ModeName + "_" + std::to_string(iToy);
          TH1D* MCHist1D = static_cast<TH1D*>(ToyDir->Get((SampleInfo[sample].Name + ProjectionSuffix).c_str()));
          ProjectionToys[iMode][sample][iToy][iDim] = M3::Clone(MCHist1D);
        }
      }
    } // end loop over samples
  } // end loop over toys

  // ByMode directory
  std::vector<TDirectory*> ModeDirectory(TotalNumberOfSamples);
  for(int iSample = 0; iSample < TotalNumberOfSamples; iSample++) {
    ModeDirectory[iSample] = SampleDirectories[iSample]->mkdir("ByMode");
  }
  // Produce By Mode Spectra
  for(int iMode = 0; iMode < NModes; iMode++) {
    auto ModeName = mode->GetMaCh3ModeName(iMode);
    ProduceSpectra(ProjectionToys[iMode], ModeDirectory, ModeName, false);
  }
  for(int iSample = 0; iSample < TotalNumberOfSamples; iSample++) {
    ModeDirectory[iSample]->Close();
    delete ModeDirectory[iSample];
  }
  file->Close(); delete file;
  if(ogdir){ ogdir->cd(); }
}

// *************************
void PredictiveThrower::ProduceSpectra(const std::vector<std::vector<std::vector<std::unique_ptr<TH1D>>>>& Toys,
                                       const std::vector<TDirectory*>& SampleDirectories,
                                       const std::string suffix,
                                       const bool DoSummary) const {
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
      if(nDims == 2 && DoSummary) {
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
          // KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
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
                                                                    const bool DebugHistograms,
                                                                    const bool WriteHist) {
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
    if(WriteHist) PostPred[sample]->Write();
  } // end loop over samples
  return PostPred;
}

// *************************
// Perform predictive analysis
void PredictiveThrower::RunPredictiveAnalysis() {
// *************************
  // Remove not useful stuff
  SanitiseInputs();

  MACH3LOG_INFO("Starting {}", __func__);
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", M3::Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", M3::Utils::getValue("MemTotal") / 1048576.0);

  TStopwatch TempClock;
  TempClock.Start();

  auto DebugHistograms = GetFromManager<bool>(fitMan->raw()["Predictive"]["DebugHistograms"], false, __FILE__, __LINE__);
  auto doByMode = GetFromManager<bool>(fitMan->raw()["Predictive"]["ByMode"], false, __FILE__, __LINE__);

  TDirectory* PredictiveDir = outputFile->mkdir("Predictive");
  std::vector<TDirectory*> SampleDirectories;
  SampleDirectories.resize(TotalNumberOfSamples+1);

  // open directory for every sample
  for (int sample = 0; sample < TotalNumberOfSamples+1; ++sample) {
    SampleDirectories[sample] = PredictiveDir->mkdir(SampleInfo[sample].Name.c_str());
  }

  // Produce Violin style spectra
  Study1DProjections(SampleDirectories);
  // Produce Post pred by each mode individually
  if(doByMode) StudyByMode1DProjections(SampleDirectories);
  // Produce posterior predictive distribution for mc
  auto PostPred_mc = MakePredictive(MC_Hist_Toy, SampleDirectories, "mc", DebugHistograms, false);
  // Produce posterior predictive distribution for w2
  auto PostPred_w2 = MakePredictive(W2_Hist_Toy, SampleDirectories, "w2", false, false);
  // Calculate Posterior Predictive LLH
  PredictiveLLH(Data_Hist, PostPred_mc, PostPred_w2, SampleDirectories);
  // Calculate Posterior Predictive $p$-value
  PosteriorPredictivepValue(PostPred_mc, SampleDirectories);
  // Check how number of events changed
  RateAnalysis(MC_Hist_Toy, SampleDirectories);

  // Close directories
  for (int sample = 0; sample < TotalNumberOfSamples+1; ++sample) {
    SampleDirectories[sample]->Close();
    delete SampleDirectories[sample];
  }

  auto StudyBeta = GetFromManager<bool>(fitMan->raw()["Predictive"]["StudyBetaParameters"], true, __FILE__, __LINE__);
  auto StudyInfoCriterion = GetFromManager<bool>(fitMan->raw()["Predictive"]["StudyInformationCriterion"], true, __FILE__, __LINE__);
  auto StudyCorr = GetFromManager<bool>(fitMan->raw()["Predictive"]["StudyCorrelations"], true, __FILE__, __LINE__);

  // Studying information criterion
  if(StudyInfoCriterion) StudyInformationCriterion(M3::kWAIC, PostPred_mc, PostPred_w2);
  // Study Prior/Posterior correlations between samples etc.
  if(StudyCorr) StudyCorrelations(PredictiveDir, MC_Hist_Toy, DebugHistograms);
  // Perform beta analysis for mc statical uncertainty
  if(StudyBeta) StudyBetaParameters(PredictiveDir);

  PredictiveDir->Close();
  delete PredictiveDir;

  outputFile->cd();

  TempClock.Stop();
  MACH3LOG_INFO("{} took {:.2f}s to finish for {} toys", __func__, TempClock.RealTime(), Ntoys);
}

// *************************
double PredictiveThrower::CalcLLH(const double data,
                                  const double mc,
                                  const double w2,
                                  const SampleHandlerInterface* SampleHandler) const {
// *************************
  double llh = SampleHandler->GetTestStatLLH(data, mc, w2);
  //KS: do times 2 because banff reports chi2
  return 2*llh;
}

// *************************
double PredictiveThrower::CalcLLH(const TH1* DatHist,
                                  const TH1* MCHist,
                                  const TH1* W2Hist,
                                  const SampleHandlerInterface* SampleHandler) const {
// *************************
  // 1D case
  if (auto h1 = dynamic_cast<const TH1D*>(DatHist)) {
    return GetLLH(h1,
                  static_cast<const TH1D*>(MCHist),
                  static_cast<const TH1D*>(W2Hist),
                  SampleHandler);
  }

  // 2D case
  if (auto h2 = dynamic_cast<const TH2D*>(DatHist)) {
    return GetLLH(h2,
                  static_cast<const TH2D*>(MCHist),
                  static_cast<const TH2D*>(W2Hist),
                  SampleHandler);
  }

  // 2D poly case
  if (auto h2p = dynamic_cast<const TH2Poly*>(DatHist)) {
    return GetLLH(h2p,
                  static_cast<const TH2Poly*>(MCHist),
                  static_cast<const TH2Poly*>(W2Hist),
                  SampleHandler);
  }

  MACH3LOG_ERROR("Unsupported histogram type in {}", __func__);
  throw MaCh3Exception(__FILE__ , __LINE__ );
}

// *************************
double PredictiveThrower::GetLLH(const TH1D* DatHist,
                                 const TH1D* MCHist,
                                 const TH1D* W2Hist,
                                 const SampleHandlerInterface* SampleHandler) const {
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
double PredictiveThrower::GetLLH(const TH2Poly* DatHist,
                                 const TH2Poly* MCHist,
                                 const TH2Poly* W2Hist,
                                 const SampleHandlerInterface* SampleHandler) const {
// *************************
  double llh = 0.0;
  for (int i = 1; i <= DatHist->GetNumberOfBins(); ++i)
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
double PredictiveThrower::GetLLH(const TH2D* DatHist,
                                 const TH2D* MCHist,
                                 const TH2D* W2Hist,
                                 const SampleHandlerInterface* SampleHandler) const {
// *************************
  double llh = 0.0;

  const int nBinsX = DatHist->GetXaxis()->GetNbins();
  const int nBinsY = DatHist->GetYaxis()->GetNbins();

  for (int i = 1; i <= nBinsX; ++i)
  {
    for (int j = 1; j <= nBinsY; ++j)
    {
      const double data = DatHist->GetBinContent(i, j);
      const double mc   = MCHist->GetBinContent(i, j);
      const double w2   = W2Hist->GetBinContent(i, j);

      llh += SampleHandler->GetTestStatLLH(data, mc, w2);
    }
  }

  // KS: do times 2 because banff reports chi2
  return 2 * llh;
}

// ****************
//KS: We have two methods how to apply statistical fluctuation standard is faster hence is default
void PredictiveThrower::MakeFluctuatedHistogram(TH1* FluctHist, TH1* Hist) {
// ****************
  // Determine which fluctuation function to call
  auto applyFluctuation = [&](auto* f, auto* h) {
    if (StandardFluctuation) {
      MakeFluctuatedHistogramStandard(f, h, random.get());
    } else {
      MakeFluctuatedHistogramAlternative(f, h, random.get());
    }
  };

  if (Hist->InheritsFrom(TH2Poly::Class())) {
    applyFluctuation(static_cast<TH2Poly*>(FluctHist), static_cast<TH2Poly*>(Hist));
  }
  else if (Hist->InheritsFrom(TH2D::Class())) {
    applyFluctuation(static_cast<TH2D*>(FluctHist), static_cast<TH2D*>(Hist));
  }
  else if (Hist->InheritsFrom(TH1D::Class())) {
    applyFluctuation(static_cast<TH1D*>(FluctHist), static_cast<TH1D*>(Hist));
  }
  else {
    MACH3LOG_ERROR("Unsupported histogram type");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
}

// *************************
void PredictiveThrower::PosteriorPredictivepValue(const std::vector<std::unique_ptr<TH1>>& PostPred_mc,
                                                  const std::vector<TDirectory*>& SampleDir) {
// *************************
  // Step 1: Initialize per-toy accumulators once
  // [Sample] [Toys]
  auto make_matrix = [&](double init = 0.0) {
    return std::vector<std::vector<double>>(
      TotalNumberOfSamples + 1,
      std::vector<double>(Ntoys, init));
  };
  auto chi2_dat       = make_matrix();
  auto chi2_mc        = make_matrix();
  auto chi2_pred      = make_matrix();
  auto chi2_rate_dat  = make_matrix();
  auto chi2_rate_mc   = make_matrix();
  auto chi2_rate_pred = make_matrix();

  // 2. Add penalty terms to global bin
  for (int iToy = 0; iToy < Ntoys; ++iToy) {
    chi2_dat[TotalNumberOfSamples][iToy] = PenaltyTerm[iToy];
    chi2_mc[TotalNumberOfSamples][iToy] = PenaltyTerm[iToy];
    chi2_pred[TotalNumberOfSamples][iToy] = PenaltyTerm[iToy];

    chi2_rate_dat[TotalNumberOfSamples][iToy] = PenaltyTerm[iToy];
    chi2_rate_mc[TotalNumberOfSamples][iToy] = PenaltyTerm[iToy];
    chi2_rate_pred[TotalNumberOfSamples][iToy] = PenaltyTerm[iToy];
  }

  /// TODO This can be multithreaded but be careful for Clone!!!
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    auto SampleHandler = SampleInfo[iSample].SamHandler;
    for (int iToy = 0; iToy < Ntoys; ++iToy) {
      // Clone histograms to avoid modifying originals
      auto DrawFluctHist = M3::Clone(MC_Hist_Toy[iSample][iToy].get());
      auto PredFluctHist = M3::Clone(PostPred_mc[iSample].get());

      // Apply fluctuations
      MakeFluctuatedHistogram(DrawFluctHist.get(), MC_Hist_Toy[iSample][iToy].get());
      MakeFluctuatedHistogram(PredFluctHist.get(), PostPred_mc[iSample].get());

      // I. SHAPE + RATE (bin-by-bin likelihood)
      chi2_dat[iSample][iToy]  = CalcLLH(Data_Hist[iSample].get(), MC_Hist_Toy[iSample][iToy].get(), W2_Hist_Toy[iSample][iToy].get(), SampleHandler);
      chi2_mc[iSample][iToy]   = CalcLLH(DrawFluctHist.get(), MC_Hist_Toy[iSample][iToy].get(), W2_Hist_Toy[iSample][iToy].get(), SampleHandler);
      chi2_pred[iSample][iToy] = CalcLLH(PredFluctHist.get(), MC_Hist_Toy[iSample][iToy].get(), W2_Hist_Toy[iSample][iToy].get(), SampleHandler);

      // II. RATE-ONLY (total normalization)
      chi2_rate_dat[iSample][iToy]  = CalcLLH(Data_Hist[iSample]->Integral(), MC_Hist_Toy[iSample][iToy]->Integral(), W2_Hist_Toy[iSample][iToy]->Integral(), SampleHandler);
      chi2_rate_mc[iSample][iToy]   = CalcLLH(DrawFluctHist->Integral(), MC_Hist_Toy[iSample][iToy]->Integral(), W2_Hist_Toy[iSample][iToy]->Integral(), SampleHandler);
      chi2_rate_pred[iSample][iToy] = CalcLLH(PredFluctHist->Integral(), MC_Hist_Toy[iSample][iToy]->Integral(), W2_Hist_Toy[iSample][iToy]->Integral(), SampleHandler);

      // III. accumulate global sums ---
      chi2_dat[TotalNumberOfSamples][iToy]  += chi2_dat[iSample][iToy];
      chi2_mc[TotalNumberOfSamples][iToy]   += chi2_mc[iSample][iToy];
      chi2_pred[TotalNumberOfSamples][iToy] += chi2_pred[iSample][iToy];

      chi2_rate_dat[TotalNumberOfSamples][iToy]  += chi2_rate_dat[iSample][iToy];
      chi2_rate_mc[TotalNumberOfSamples][iToy]   += chi2_rate_mc[iSample][iToy];
      chi2_rate_pred[TotalNumberOfSamples][iToy] += chi2_rate_pred[iSample][iToy];
    }
  }

  // 4. Produce pValue plots
  // Shape+rate posterior predictive checks
  MakeChi2Plots(chi2_mc,   "-2LLH (Draw Fluc, Draw)", chi2_dat, "-2LLH (Data, Draw)", SampleDir, "_drawfluc_draw");
  MakeChi2Plots(chi2_pred, "-2LLH (Pred Fluc, Draw)", chi2_dat, "-2LLH (Data, Draw)", SampleDir, "_predfluc_draw");

  // Rate-only posterior predictive checks
  MakeChi2Plots(chi2_rate_mc,   "-2LLH (Rate Draw Fluc, Draw)", chi2_rate_dat, "-2LLH (Rate Data, Draw)", SampleDir, "_rate_drawfluc_draw");
  MakeChi2Plots(chi2_rate_pred, "-2LLH (Rate Pred Fluc, Draw)", chi2_rate_dat, "-2LLH (Rate Data, Draw)", SampleDir, "_rate_predfluc_draw");
}

// *************************
void PredictiveThrower::PredictiveLLH(const std::vector<std::unique_ptr<TH1>>& Data_histogram,
                                               const std::vector<std::unique_ptr<TH1>>& PostPred_mc,
                                               const std::vector<std::unique_ptr<TH1>>& PostPred_w,
                                               const std::vector<TDirectory*>& SampleDir) {
// *************************
  MACH3LOG_INFO("{:<55} {:<10} {:<10} {:<10}", "Sample", "DataInt", "MCInt", "-2LLH");
  MACH3LOG_INFO("{:-<55} {:-<10} {:-<10} {:-<10}", "", "", "", "");
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    SampleDir[iSample]->cd();
    ExtractLLH(Data_histogram[iSample].get(), PostPred_mc[iSample].get(), PostPred_w[iSample].get(), SampleInfo[iSample].SamHandler);
    PostPred_mc[iSample]->Write();
  }
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
      chi2_y_sample[iToy] = Chi2_y[iSample][iToy];
      chi2_x_per_sample[iToy]  = Chi2_x[iSample][iToy];
    }

    const double min_val = std::min(*std::min_element(chi2_y_sample.begin(), chi2_y_sample.end()),
                              *std::min_element(chi2_x_per_sample.begin(), chi2_x_per_sample.end()));
    const double max_val = std::max(*std::max_element(chi2_y_sample.begin(), chi2_y_sample.end()),
                              *std::max_element(chi2_x_per_sample.begin(), chi2_x_per_sample.end()));

    auto chi2_hist = std::make_unique<TH2D>((SampleInfo[iSample].Name+ Title).c_str(),
                                            (SampleInfo[iSample].Name+ Title).c_str(),
                                            50, min_val, max_val, 50, min_val, max_val);
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

// ****************
// Study Prior/Posterior correlations between samples etc.
void PredictiveThrower::StudyCorrelations(TDirectory* PredictiveDir,
                                          const std::vector<std::vector<std::unique_ptr<TH1>>>& Toys,
                                          const bool DebugHistograms) const {
// ****************
  MACH3LOG_INFO("Startin {}", __func__);

  // Make a new directory
  TDirectory *CorrDir = PredictiveDir->mkdir("Correlations");
  CorrDir->cd();

  std::vector<double> minVals(TotalNumberOfSamples, std::numeric_limits<double>::max());
  std::vector<double> maxVals(TotalNumberOfSamples, std::numeric_limits<double>::lowest());
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < TotalNumberOfSamples; ++i)
  {
    for (const auto& toyHist : Toys[i])
    {
      const double val = toyHist->Integral();
      if (val < minVals[i]) minVals[i] = val;
      if (val > maxVals[i]) maxVals[i] = val;
    }
  }
  auto hSamCorr = std::make_unique<TH2D>("Sample Correlation", "Sample Correlation", TotalNumberOfSamples, 0,
                                         TotalNumberOfSamples, TotalNumberOfSamples, 0, TotalNumberOfSamples);
  hSamCorr->SetDirectory(nullptr);
  hSamCorr->GetZaxis()->SetTitle("Correlation");
  hSamCorr->SetMinimum(-1);
  hSamCorr->SetMaximum(1);
  hSamCorr->GetXaxis()->SetLabelSize(0.015);
  hSamCorr->GetYaxis()->SetLabelSize(0.015);
  // Loop over the Covariance matrix entries
  for (int i = 0; i < TotalNumberOfSamples; ++i) {
    hSamCorr->SetBinContent(i+1, i+1, 1.0);
    hSamCorr->GetXaxis()->SetBinLabel(i+1, SampleInfo[i].Name.c_str());
    for (int j = 0; j < TotalNumberOfSamples; ++j) {
      hSamCorr->GetYaxis()->SetBinLabel(j+1, SampleInfo[j].Name.c_str());
    }
  }

  std::vector<std::vector<std::unique_ptr<TH2D>>> SamCorr(TotalNumberOfSamples);
  for (int i = 0; i < TotalNumberOfSamples; ++i)
  {
    SamCorr[i].resize(TotalNumberOfSamples);
    const double Min_i = minVals[i];
    const double Max_i = maxVals[i];
    for (int j = 0; j < TotalNumberOfSamples; ++j)
    {
      const double Min_j = minVals[j];
      const double Max_j = maxVals[j];
      // TH2D to hold the Correlation
      std::string name  = "SamCorr_" + std::to_string(i) + "_" + std::to_string(j);
      SamCorr[i][j] = std::make_unique<TH2D>(name.c_str(), name.c_str(), 70, Min_i, Max_i, 70, Min_j, Max_j);
      SamCorr[i][j]->SetDirectory(nullptr);
      SamCorr[i][j]->SetMinimum(0);
      SamCorr[i][j]->GetXaxis()->SetTitle(SampleInfo[i].Name.c_str());
      SamCorr[i][j]->GetYaxis()->SetTitle(SampleInfo[j].Name.c_str());
      SamCorr[i][j]->GetZaxis()->SetTitle("Events");
    }
  }

  // Now we are sure we have the diagonal elements, let's make the off-diagonals
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < TotalNumberOfSamples; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      // Skip the diagonal elements which we've already done above
      if (j == i) continue;

      for (int iToy = 0; iToy < Ntoys; ++iToy)
      {
        SamCorr[i][j]->Fill(Toys[i][iToy]->Integral(), Toys[j][iToy]->Integral());
      }
      SamCorr[i][j]->Smooth();

      // The value of the Covariance
      const double corr = SamCorr[i][j]->GetCorrelationFactor();
      hSamCorr->SetBinContent(i+1, j+1, corr);
      hSamCorr->SetBinContent(j+1, i+1, corr);
    }// End j loop
  }// End i loop

  hSamCorr->Draw("colz");
  hSamCorr->Write("Sample_Corr");

  if(DebugHistograms) {
    for (int i = 0; i < TotalNumberOfSamples; ++i){
      for (int j = 0; j <= i; ++j) {
        // Skip the diagonal elements which we've already done above
        if (j == i) continue;
        SamCorr[i][j]->Write();
      }// End j loop
    }// End i loop
  } // end if debugHist

  PredictiveDir->cd();
}

// ****************
// Calculate the LLH for TH1, set the LLH to title of MCHist
void PredictiveThrower::ExtractLLH(TH1*  DatHist, TH1* MCHist, TH1* W2Hist, const SampleHandlerInterface* SampleHandler) const {
// ****************
  const double llh = CalcLLH(DatHist, MCHist, W2Hist, SampleHandler);
  std::stringstream ss;
  ss << "_2LLH=" << llh;
  MCHist->SetTitle((std::string(MCHist->GetTitle())+ss.str()).c_str());
  MACH3LOG_INFO("{:<55} {:<10.2f} {:<10.2f} {:<10.2f}", MCHist->GetName(), DatHist->Integral(), MCHist->Integral(), llh);
}

// ****************
// Make the 1D Event Rate Hist
void PredictiveThrower::MakeCutEventRate(TH1D *Histogram, const double DataRate) const {
// ****************
  // Open the ROOT file
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  // For the event rate histogram add a TLine to the data rate
  auto TempLine = std::make_unique<TLine>(DataRate, Histogram->GetMinimum(), DataRate, Histogram->GetMaximum());
  TempLine->SetLineColor(kRed);
  TempLine->SetLineWidth(2);
  // Also fit a Gaussian because why not?
  auto Fitter = std::make_unique<TF1>("Fit", "gaus", Histogram->GetBinLowEdge(1), Histogram->GetBinLowEdge(Histogram->GetNbinsX()+1));
  Histogram->Fit(Fitter.get(), "RQ");
  Fitter->SetLineColor(kRed-5);
  // Calculate a p-value
  double Above = 0.0;
  for (int z = 0; z < Histogram->GetNbinsX(); ++z) {
    const double xvalue = Histogram->GetBinCenter(z+1);
    if (xvalue >= DataRate) {
      Above += Histogram->GetBinContent(z+1);
    }
  }
  const double pvalue = Above/Histogram->Integral();
  TLegend Legend(0.4, 0.75, 0.98, 0.90);
  Legend.SetFillColor(0);
  Legend.SetFillStyle(0);
  Legend.SetLineWidth(0);
  Legend.SetLineColor(0);
  Legend.AddEntry(TempLine.get(), Form("Data, %.0f, p-value=%.2f", DataRate, pvalue), "l");
  Legend.AddEntry(Histogram, Form("MC, #mu=%.1f#pm%.1f", Histogram->GetMean(), Histogram->GetRMS()), "l");
  Legend.AddEntry(Fitter.get(), Form("Gauss, #mu=%.1f#pm%.1f", Fitter->GetParameter(1), Fitter->GetParameter(2)), "l");
  std::string TempTitle = std::string(Histogram->GetName());
  TempTitle += "_canv";
  TCanvas TempCanvas(TempTitle.c_str(), TempTitle.c_str(), 1024, 1024);
  TempCanvas.SetGridx();
  TempCanvas.SetGridy();
  TempCanvas.SetRightMargin(0.03);
  TempCanvas.SetBottomMargin(0.08);
  TempCanvas.SetLeftMargin(0.10);
  TempCanvas.SetTopMargin(0.06);
  TempCanvas.cd();
  Histogram->Draw();
  TempLine->Draw("same");
  Fitter->Draw("same");
  Legend.Draw("same");
  TempCanvas.Write();
  Histogram->Write();
  gErrorIgnoreLevel = originalErrorWarning;
}

// *************************
void PredictiveThrower::RateAnalysis(const std::vector<std::vector<std::unique_ptr<TH1>>>& Toys,
                                     const std::vector<TDirectory*>& SampleDirectories) const {
// *************************
  std::vector<std::unique_ptr<TH1D>> EventHist(TotalNumberOfSamples+1);
  for (int iSample = 0; iSample < TotalNumberOfSamples+1; ++iSample) {
    std::string Title = "EventHist: ";
    if (iSample == TotalNumberOfSamples) {
      Title = "Total";
    } else {
      Title = SampleInfo[iSample].Name;
    }
    Title += "_sum";
    //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
    // https://root.cern.ch/doc/master/classTH1.html#auto-bin
    EventHist[iSample] = std::make_unique<TH1D>(Title.c_str(), Title.c_str(), 100, 1, -1);
    EventHist[iSample]->SetDirectory(nullptr);
    EventHist[iSample]->GetXaxis()->SetTitle("Total event rate");
    EventHist[iSample]->GetYaxis()->SetTitle("Counts");
    EventHist[iSample]->SetLineWidth(2);
  }

  // First fill per-sample histograms
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    for (int iToy = 0; iToy < Ntoys; ++iToy) {
      double Count = Toys[iSample][iToy]->Integral();
      EventHist[iSample]->Fill(Count);
    }
  }

  // Now fill total histogram properly (per toy)
  for (int iToy = 0; iToy < Ntoys; ++iToy) {
    double TotalCount = 0.0;
    for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
      TotalCount += Toys[iSample][iToy]->Integral();
    }
    EventHist[TotalNumberOfSamples]->Fill(TotalCount);
  }

  double DataRate = 0.0;
  std::vector<double> DataRates(TotalNumberOfSamples+1);
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:DataRate)
  #endif
  for (int i = 0; i < TotalNumberOfSamples; ++i) {
    DataRates[i] = Data_Hist[i]->Integral();
    DataRate += DataRates[i];
  }
  DataRates[TotalNumberOfSamples] = DataRate;

  for (int SampleNum = 0; SampleNum < TotalNumberOfSamples+1; ++SampleNum) {
    SampleDirectories[SampleNum]->cd();
    //Make fancy event rate histogram
    MakeCutEventRate(EventHist[SampleNum].get(), DataRates[SampleNum]);
  }
}


// ****************
void PredictiveThrower::StudyInformationCriterion(M3::kInfCrit Criterion,
                                                  const std::vector<std::unique_ptr<TH1>>& PostPred_mc,
                                                  const std::vector<std::unique_ptr<TH1>>& PostPred_w) {
// ****************
  MACH3LOG_INFO("******************************");
  switch(Criterion) {
    case M3::kInfCrit::kBIC:
      // Study Bayesian Information Criterion
      StudyBIC(PostPred_mc, PostPred_w);
      break;
    case M3::kInfCrit::kDIC:
      // Study Deviance Information Criterion
      StudyDIC(PostPred_mc, PostPred_w);
      break;
    case M3::kInfCrit::kWAIC:
      // Study Watanabe-Akaike information criterion (WAIC)
      StudyWAIC();
      break;
    case M3::kInfCrit::kInfCrits:
      MACH3LOG_ERROR("kInfCrits is not a valid kInfCrit!");
      throw MaCh3Exception(__FILE__, __LINE__);
    default:
      MACH3LOG_ERROR("UNKNOWN Information Criterion SPECIFIED!");
      MACH3LOG_ERROR("You gave {}", static_cast<int>(Criterion));
      throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  MACH3LOG_INFO("******************************");
}

// ****************
void PredictiveThrower::StudyBIC(const std::vector<std::unique_ptr<TH1>>& PostPred_mc,
                                 const std::vector<std::unique_ptr<TH1>>& PostPred_w) {
// ****************
  //make fancy event rate histogram
  double DataRate = 0.0;
  double BinsRate = 0.0;
  double TotalLLH = 0.0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:DataRate, BinsRate, TotalLLH)
  #endif
  for (int i = 0; i < TotalNumberOfSamples; ++i)
  {
    auto SampleHandler = SampleInfo[i].SamHandler;
    auto* h = Data_Hist[i].get();
    DataRate += h->Integral();
    if (auto h1 = dynamic_cast<TH1D*>(h)) {
      BinsRate += h1->GetNbinsX();
    } else if (auto h2 = dynamic_cast<TH2D*>(h)) {
      BinsRate += h2->GetNbinsX() * h2->GetNbinsY();
    } else if (auto h2poly = dynamic_cast<TH2Poly*>(h)) {
      BinsRate += h2poly->GetNumberOfBins();
    } else {
      MACH3LOG_WARN("Unknown histogram type in DataHist[{}]", i);
    }
    TotalLLH += CalcLLH(Data_Hist[i].get(), PostPred_mc[i].get(), PostPred_w[i].get(), SampleHandler);
  }

  const double EventRateBIC = GetBIC(TotalLLH, DataRate, NModelParams);
  const double BinBasedBIC = GetBIC(TotalLLH, BinsRate, NModelParams);
  MACH3LOG_INFO("Calculated Bayesian Information Criterion using global number of events: {:.2f}", EventRateBIC);
  MACH3LOG_INFO("Calculated Bayesian Information Criterion using global number of bins: {:.2f}", BinBasedBIC);
  MACH3LOG_INFO("Additional info: NModelParams: {}, DataRate: {:.2f}, BinsRate: {:.2f}", NModelParams, DataRate, BinsRate);
}

// ****************
// Get the Deviance Information Criterion (DIC)
void PredictiveThrower::StudyDIC(const std::vector<std::unique_ptr<TH1>>& PostPred_mc,
                                 const std::vector<std::unique_ptr<TH1>>& PostPred_w) {
// ****************
  //The posterior mean of the deviance
  double Dbar = 0.;
  double TotalLLH = 0.0;

  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:Dbar)
  #endif
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample)
  {
    auto SampleHandler = SampleInfo[iSample].SamHandler;
    TotalLLH += CalcLLH(Data_Hist[iSample].get(), PostPred_mc[iSample].get(), PostPred_w[iSample].get(), SampleHandler);
    double LLH_temp = 0.;
    for (int iToy = 0; iToy < Ntoys; ++iToy)
    {
      LLH_temp += CalcLLH(Data_Hist[iSample].get(), MC_Hist_Toy[iSample][iToy].get(), W2_Hist_Toy[iSample][iToy].get(), SampleHandler);
    }
    Dbar += LLH_temp;
  }
  Dbar = Dbar / Ntoys;

  // A point estimate of the deviance
  const double Dhat = TotalLLH;

  //Effective number of parameters
  const double p_D = std::fabs(Dbar - Dhat);

  //Actual test stat
  const double DIC_stat = Dhat + 2 * p_D;
  MACH3LOG_INFO("Effective number of parameters following DIC formalism is equal to: {:.2f}", p_D);
  MACH3LOG_INFO("DIC test statistic = {:.2f}", DIC_stat);
}


// ****************
// Helper: update WAIC accumulators for a single toy/bin
void AccumulateWAICToy(const double neg_LLH_temp,
                       double& mean_llh,
                       double& mean_llh_squared,
                       double& sum_exp_llh) {
// ****************
  // Negate the negative log-likelihood to get the actual log-likelihood
  double LLH_temp = -neg_LLH_temp;

  mean_llh += LLH_temp;
  mean_llh_squared += LLH_temp * LLH_temp;
  sum_exp_llh += std::exp(LLH_temp);
}

// ****************
// Helper function to finalize WAIC contributions for one bin
void AccumulateWAICBin(double& mean_llh, double& mean_llh_squared, double& sum_exp_llh,
                       const unsigned int Ntoys, double& lppd, double& p_WAIC) {
// ****************
  // Compute the mean log-likelihood and the squared mean
  mean_llh /= Ntoys;
  mean_llh_squared /= Ntoys;
  sum_exp_llh /= Ntoys;
  sum_exp_llh = std::log(sum_exp_llh);

  // Log pointwise predictive density based on Eq. 4 in Gelman2014
  lppd += sum_exp_llh;

  // Compute the effective number of parameters for WAIC
  p_WAIC += mean_llh_squared - (mean_llh * mean_llh);
}

// ****************
// Get the Watanabe-Akaike information criterion (WAIC)
void PredictiveThrower::StudyWAIC() {
// ****************
  // log pointwise predictive density
  double lppd = 0.;
  // effective number of parameters
  double p_WAIC = 0.;

  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:lppd, p_WAIC)
  #endif
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
    auto SampleHandler = SampleInfo[iSample].SamHandler;
    auto* hData = Data_Hist[iSample].get();

    if (auto h2poly = dynamic_cast<TH2Poly*>(hData)) {
      // TH2Poly: irregular bins, linear indexing
      for (int i = 1; i <= h2poly->GetNumberOfBins(); ++i) {
        const double data = Data_Hist[iSample]->GetBinContent(i);
        double mean_llh = 0.;
        double sum_exp_llh = 0;
        double mean_llh_squared = 0.;

        for (int iToy = 0; iToy < Ntoys; ++iToy) {
          const double mc = MC_Hist_Toy[iSample][iToy]->GetBinContent(i);
          const double w2 = W2_Hist_Toy[iSample][iToy]->GetBinContent(i);
          // Get the -log-likelihood for this sample and bin
          double neg_LLH_temp = SampleHandler->GetTestStatLLH(data, mc, w2);
          AccumulateWAICToy(neg_LLH_temp, mean_llh, mean_llh_squared, sum_exp_llh);
        }
        AccumulateWAICBin(mean_llh, mean_llh_squared, sum_exp_llh, Ntoys, lppd, p_WAIC);
      }
    } else if (auto h2 = dynamic_cast<TH2D*>(hData)) {
      // TH2D: nested loops over X and Y
      for (int ix = 1; ix <= h2->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h2->GetNbinsY(); ++iy) {
          const double data = hData->GetBinContent(ix, iy);
          double mean_llh = 0.;
          double mean_llh_squared = 0.;
          double sum_exp_llh = 0.;
          for (int iToy = 0; iToy < Ntoys; ++iToy) {
            const double mc = MC_Hist_Toy[iSample][iToy]->GetBinContent(ix, iy);
            const double w2 = W2_Hist_Toy[iSample][iToy]->GetBinContent(ix, iy);
            // Get the -log-likelihood for this sample and bin
            double neg_LLH_temp = SampleHandler->GetTestStatLLH(data, mc, w2);
            AccumulateWAICToy(neg_LLH_temp, mean_llh, mean_llh_squared, sum_exp_llh);
          }
          AccumulateWAICBin(mean_llh, mean_llh_squared, sum_exp_llh, Ntoys, lppd, p_WAIC);
        }
      }
    } else if (auto h1 = dynamic_cast<TH1D*>(hData)) {
      // TH1D: 1D histogram
      for (int iBin = 1; iBin <= h1->GetNbinsX(); ++iBin) {
        const double data = hData->GetBinContent(iBin);
        double mean_llh = 0.;
        double mean_llh_squared = 0.;
        double sum_exp_llh = 0.;
        for (int iToy = 0; iToy < Ntoys; ++iToy) {
          const double mc = MC_Hist_Toy[iSample][iToy]->GetBinContent(iBin);
          const double w2 = W2_Hist_Toy[iSample][iToy]->GetBinContent(iBin);

          // Get the -log-likelihood for this sample and bin
          double neg_LLH_temp = SampleHandler->GetTestStatLLH(data, mc, w2);
          AccumulateWAICToy(neg_LLH_temp, mean_llh, mean_llh_squared, sum_exp_llh);
        }
        AccumulateWAICBin(mean_llh, mean_llh_squared, sum_exp_llh, Ntoys, lppd, p_WAIC);
      }
    }
  }

  // Compute WAIC, see Eq. 13 in Gelman2014
  double WAIC = -2 * (lppd - p_WAIC);
  MACH3LOG_INFO("Effective number of parameters following WAIC formalism is equal to: {:.2f}", p_WAIC);
  MACH3LOG_INFO("WAIC = {:.2f}", WAIC);
}
