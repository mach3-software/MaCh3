#include "PredictiveThrower.h"
#include "Samples/SampleHandlerFD.h"
#include "Parameters/ParameterHandlerGeneric.h"

// *************************
PredictiveThrower::PredictiveThrower(manager *man) : FitterBase(man) {
// *************************
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

  if(ModelSystematic && ParameterGroupsNotVaried.size() > 0) ModelSystematic->SetGroupOnlyParameters(ParameterGroupsNotVaried);

  /// Alternatively vary only selected params
  if (ModelSystematic && !ParameterOnlyToVary.empty()) {
    for (int i = 0; i < ModelSystematic->GetNumParams(); ++i) {
      if (ParameterOnlyToVary.find(i) == ParameterOnlyToVary.end()) {
        ModelSystematic->SetParProp(i, ModelSystematic->GetParInit(i));
      }
    }
  }
}

// *************************
// Produce MaCh3 toys:
void PredictiveThrower::SetupToyGeneration() {
// *************************
  int counter = 0;
  TotalNumberOfSamples = 0;

  for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
  {
    auto* MaCh3Sample = dynamic_cast<SampleHandlerFD*>(samples[iPDF]);
    if (!MaCh3Sample) {
      MACH3LOG_ERROR("Sample {} do not inherit from SampleHandlerFD this is not implemented", samples[iPDF]->GetTitle());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    TotalNumberOfSamples += samples[iPDF]->GetNsamples();
    if(samples[iPDF]->GetNsamples() > 1){
      MACH3LOG_ERROR("Sample has more than one sample {} ::", samples[iPDF]->GetNsamples());
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }
  for (size_t s = 0; s < systematics.size(); ++s) {
    auto* MaCh3Params = dynamic_cast<ParameterHandlerGeneric*>(systematics[s]);
    if(MaCh3Params) {
      ModelSystematic = MaCh3Params;
      counter++;
    }
  }
  MC_Hist_Toy.resize(TotalNumberOfSamples);
  W2_Hist_Toy.resize(TotalNumberOfSamples);
  Data_Hist.resize(TotalNumberOfSamples);
  SampleObjectMap.resize(TotalNumberOfSamples);
  SampleNames.resize(TotalNumberOfSamples+1);

  int currentIndex = 0;
  for (size_t iPDF = 0; iPDF < samples.size(); ++iPDF) {
    for (int subSampleIndex = 0; subSampleIndex < samples[iPDF]->GetNsamples(); ++subSampleIndex) {
      SampleObjectMap[currentIndex] = static_cast<int>(iPDF); // map the current global sample index to this sample object
      ++currentIndex;
    }
  }

  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    MC_Hist_Toy[sample].resize(Ntoys);
    W2_Hist_Toy[sample].resize(Ntoys);
    SampleNames[sample] = samples[sample]->GetTitle();
  }
  SampleNames[TotalNumberOfSamples] = "Total";

  if(Is_PriorPredictive) {
    MACH3LOG_INFO("You've chosen to run Prior Predictive Distribution");
  } else {
    auto PosteriorFileName = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);
    //KS: We use MCMCProcessor to get names of covariances that were actually used to produce given chain
    MCMCProcessor Processor(PosteriorFileName);
    Processor.Initialise();

    ///Let's ask the manager what are the file with covariance matrix
    YAML::Node ConfigInChain = Processor.GetCovConfig(kXSecPar);
    if(ModelSystematic){
      YAML::Node ConfigNow = ModelSystematic->GetConfig();
      if (!compareYAMLNodes(ConfigNow, ConfigInChain))
      {
        MACH3LOG_ERROR("Yaml configs in previous chain and current one are different", PosteriorFileName);
        throw MaCh3Exception(__FILE__ , __LINE__ );
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
// Produce MaCh3 toys:
void PredictiveThrower::ProduceToys() {
// *************************
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

  TDirectory* ToyDirectory = outputFile->mkdir("Toys");
  ToyDirectory->cd();

  for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
  {
    auto* MaCh3Sample = dynamic_cast<SampleHandlerFD*>(samples[iPDF]);
    // Get nominal spectra and event rates
    TH1D* DataHist1D = static_cast<TH1D*>(MaCh3Sample->GetDataHist(1));
    Data_Hist[iPDF] = M3::Clone(DataHist1D, MaCh3Sample->GetTitle() + "_data");
    Data_Hist[iPDF]->Write((MaCh3Sample->GetTitle() + "_data").c_str());

    TH1D* MCHist1D = static_cast<TH1D*>(MaCh3Sample->GetMCHist(1));
    MCHist1D->Write((MaCh3Sample->GetTitle() + "_mc").c_str());

    TH1D* W2Hist1D = static_cast<TH1D*>(MaCh3Sample->GetW2Hist(1));
    W2Hist1D->Write((MaCh3Sample->GetTitle() + "_w2").c_str());
    delete W2Hist1D;
  }

  /// this store value of parameters sampled from a chain
  std::vector<std::vector<double>> branch_vals(systematics.size());
  std::vector<std::vector<std::string>> branch_name(systematics.size());

  TChain* PosteriorFile = new TChain("posteriors");
  PosteriorFile->Add(PosteriorFileName.c_str());
  int Step = 0;
  PosteriorFile->SetBranchAddress("step", &Step);

  if (PosteriorFile->GetBranch("Weight")) {
    PosteriorFile->SetBranchStatus("Weight", true);
    PosteriorFile->SetBranchAddress("Weight", &Weight);
  } else {
    MACH3LOG_WARN("Not applying reweighting weight");
    Weight = 1.0;
  }

  for (size_t s = 0; s < systematics.size(); ++s) {
    systematics[s]->MatchMaCh3OutputBranches(PosteriorFile, branch_vals[s], branch_name[s]);
  }
  //Get the burn-in from the config
  auto burn_in = Get<int>(fitMan->raw()["Predictive"]["BurnInSteps"], __FILE__, __LINE__);

  TStopwatch TempClock;
  TempClock.Start();
  for(int i = 0; i < Ntoys; i++)
  {
    if( i % (Ntoys/10) == 0) {
      MaCh3Utils::PrintProgressBar(i, Ntoys);
    }
    int entry = 0;
    Step = -999;

    //YSP: Ensures you get an entry from the mcmc even when burn_in is set to zero (Although not advised :p ).
    //Take 200k burn in steps, WP: Eb C in 1st peaky
    // If we have combined chains by hadd need to check the step in the chain
    // Note, entry is not necessarily same as step due to merged ROOT files, so can't choose entry in the range BurnIn - nEntries :(
    while(Step < burn_in){
      entry = random->Integer(static_cast<unsigned int>(PosteriorFile->GetEntries()));
      PosteriorFile->GetEntry(entry);
    }
    if(!Is_PriorPredictive) Draw = entry;
    for (size_t s = 0; s < systematics.size(); ++s)
    {
      systematics[s]->SetParameters(branch_vals[s]);

      //KS: Below line can help you get prior predictive distributions which are helpful for getting pre and post ND fit spectra
      //YSP: If not set in the config, the code runs SK Posterior Predictive distributions by default. If true, then the code runs SK prior predictive.
      if(Is_PriorPredictive) systematics[s]->ThrowParameters();
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

    for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
    {
      auto* MaCh3Sample = dynamic_cast<SampleHandlerFD*>(samples[iPDF]);

      TH1D* MCHist1D = static_cast<TH1D*>(MaCh3Sample->GetMCHist(1));
      MC_Hist_Toy[iPDF][i] = M3::Clone(MCHist1D, MaCh3Sample->GetTitle() + "_mc_" + std::to_string(i));
      MC_Hist_Toy[iPDF][i]->Write();

      TH1D* W2Hist1D = static_cast<TH1D*>(MaCh3Sample->GetW2Hist(1));
      W2_Hist_Toy[iPDF][i] = M3::Clone(W2Hist1D, MaCh3Sample->GetTitle() + "_w2_" + std::to_string(i));
      W2_Hist_Toy[iPDF][i]->Write();
      delete W2Hist1D;
    }
    ToyTree->Fill();
  }//end of toys loop
  TempClock.Stop();

  delete PosteriorFile;
  ToyDirectory->Close();
  delete ToyDirectory;

  outputFile->cd();
  ToyTree->Write();
  delete ToyTree;

  MACH3LOG_INFO("{} took {:.2f}s to finish for {} toys", __func__, TempClock.RealTime(), Ntoys);
}

// *************************
std::vector<std::unique_ptr<TH2D>> PredictiveThrower::ProduceSpectra(const std::vector<std::vector<std::unique_ptr<TH1D>>>& Toys,
                                                                     const std::string suffix) {
// *************************
  std::vector<double> MaxValue(TotalNumberOfSamples);
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    for (int toy = 0; toy < Ntoys; ++toy) {
      double max_val = Toys[sample][toy]->GetMaximum();
      MaxValue[sample] = std::max(MaxValue[sample], max_val);
    }
  }

  std::vector<std::unique_ptr<TH2D>> Spectra(TotalNumberOfSamples);
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    // Get MC histogram x-axis binning
    TH1D* refHist = Toys[sample][0].get();

    const int n_bins_x = refHist->GetNbinsX();
    std::vector<double> x_bin_edges(n_bins_x + 1);
    for (int b = 0; b <= n_bins_x; ++b) {
      x_bin_edges[b] = refHist->GetXaxis()->GetBinLowEdge(b + 1); // ROOT bins start at 1
    }
    // Last edge is upper edge of last bin:
    x_bin_edges[n_bins_x] = refHist->GetXaxis()->GetBinUpEdge(n_bins_x);

    constexpr int n_bins_y = 400;
    constexpr double y_min = 0.0;
    const double y_max = MaxValue[sample] * 1.05;

    // Create TH2D with variable binning on x axis
    Spectra[sample] = std::make_unique<TH2D>(
      (SampleNames[sample] + "_" + suffix).c_str(),   // name
      (SampleNames[sample] + "_" + suffix).c_str(),   // title
      n_bins_x, x_bin_edges.data(),                   // x axis bins
      n_bins_y, y_min, y_max                          // y axis bins
    );

    Spectra[sample]->SetDirectory(nullptr);
    Spectra[sample]->Sumw2(true);
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    for (int toy = 0; toy < Ntoys; ++toy) {
      FastViolinFill(Spectra[sample].get(), Toys[sample][toy].get());
    }
  }

  return Spectra;
}

// *************************
std::unique_ptr<TH1D> PredictiveThrower::MakePredictive(const std::vector<std::unique_ptr<TH1D>>& Toys,
                                                        const std::string& Sample_Name,
                                                        const std::string& suffix,
                                                        const bool DebugHistograms) {
// *************************
  constexpr int nXBins = 100;
  int nbinsx = Toys[0]->GetNbinsX();

  auto PredictiveHist = std::make_unique<TH1D>((Sample_Name + "_" + suffix + "_PostPred").c_str(),
                                               (Sample_Name + "_" + suffix + "_PostPred").c_str(),
                                                nbinsx, Toys[0]->GetXaxis()->GetXbins()->GetArray());
  PredictiveHist->SetDirectory(nullptr);

  for (int i = 1; i <= nbinsx; ++i) {
    double x_low  = Toys[0]->GetXaxis()->GetBinLowEdge(i);
    double x_high = Toys[0]->GetXaxis()->GetBinUpEdge(i);
    TString projName = TString::Format("%s %s X: [%.3f, %.3f]", Sample_Name.c_str(), suffix.c_str(), x_low, x_high);
    auto PosteriorHist = std::make_unique<TH1D>(projName, projName, nXBins, 1, -1);
    PosteriorHist->SetDirectory(nullptr);

    for (size_t iToy = 0; iToy < Toys.size(); ++iToy) {
      const double Content = Toys[iToy]->GetBinContent(i);
      PosteriorHist->Fill(Content, ReweightWeight[iToy]);
    }

    if(DebugHistograms) PosteriorHist->Write();

    const double nMean = PosteriorHist->GetMean();
    const double nMeanError = PosteriorHist->GetRMS();

    PredictiveHist->SetBinContent(i, nMean);
    PredictiveHist->SetBinError(i, nMeanError);
  }

  PredictiveHist->Write();
  return PredictiveHist;
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

  for (int sample = 0; sample < TotalNumberOfSamples+1; ++sample) {
    SampleDirectories[sample] = PredictiveDir->mkdir(SampleNames[sample].c_str());
  }

  std::vector<std::unique_ptr<TH2D>> Spectra_mc = ProduceSpectra(MC_Hist_Toy, "mc");
  //std::vector<std::unique_ptr<TH2D>> Spectra_w2 = ProduceSpectra(W2_Hist_Toy, "w2");
  std::vector<std::unique_ptr<TH1D>> PostPred_mc(TotalNumberOfSamples);
  //std::vector<std::unique_ptr<TH1D>> PostPred_w2(TotalNumberOfSamples);
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    SampleDirectories[sample]->cd();
    Spectra_mc[sample]->Write();
    //Spectra_w2[sample]->Write();

    PostPred_mc[sample] = MakePredictive(MC_Hist_Toy[sample], SampleNames[sample], "mc", DebugHistograms);
    //PostPred_w2[sample] = MakePredictive(W2_Hist_Toy[sample], SampleNames[sample], "w2", DebugHistograms);
  }

  PosteriorPredictivepValue(PostPred_mc,
                            //PostPred_w2,
                            SampleDirectories);

  // Close directories
  for (int sample = 0; sample < TotalNumberOfSamples+1; ++sample) {
    SampleDirectories[sample]->Close();
    delete SampleDirectories[sample];
  }

  PredictiveDir->Close();
  delete PredictiveDir;

  outputFile->cd();

  TempClock.Stop();
  MACH3LOG_INFO("{} took {:.2f}s to finish for {} toys", __func__, TempClock.RealTime(), Ntoys);
}

// *************************
double PredictiveThrower::GetLLH(const std::unique_ptr<TH1D>& DatHist,
                                 const std::unique_ptr<TH1D>& MCHist,
                                 const std::unique_ptr<TH1D>& W2Hist,
                                 SampleHandlerBase* SampleHandler) {
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
void PredictiveThrower::PosteriorPredictivepValue(const std::vector<std::unique_ptr<TH1D>>& PostPred_mc,
                                                  //const std::vector<std::unique_ptr<TH1D>>& PostPred_w2,
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
      auto DrawFluctHist = M3::Clone(MC_Hist_Toy[iSample][iToy].get());
      MakeFluctuatedHistogramAlternative(DrawFluctHist.get(), MC_Hist_Toy[iSample][iToy].get(), random.get());

      auto PredFluctHist = M3::Clone(PostPred_mc[iSample].get());
      MakeFluctuatedHistogramAlternative(PredFluctHist.get(), PostPred_mc[iSample].get(), random.get());

      // Okay now we can do our chi2 calculation for our sample
      chi2_dat_vec[iToy][iSample]  = GetLLH(Data_Hist[iSample], MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], samples[SampleObjectMap[iSample]]);
      chi2_mc_vec[iToy][iSample]   = GetLLH(DrawFluctHist, MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], samples[SampleObjectMap[iSample]]);
      chi2_pred_vec[iToy][iSample] = GetLLH(PredFluctHist, MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], samples[SampleObjectMap[iSample]]);

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
                   const std::string Tittle) {
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

    double min_val = std::min(*std::min_element(chi2_y_sample.begin(), chi2_y_sample.end()),
                              *std::min_element(chi2_x_per_sample.begin(), chi2_x_per_sample.end()));
    double max_val = std::max(*std::max_element(chi2_y_sample.begin(), chi2_y_sample.end()),
                              *std::max_element(chi2_x_per_sample.begin(), chi2_x_per_sample.end()));

    auto chi2_hist = std::make_unique<TH2D>((SampleNames[iSample] + Tittle).c_str(),
                                            (SampleNames[iSample] + Tittle).c_str(),
                                            100, min_val, max_val, 100, min_val, max_val);

    chi2_hist->GetXaxis()->SetTitle(Chi2_x_title.c_str());
    chi2_hist->GetYaxis()->SetTitle(Chi2_y_title.c_str());

    for (int iToy = 0; iToy < Ntoys; ++iToy) {
      chi2_hist->Fill(chi2_x_per_sample[iToy], chi2_y_sample[iToy]);
    }

    Get2DBayesianpValue(chi2_hist.get());
    chi2_hist->Write();
  }
}
