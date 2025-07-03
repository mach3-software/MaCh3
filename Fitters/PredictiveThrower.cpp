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

  Is_PriorPredictive = false;
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

  /// TODO WARNING BLARB Add fixing single param
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
      MACH3LOG_ERROR(":: Sample {} do not inherit from  SampleHandlerFD this is not implemented::", samples[iPDF]->GetTitle());
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
  SampleNames.resize(TotalNumberOfSamples);
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    MC_Hist_Toy[sample].resize(Ntoys);
    W2_Hist_Toy[sample].resize(Ntoys);
    SampleNames[sample] = samples[sample]->GetTitle();
  }

  Is_PriorPredictive = Get<bool>(fitMan->raw()["Predictive"]["PriorPredictive"], __FILE__, __LINE__);

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


  if(ModelSystematic) {
    auto ThrowParamGroupOnly = GetFromManager<std::vector<std::string>>(fitMan->raw()["Predictive"]["ThrowParamGroupOnly"], {}, __FILE__, __LINE__);
    auto UniqueParamGroup = ModelSystematic->GetUniqueParameterGroups();

    MACH3LOG_INFO("I have following parameter groups: {}", fmt::join(UniqueParamGroup, ", "));
    if(ThrowParamGroupOnly.size() == 0) {
      MACH3LOG_INFO("I will vary all");
    } else {
      // Compute UniqueParamGroup - ThrowParamGroupOnly
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

  for (size_t s = 0; s < systematics.size(); ++s)
  {
    branch_vals[s].resize(systematics[s]->GetNumParams(), M3::_BAD_DOUBLE_);
    branch_name[s].resize(systematics[s]->GetNumParams());
    for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
      branch_name[s][i] = systematics[s]->GetParName(i);
      PosteriorFile->SetBranchAddress(branch_name[s][i].c_str(), &branch_vals[s][i]);
    }
    PosteriorFile->GetEntry(PosteriorFile->GetEntries()-1);

    for (int i = 0; i < systematics[s]->GetNumParams(); ++i) {
      if(branch_vals[s][i] == M3::_BAD_DOUBLE_)
      {
        MACH3LOG_ERROR("Parameter {} is unvitalised with value {}", i, branch_vals[s][i]);
        MACH3LOG_ERROR("Please check more precisely chain you passed {}", PosteriorFileName);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
    }
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
    if(FullLLH){
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
                                                        const std::string& suffix) {
// *************************
  const int nXBins = 500;
  int nbinsx = Toys[0]->GetNbinsX();

  auto hist1d = std::make_unique<TH1D>((Sample_Name + "_" + suffix + "_hist_1d").c_str(),
                                       (Sample_Name + "_" + suffix + "_hist_1d").c_str(),
                                       nbinsx, Toys[0]->GetXaxis()->GetXbins()->GetArray());
  hist1d->SetDirectory(nullptr);

  for (int i = 1; i <= nbinsx; ++i) {
    TString projName = TString::Format("%s_%s_py%d", Sample_Name.c_str(), suffix.c_str(), i);
    auto PosteriorHist = std::make_unique<TH1D>(projName, projName, nXBins, 1, -1);
    PosteriorHist->SetDirectory(nullptr);

    for (size_t iToy = 0; iToy < Toys.size(); ++iToy) {
      const double Content = Toys[iToy]->GetBinContent(i);
      PosteriorHist->Fill(Content);
    }

    PosteriorHist->Write();

    const double nMean = PosteriorHist->GetMean();
    const double nMeanError = PosteriorHist->GetRMS();

    hist1d->SetBinContent(i, nMean);
    hist1d->SetBinError(i, nMeanError);
  }

  hist1d->Write();
  return hist1d;
                                                        }


// *************************
std::unique_ptr<TH1D> PredictiveThrower::MakeSpectra(const std::unique_ptr<TH2D>& Spectra,
                                                     const std::string& Sample_Name,
                                                     const std::string& suffix) {
// *************************
  int nbinsx = Spectra->GetNbinsX();
  std::vector<double> x(nbinsx), y(nbinsx), ex(nbinsx), ey(nbinsx);
  for(int i = 1; i <= nbinsx; i++)
  {
    TString projName = TString::Format("%s_%s_py%d", Sample_Name.c_str(), suffix.c_str(), i);
    TH1D* proj = Spectra->ProjectionY(projName, i, i);
    double x_low  = Spectra->GetXaxis()->GetBinLowEdge(i);
    double x_high = Spectra->GetXaxis()->GetBinUpEdge(i);

    proj->SetTitle(TString::Format("%s %s X: [%.3f, %.3f]",
                                   Sample_Name.c_str(), suffix.c_str(), x_low, x_high));
    if(proj->GetEntries() > 0)
    {
      int first_bin = 99999;
      int last_bin = 0;

      for(int bin_i = 1 ; bin_i < proj->GetNbinsX() - 1 ; bin_i++){
        double val = proj->GetBinContent(bin_i);
        if(val > 0 && bin_i < first_bin){
          first_bin = bin_i;
        }

        if(val > 0){
          last_bin = bin_i;
        }
      }

      // ETA - this is an attempt to automate
      // the bin width in the 1D projections.
      // Often we fit 1D gaussians to these projections
      // and for this the number of bins should really be
      // at least 10 across the range of the projection.

      //Get the range the filled bins cover
      double diff = proj->GetBinCenter(last_bin) - proj->GetBinCenter(first_bin);
      double opt_bin_width = diff/10;
      double current_bin_width = proj->GetBinWidth(last_bin);
      int rebin = 1;

      if(diff > 0) { rebin = static_cast<int>(opt_bin_width / current_bin_width); }
      if(rebin==0){rebin+=1;}
      //If the number of bins in the rebinning isn't a factor
      //of the total number of bins keep reducing it until it is
      while((proj->GetNbinsX() % rebin) != 0){
        rebin--;
      }

      //Just in case the above loop does a bad job
      if(rebin > 0){
        proj->Rebin(rebin);
      }

      //Get the bin centre array
      x[i-1] = Spectra->GetXaxis()->GetBinCenter(i);
      //Set the error on as being half the bin width
      ex[i-1] = Spectra->GetXaxis()->GetBinWidth(i)/2.0;
      //Set the centre of the bin to be
      y[i-1]=proj->GetMean();
      //Set the error on y to be RMS of the bin
      ey[i-1]=proj->GetRMS();
    }
    proj->Write();
    delete proj;
  }

  // Make graph for this sample
  auto errorbars = std::make_unique<TGraphErrors>(nbinsx, x.data(), y.data(), ex.data(), ey.data());
  errorbars->SetName((Sample_Name + "_" + suffix + "_MeanRMS").c_str());
  errorbars->SetTitle("Mean and RMS vs X bin");
  errorbars->GetXaxis()->SetTitle(Spectra->GetXaxis()->GetTitle());
  errorbars->GetYaxis()->SetTitle("Mean +/- RMS");

  // Optionally write graph
  errorbars->Write();

  TAxis* xaxis = Spectra->GetXaxis();
  std::vector<double> x_edges(nbinsx + 1);
  for (int i = 0; i <= nbinsx; ++i) {
    x_edges[i] = xaxis->GetBinLowEdge(i + 1);
  }

  auto hist1d = std::make_unique<TH1D>((Sample_Name + "_" + suffix + "_hist_1d").c_str(), (Sample_Name + "_" + suffix + "_hist_1d").c_str(), nbinsx, x_edges.data());
  hist1d->SetDirectory(nullptr);
  // Fill histogram bin contents and errors (bins start from 1 to nbinsx)
  for (int bin = 1; bin <= nbinsx; ++bin) {
    hist1d->SetBinContent(bin, y[bin - 1]);
    hist1d->SetBinError(bin, ey[bin - 1]);
  }
  hist1d->Write();
  return hist1d;
}

// *************************
// Perform predictive analysis
void PredictiveThrower::RunPredictiveAnalysis() {
// *************************
  MACH3LOG_INFO("Starting {}", __func__);
  TStopwatch TempClock;
  TempClock.Start();

  TDirectory* PredictiveDir = outputFile->mkdir("Predictive");
  std::vector<TDirectory*> SampleDirectories;
  SampleDirectories.resize(TotalNumberOfSamples);

  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
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
    PostPred_mc[sample] = MakeSpectra(Spectra_mc[sample], SampleNames[sample], "mc");
    //PostPred_w2[sample] = MakeSpectra(Spectra_w2[sample], SampleNames[sample], "w2");

    MakePredictive(MC_Hist_Toy[sample], SampleNames[sample], "mc");
  }

  PosteriorPredictivepValue(PostPred_mc,
                            //PostPred_w2,
                            SampleDirectories);

  // Close directories
  for (int sample = 0; sample < TotalNumberOfSamples; ++sample) {
    SampleDirectories[sample]->Close();
    delete SampleDirectories[sample];
  }

  PredictiveDir->Close();
  delete PredictiveDir;

  outputFile->cd();

  TempClock.Stop();
  MACH3LOG_INFO("{} took {:.2f}s to finish for {} toys", __func__, TempClock.RealTime(), Ntoys);
}



// ****************
double PredictiveThrower::GetLLH(const std::unique_ptr<TH1D>& DatHist,
                                 const std::unique_ptr<TH1D>& MCHist,
                                 const std::unique_ptr<TH1D>& W2Hist,
                                 SampleHandlerBase* SampleHandler) {
// ****************
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

  std::vector<std::vector<double>> chi2_dat_vec(Ntoys);
  std::vector<std::vector<double>> chi2_mc_vec(Ntoys);
  std::vector<std::vector<double>> chi2_pred_vec(Ntoys);

  /// TODO add penalty term!!!
  for(int iToy = 0; iToy < Ntoys; iToy++) {
    chi2_dat_vec[iToy].resize(TotalNumberOfSamples, 0.0);
    chi2_mc_vec[iToy].resize(TotalNumberOfSamples, 0.0);
    chi2_pred_vec[iToy].resize(TotalNumberOfSamples, 0.0);

    /// TODO This can be multithreaded
    for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
      auto DrawFluctHist = M3::Clone(MC_Hist_Toy[iSample][iToy].get());
      MakeFluctuatedHistogramAlternative(DrawFluctHist.get(), MC_Hist_Toy[iSample][iToy].get(), random.get());

      auto PredFluctHist = M3::Clone(PostPred_mc[iSample].get());
      MakeFluctuatedHistogramAlternative(PredFluctHist.get(), PostPred_mc[iSample].get(), random.get());

      // Okay now we can do our chi2 calculation for our sample
      chi2_dat_vec[iToy][iSample]  = GetLLH(Data_Hist[iSample], MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], samples[iSample]);
      chi2_mc_vec[iToy][iSample]   = GetLLH(DrawFluctHist, MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], samples[iSample]);
      chi2_pred_vec[iToy][iSample] = GetLLH(PredFluctHist, MC_Hist_Toy[iSample][iToy], W2_Hist_Toy[iSample][iToy], samples[iSample]);

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
  for (int iSample = 0; iSample < TotalNumberOfSamples; ++iSample) {
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
