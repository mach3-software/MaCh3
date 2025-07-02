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
}

// *************************
// Destructor:
PredictiveThrower::~PredictiveThrower() {
// *************************

}



// *************************
void PredictiveThrower::SetParamters() {
// *************************
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

  /// TODO this need to be changed blarb
  auto DoNotThrowParamGroup = Get<std::vector<std::string>>(fitMan->raw()["Predictive"]["DoNotThrowParamGroup"], __FILE__, __LINE__);
  if(ModelSystematic) ModelSystematic->SetGroupOnlyParameters(DoNotThrowParamGroup);
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

  if(counter > 1) {
    MACH3LOG_ERROR("Found {} ParmaterHandler inheriting from ParameterHandlerGeneric, I can accept at most 1", counter);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (size_t s = 0; s < systematics.size(); ++s) {
    NModelParams += systematics[s]->GetNumParams();
  }
}


// *************************
// Produce MaCh3 toys:
void PredictiveThrower::ProduceToys() {
// *************************
  SetupToyGeneration();

  bool Is_PriorPredictive = Get<bool>(fitMan->raw()["Predictive"]["PriorPredictive"], __FILE__, __LINE__);
  auto PosteriorFileName = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);

  MACH3LOG_INFO("Starting {}", __func__);

  if(Is_PriorPredictive) {
    MACH3LOG_INFO("You've chosen to run Prior Predictive Distribution");
  } else {
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
    if (!MaCh3Sample) {
      MACH3LOG_ERROR(":: Sample {} do not inherit from  SampleHandlerFD this is not implemented::", samples[iPDF]->GetTitle());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Get nominal spectra and event rates
    TH1D* DataHist1D = static_cast<TH1D*>(MaCh3Sample->GetDataHist(1));
    DataHist1D->Write((MaCh3Sample->GetTitle() + "_data").c_str());

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
  int Ntoys = Get<int>(fitMan->raw()["Predictive"]["Ntoy"], __FILE__, __LINE__);

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

    for (size_t iPDF = 0; iPDF < samples.size(); iPDF++) {
      samples[iPDF]->Reweight();
    }

    for (size_t iPDF = 0; iPDF < samples.size(); iPDF++)
    {
      auto* MaCh3Sample = dynamic_cast<SampleHandlerFD*>(samples[iPDF]);

      TH1D* MCHist1D = static_cast<TH1D*>(MaCh3Sample->GetMCHist(1));
      MCHist1D->Write((MaCh3Sample->GetTitle() + "_mc_" + std::to_string(i)).c_str());

      TH1D* W2Hist1D = static_cast<TH1D*>(MaCh3Sample->GetW2Hist(1));
      W2Hist1D->Write((MaCh3Sample->GetTitle() + "_w2_" + std::to_string(i)).c_str());
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
