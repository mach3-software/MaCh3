#include "OscProcessor.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TArrow.h"
_MaCh3_Safe_Include_End_ //}

#pragma GCC diagnostic ignored "-Wfloat-conversion"

// ****************************
OscProcessor::OscProcessor(const std::string &InputFile) : MCMCProcessor(InputFile) {
// ****************************
  //KS: WARNING this only work when you project from Chain, will nor work when you try SetBranchAddress etc. Turn it on only if you know how to use it
  PlotJarlskog = false;

  /// @todo Here where we should add all unitarity triangles, fancy Jarlskog studies and other hacky things that only make sense for oscitations
  Sin2Theta13Index = M3::_BAD_INT_;
  Sin2Theta12Index = M3::_BAD_INT_;
  Sin2Theta23Index = M3::_BAD_INT_;
  DeltaCPIndex = M3::_BAD_INT_;
  DeltaM2_23Index = M3::_BAD_INT_;
}

// ****************************
// The destructor
OscProcessor::~OscProcessor() {
// ****************************

}

// ***************
// Read the Osc cov file and get the input central values and errors
void OscProcessor::LoadAdditionalInfo() {
// ***************
  // KS: Check if OscParams were enabled, in future we will also get
  for(size_t i = 0; i < ParameterGroup.size(); i++) {
    if(ParameterGroup[i] == "Osc"){
      OscEnabled = true;
      break;
    }
  }

  if(OscEnabled)
  {
    for (int i = 0; i < nDraw; ++i)
    {
      //Those keep which parameter type we run currently and relative number
      const int ParamEnum = ParamType[i];
      const int ParamNo = i - ParamTypeStartPos[ParameterEnum(ParamEnum)];
      const std::string CurrentName = ParamNames[ParamEnum][ParamNo].Data();

      /// @todo remove this hardcoding (e.g., use a map or enum-to-name function)
      if (CurrentName == "sin2th_13") {
        Sin2Theta13Index = i;
        Sin2Theta13Name  = CurrentName;
      } else if (CurrentName == "sin2th_12") {
        Sin2Theta12Index = i;
        Sin2Theta12Name  = CurrentName;
      } else if (CurrentName == "sin2th_23") {
        Sin2Theta23Index = i;
        Sin2Theta23Name  = CurrentName;
      } else if (CurrentName == "delta_cp") {
        DeltaCPIndex     = i;
        DeltaCPName      = CurrentName;
      } else if (CurrentName == "delm2_23") {
        DeltaM2_23Index  = i;
        DeltaM2_23Name   = CurrentName;
      }
    }
  } else{
    MACH3LOG_WARN("Didn't find oscillation parameters");
  }

  if(PlotJarlskog && OscEnabled)
  {
    Chain->SetAlias("J_cp", "TMath::Sqrt(sin2th_13)*TMath::Sqrt(1.-sin2th_13)*TMath::Sqrt(1.-sin2th_13)*TMath::Sqrt(sin2th_12)*TMath::Sqrt(1.-sin2th_12)*TMath::Sqrt(sin2th_23)*TMath::Sqrt(1.-sin2th_23)*TMath::Sin(delta_cp)");
    BranchNames.push_back("J_cp");
    ParamType.push_back(kXSecPar);
    nParam[kXSecPar]++;
    nDraw++;

    /// @todo we should actually calculate central value and prior error but leave it for now...
    ParamNom[kXSecPar].push_back( 0. );
    ParamCentral[kXSecPar].push_back( 0. );
    ParamErrors[kXSecPar].push_back( 1. );
    // Push back the name
    ParamNames[kXSecPar].push_back("J_cp");
    ParamFlat[kXSecPar].push_back( false );
  } else if(PlotJarlskog && !OscEnabled) {
    MACH3LOG_ERROR("Trying to enable Jarlskog without oscillations");
    throw MaCh3Exception(__FILE__,__LINE__);
  }
}

// ***************
// Calculate Jarlskog Invariant using oscillation parameters
double OscProcessor::CalcJarlskog(const double s2th13, const double s2th23, const double s2th12, const double dcp) const {
// ***************
  const double s13  = std::sqrt(s2th13);
  const double s23  = std::sqrt(s2th23);
  const double s12  = std::sqrt(s2th12);
  const double sdcp = std::sin(dcp);
  const double c13  = std::sqrt(1.-s2th13);
  const double c12  = std::sqrt(1.-s2th12);
  const double c23  = std::sqrt(1.-s2th23);

  const double j = s13*c13*c13*s12*c12*s23*c23*sdcp;

  return j;
}

// ***************
double OscProcessor::SamplePriorForParam(const int paramIndex, const std::unique_ptr<TRandom3>& randGen, const std::vector<double>& FlatBounds) const {
// ***************
  TString Title = "";
  double Prior = 1.0, PriorError = 1.0;
  bool FlatPrior = false;

  // Get info for this parameter
  GetNthParameter(paramIndex, Prior, PriorError, Title);

  ParameterEnum ParType = ParamType[paramIndex];
  int ParamTemp = paramIndex - ParamTypeStartPos[ParType];
  FlatPrior = ParamFlat[ParType][ParamTemp];

  if (FlatPrior) {
    return randGen->Uniform(FlatBounds[0], FlatBounds[1]);
  } else {
    // Gaussian prior centered at Prior with width PriorError
    return randGen->Gaus(Prior, PriorError);
  }
}

// ***************
// Perform Several Jarlskog Plotting
void OscProcessor::PerformJarlskogAnalysis() {
// ***************
  if(!OscEnabled ||
    Sin2Theta13Index == M3::_BAD_INT_ ||
    Sin2Theta12Index == M3::_BAD_INT_ ||
    Sin2Theta23Index == M3::_BAD_INT_ ||
    DeltaCPIndex == M3::_BAD_INT_||
    DeltaM2_23Index == M3::_BAD_INT_)
  {
    MACH3LOG_WARN("Will not {}, as oscillation parameters are missing", __func__);
    return;
  }
  MACH3LOG_INFO("Starting {}", __func__);

  bool DoReweight = false;

  double s2th13, s2th23, s2th12, dcp, dm2 = M3::_BAD_DOUBLE_;
  double weight = 1.0;
  std::pair<double, double> Sin13_NewPrior;

  // Now read the MCMC file
  TFile *TempFile = new TFile((MCMCFile + ".root").c_str(), "open");

  // Get the settings for the MCMC
  TMacro *Config = TempFile->Get<TMacro>("Reweight_Config");

  if (Config != nullptr) {
    YAML::Node Settings = TMacroToYAML(*Config);
    if(CheckNodeExists(Settings, "Weight", Sin2Theta13Name)) {
      if(Settings["ReweightType"].as<std::string>() == "Gaussian") {
      Sin13_NewPrior = Get<std::pair<double, double>>(Settings["Weight"][Sin2Theta13Name], __FILE__, __LINE__);
      MACH3LOG_INFO("Found Weight in chain, using RC reweighting with new priors {} +- {}", Sin13_NewPrior.first, Sin13_NewPrior.second);
      DoReweight = true;
      } else {
        MACH3LOG_ERROR("ReweightType {} not supported for Jarlskog reweighting", Settings["ReweightType"].as<std::string>());
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }

  TempFile->Close();
  delete TempFile;

  TDirectory *JarlskogDir = OutputFile->mkdir("Jarlskog");
  JarlskogDir->cd();

  unsigned int step = 0;
  Chain->SetBranchStatus("*", false);

  Chain->SetBranchStatus(Sin2Theta13Name.c_str(), true);
  Chain->SetBranchAddress(Sin2Theta13Name.c_str(), &s2th13);

  Chain->SetBranchStatus(Sin2Theta23Name.c_str(), true);
  Chain->SetBranchAddress(Sin2Theta23Name.c_str(), &s2th23);

  Chain->SetBranchStatus(Sin2Theta12Name.c_str(), true);
  Chain->SetBranchAddress(Sin2Theta12Name.c_str(), &s2th12);

  Chain->SetBranchStatus(DeltaCPName.c_str(), true);
  Chain->SetBranchAddress(DeltaCPName.c_str(), &dcp);

  Chain->SetBranchStatus(DeltaM2_23Name.c_str(), true);
  Chain->SetBranchAddress(DeltaM2_23Name.c_str(), &dm2);

  Chain->SetBranchStatus("step", true);
  Chain->SetBranchAddress("step", &step);

  if(DoReweight) {
    Chain->SetBranchStatus("Weight", true);
    Chain->SetBranchAddress("Weight", &weight);
  } else {
    MACH3LOG_WARN("Not applying reweighting weight");
    weight = 1.0;
  }

  // Original histograms
  auto jarl = std::make_unique<TH1D>("jarl", "jarl", 1000, -0.05, 0.05);
  jarl->SetDirectory(nullptr);
  auto jarl_th23 = std::make_unique<TH2D>("jarl_th23", "jarl_th23", 500, -0.05, 0.05, 500, 0.3, 0.7);
  jarl_th23->SetDirectory(nullptr);
  auto jarl_dcp = std::make_unique<TH2D>("jarl_dcp", "jarl_dcp", 500, -0.05, 0.05, 500, -1. * TMath::Pi(), TMath::Pi());
  jarl_dcp->SetDirectory(nullptr);

  jarl->SetTitle("Jarlskog Invariant;J #equiv s_{13}c_{13}^{2}s_{12}c_{12}s_{23}c_{23}sin#delta;Posterior probability");
  jarl_th23->SetTitle("Jarlskog Invariant;J #equiv s_{13}c_{13}^{2}s_{12}c_{12}s_{23}c_{23}sin#delta;Posterior probability");

  // Clones
  auto jarl_IH  = M3::Clone(jarl.get(), "jarl_IH");
  auto jarl_NH  = M3::Clone(jarl.get(), "jarl_NH");

  auto jarl_th23_IH  = M3::Clone(jarl_th23.get(), "jarl_th23_IH");
  auto jarl_th23_NH  = M3::Clone(jarl_th23.get(), "jarl_th23_NH");

  auto jarl_dcp_IH  = M3::Clone(jarl_dcp.get(), "jarl_dcp_IH");
  auto jarl_dcp_NH  = M3::Clone(jarl_dcp.get(), "jarl_dcp_NH");

  auto jarl_flatsindcp          = M3::Clone(jarl.get(), "jarl_flatsindcp");
  auto jarl_IH_flatsindcp       = M3::Clone(jarl.get(), "jarl_IH_flatsindcp");
  auto jarl_NH_flatsindcp       = M3::Clone(jarl.get(), "jarl_NH_flatsindcp");

  auto jarl_th23_flatsindcp     = M3::Clone(jarl_th23.get(), "jarl_th23_flatsindcp");
  auto jarl_th23_IH_flatsindcp  = M3::Clone(jarl_th23.get(), "jarl_th23_IH_flatsindcp");
  auto jarl_th23_NH_flatsindcp  = M3::Clone(jarl_th23.get(), "jarl_th23_NH_flatsindcp");

  auto jarl_prior               = M3::Clone(jarl.get(), "jarl_prior");
  auto jarl_prior_flatsindcp    = M3::Clone(jarl.get(), "jarl_prior_flatsindcp");
  std::unique_ptr<TH1D> jarl_wRC_prior, jarl_wRC_prior_flatsindcp, jarl_wRC_prior_t2kth23;
  // Only use this if chain has reweigh weight [mostly coming from Reactor Constrains]
  if(DoReweight){
    jarl_wRC_prior            = M3::Clone(jarl.get(), "jarl_wRC_prior");
    jarl_wRC_prior_flatsindcp = M3::Clone(jarl.get(), "jarl_wRC_prior_flatsindcp");
    jarl_wRC_prior_t2kth23    = M3::Clone(jarl.get(), "jarl_wRC_prior_flatsindcp");
  }

  // to apply a prior that is flat in sin(dcp) intead of dcp
  auto prior3 = std::make_unique<TF1>("prior3", "TMath::Abs(TMath::Cos(x))");

  // T2K prior is flat (and uncorrelated) in dcp, sin^2(th13), sin^2(th23)
  auto randGen = std::make_unique<TRandom3>(0);
  const Long64_t countwidth = nEntries/5;

  for(int i = 0; i < nEntries; ++i) {
    if (i % countwidth == 0) {
      MaCh3Utils::PrintProgressBar(i, nEntries);
      MaCh3Utils::EstimateDataTransferRate(Chain, i);
    } else {
      Chain->GetEntry(i);
    }

    if(step < BurnInCut) continue; // burn-in cut

    const double j = CalcJarlskog(s2th13, s2th23, s2th12, dcp);
    const double prior_weight = prior3->Eval(dcp);

    jarl->Fill(j, weight);
    jarl_th23->Fill(j, s2th23, weight);
    jarl_dcp->Fill(j, dcp, weight);

    jarl_flatsindcp->Fill(j, prior_weight*weight);
    jarl_th23_flatsindcp->Fill(j, s2th23, prior_weight*weight);

    const double prior_s2th13 = SamplePriorForParam(Sin2Theta13Index, randGen, {0.,1.});
    const double prior_s2th23 = SamplePriorForParam(Sin2Theta23Index, randGen, {0.,1.});
    const double prior_s2th12 = SamplePriorForParam(Sin2Theta12Index, randGen, {0.,1.});
    const double prior_dcp = SamplePriorForParam(DeltaCPIndex, randGen, {-1.*TMath::Pi(),TMath::Pi()});
    // KS: This is hardcoded but we always assume flat in delta CP so probably fine
    const double prior_sindcp = randGen->Uniform(-1., 1.);

    const double prior_s13          = std::sqrt(prior_s2th13);
    const double prior_s23          = std::sqrt(prior_s2th23);
    const double prior_s12          = std::sqrt(prior_s2th12);
    const double prior_sdcp         = std::sin(prior_dcp);
    const double prior_c13          = std::sqrt(1.-prior_s2th13);
    const double prior_c12          = std::sqrt(1.-prior_s2th12);
    const double prior_c23          = std::sqrt(1.-prior_s2th23);
    const double prior_j            = prior_s13*prior_c13*prior_c13*prior_s12*prior_c12*prior_s23*prior_c23*prior_sdcp;
    const double prior_flatsindcp_j = prior_s13*prior_c13*prior_c13*prior_s12*prior_c12*prior_s23*prior_c23*prior_sindcp;

    jarl_prior->Fill(prior_j);
    jarl_prior_flatsindcp->Fill(prior_flatsindcp_j);

    if(DoReweight) {
      const double prior_wRC_s2th13       = randGen->Gaus(Sin13_NewPrior.first, Sin13_NewPrior.second);
      const double prior_wRC_s13          = std::sqrt(prior_wRC_s2th13);
      const double prior_wRC_c13          = std::sqrt(1.-prior_wRC_s2th13);
      const double prior_wRC_j            = prior_wRC_s13*prior_wRC_c13*prior_wRC_c13*prior_s12*prior_c12*prior_s23*prior_c23*prior_sdcp;
      const double prior_wRC_flatsindcp_j = prior_wRC_s13*prior_wRC_c13*prior_wRC_c13*prior_s12*prior_c12*prior_s23*prior_c23*prior_sindcp;
      const double s23                    = std::sqrt(s2th23);
      const double c23                    = std::sqrt(1.-s2th23);

      jarl_wRC_prior->Fill(prior_wRC_j);
      jarl_wRC_prior_flatsindcp->Fill(prior_wRC_flatsindcp_j);
      jarl_wRC_prior_t2kth23->Fill(prior_wRC_s13*prior_wRC_c13*prior_wRC_c13*prior_s12*prior_c12*s23*c23*prior_sdcp);
    }

    if(dm2 > 0.) {
      jarl_NH->Fill(j, weight);
      jarl_th23_NH->Fill(j, s2th23, weight);
      jarl_dcp_NH->Fill(j, dcp, weight);
      jarl_NH_flatsindcp->Fill(j, prior_weight*weight);
      jarl_th23_NH_flatsindcp->Fill(j, s2th23, prior_weight*weight);
    }
    else if(dm2 < 0.) {
      jarl_IH->Fill(j, weight);
      jarl_th23_IH->Fill(j, s2th23, weight);
      jarl_dcp_IH->Fill(j, dcp, weight);
      jarl_IH_flatsindcp->Fill(j, prior_weight*weight);
      jarl_th23_IH_flatsindcp->Fill(j, s2th23, prior_weight*weight);
    }
  }

  jarl->Write("jarlskog_both");
  jarl_NH->Write("jarlskog_NH");
  jarl_IH->Write("jarlskog_IH");
  jarl_th23->Write("jarlskog_th23_both");
  jarl_th23_NH->Write("jarlskog_th23_NH");
  jarl_th23_IH->Write("jarlskog_th23_IH");

  jarl_dcp->Write("jarlskog_dcp_both");
  jarl_dcp_NH->Write("jarlskog_dcp_NH");
  jarl_dcp_IH->Write("jarlskog_dcp_IH");


  jarl_flatsindcp->Write("jarlskog_both_flatsindcp");
  jarl_NH_flatsindcp->Write("jarlskog_NH_flatsindcp");
  jarl_IH_flatsindcp->Write("jarlskog_IH_flatsindcp");
  jarl_th23_flatsindcp->Write("jarlskog_th23_both_flatsindcp");
  jarl_th23_NH_flatsindcp->Write("jarlskog_th23_NH_flatsindcp");
  jarl_th23_IH_flatsindcp->Write("jarlskog_th23_IH_flatsindcp");

  jarl_prior->Write("jarl_prior");
  jarl_prior_flatsindcp->Write("jarl_prior_flatsindcp");
  if(DoReweight) {
    jarl_wRC_prior->Write("jarl_wRC_prior");
    jarl_wRC_prior_flatsindcp->Write("jarl_wRC_prior_flatsindcp");
    jarl_wRC_prior_t2kth23->Write("jarl_wRC_prior_t2kth23");
  }

  MakeJarlskogPlot(jarl, jarl_flatsindcp,
                   jarl_NH, jarl_NH_flatsindcp,
                   jarl_IH, jarl_IH_flatsindcp);

  // Perform Savage Dickey analysis
  if(DoReweight) {
    SavageDickeyPlot(jarl, jarl_wRC_prior, "Jarlskog flat #delta_{CP}", 0);
    SavageDickeyPlot(jarl_flatsindcp, jarl_wRC_prior_flatsindcp, "Jarlskog flat sin#delta_{CP}", 0);
  } else {
    SavageDickeyPlot(jarl, jarl_prior, "Jarlskog flat #delta_{CP}", 0);
    SavageDickeyPlot(jarl_flatsindcp, jarl_prior_flatsindcp, "Jarlskog flat sin#delta_{CP}", 0);
  }

  JarlskogDir->Close();
  delete JarlskogDir;

  Chain->SetBranchStatus("*", true);
  OutputFile->cd();
}


// ***************
void OscProcessor::MakeJarlskogPlot(const std::unique_ptr<TH1D>& jarl,
                                    const std::unique_ptr<TH1D>& jarl_flatsindcp,
                                    const std::unique_ptr<TH1D>& jarl_NH,
                                    const std::unique_ptr<TH1D>& jarl_NH_flatsindcp,
                                    const std::unique_ptr<TH1D>& jarl_IH,
                                    const std::unique_ptr<TH1D>& jarl_IH_flatsindcp) {
// ***************
  MACH3LOG_INFO("Starting {}", __func__);
  int originalErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  // 1-->NH, 0-->both, -1-->IH
  for(int hierarchy = -1; hierarchy <= 1; hierarchy++)
  {
    std::unique_ptr<TH1D> j_hist;
    std::unique_ptr<TH1D> j_hist_sdcp;
    if(hierarchy == 1) {
      j_hist = M3::Clone(jarl_NH.get(), "");
      j_hist_sdcp = M3::Clone(jarl_NH_flatsindcp.get(), "");
      j_hist->SetTitle(";J_{CP} #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta_{CP};Posterior probability");
    } else if(hierarchy == 0) {
      j_hist = M3::Clone(jarl.get(), "");
      j_hist_sdcp = M3::Clone(jarl_flatsindcp.get(), "");
      j_hist->SetTitle(";J_{CP} #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta_{CP};Posterior probability");
    } else if(hierarchy == -1) {
      j_hist = M3::Clone(jarl_IH.get(), "");
      j_hist_sdcp = M3::Clone(jarl_IH_flatsindcp.get(), "");
      j_hist->SetTitle(";J_{CP} #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta_{CP};Posterior probability");
    } else {
      MACH3LOG_ERROR("Invalid hierarchy option. 1 for NH, 0 for both, -1 for IH");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    j_hist->Rebin(7);
    j_hist_sdcp->Rebin(7);

    j_hist->SetLineColor(kAzure-2);
    j_hist_sdcp->SetLineColor(kOrange+1);
    j_hist->SetLineWidth(2);
    j_hist_sdcp->SetLineWidth(2);

    auto StyleAxis = [](TH1* h) {
      auto xAxis = h->GetXaxis();
      auto yAxis = h->GetYaxis();

      xAxis->SetLabelSize(0.04);
      xAxis->SetLabelFont(132);
      xAxis->SetTitleSize(0.04);
      xAxis->SetTitleOffset(0.80);
      xAxis->SetTitleFont(132);
      xAxis->SetNdivisions(505);
      xAxis->SetTickSize(0.04);

      yAxis->SetLabelSize(0.04);
      yAxis->SetLabelFont(132);
      yAxis->SetTitleSize(0.04);
      yAxis->SetTitleOffset(1.2);
      yAxis->SetTitleFont(132);
      yAxis->SetNdivisions(505);
      yAxis->SetTickSize(0.04);
    };

    StyleAxis(j_hist.get());

    j_hist->GetXaxis()->SetRangeUser(-0.04,0.04);
    j_hist->Scale(1./j_hist->Integral());
    j_hist_sdcp->Scale(1./j_hist_sdcp->Integral());

    std::unique_ptr<TH1D> j_hist_copy = M3::Clone(j_hist.get(), "j_hist_copy");
    std::unique_ptr<TH1D> j_hist_1sig = M3::Clone(j_hist.get(), "j_hist_1sig");
    std::unique_ptr<TH1D> j_hist_2sig = M3::Clone(j_hist.get(), "j_hist_2sig");
    std::unique_ptr<TH1D> j_hist_3sig = M3::Clone(j_hist.get(), "j_hist_3sig");

    //upper and lower edges
    double j_bf = j_hist_copy->GetXaxis()->GetBinCenter(j_hist_copy->GetMaximumBin());
    double j_1sig_low = 9999999.;
    double j_1sig_up  = -9999999.;
    double j_2sig_low = 9999999.;;
    double j_2sig_up  = -9999999.;
    double j_3sig_low = 9999999.;;
    double j_3sig_up  = -9999999.;


    std::unique_ptr<TH1D> j_hist_sdcp_copy = M3::Clone(j_hist_sdcp.get(), "j_hist_sdcp_copy");
    std::unique_ptr<TH1D> j_hist_sdcp_1sig = M3::Clone(j_hist_sdcp.get(), "j_hist_sdcp_1sig");
    std::unique_ptr<TH1D> j_hist_sdcp_2sig = M3::Clone(j_hist_sdcp.get(), "j_hist_sdcp_2sig");
    std::unique_ptr<TH1D> j_hist_sdcp_3sig = M3::Clone(j_hist_sdcp.get(), "j_hist_sdcp_3sig");

    //upper and lower edges
    double j_sdcp_1sig_low = 9999999.;
    double j_sdcp_1sig_up  = -9999999.;
    double j_sdcp_2sig_low = 9999999.;;
    double j_sdcp_2sig_up  = -9999999.;
    double j_sdcp_3sig_low = 9999999.;;
    double j_sdcp_3sig_up  = -9999999.;

    double contlevel1 = 0.68;
    double contlevel2 = 0.90;
    double contlevel4 = 0.99;
    double contlevel5 = 0.9973;
    double integral, tsum = 0.;

    integral = j_hist_copy->Integral();

    while((tsum/integral)<contlevel5) {
      double tmax = j_hist_copy->GetMaximum();
      int bin = j_hist_copy->GetMaximumBin();
      double xval = j_hist_copy->GetXaxis()->GetBinCenter(bin);
      double xwidth = j_hist_copy->GetXaxis()->GetBinWidth(bin);
      if((tsum/integral)<contlevel1) {
        j_hist_copy->SetBinContent(bin,-1.0);
        j_hist_1sig->SetBinContent(bin,0.);
        j_hist_2sig->SetBinContent(bin,0.);
        j_hist_3sig->SetBinContent(bin,0.);
        if(xval<j_1sig_low && xval<j_bf) j_1sig_low = xval - xwidth/2.;
        if(xval>j_1sig_up && xval>j_bf) j_1sig_up = xval + xwidth/2.;
      }
      if((tsum/integral)<contlevel2  && (tsum / integral > contlevel1) ) {
        j_hist_copy->SetBinContent(bin,-5.0);
        j_hist_2sig->SetBinContent(bin,0.);
        j_hist_3sig->SetBinContent(bin,0.);
        if(xval<j_2sig_low && xval<j_bf) j_2sig_low = xval - xwidth/2.;
        if(xval>j_2sig_up && xval>j_bf) j_2sig_up = xval + xwidth/2.;
      }
      if((tsum/integral)<contlevel4  && (tsum / integral > contlevel1) ) {
        j_hist_copy->SetBinContent(bin,-9.0);
        j_hist_3sig->SetBinContent(bin,0.);
        if(xval < j_3sig_low && xval <j_bf) j_3sig_low = xval - xwidth/2.;
        if(xval > j_3sig_up && xval > j_bf) j_3sig_up = xval + xwidth/2.;
      }
      tsum+=tmax;
    }

    integral = j_hist_sdcp_copy->Integral();
    tsum = 0.;

    while((tsum/integral)<contlevel5) {
      double tmax = j_hist_sdcp_copy->GetMaximum();
      int bin = j_hist_sdcp_copy->GetMaximumBin();
      double xval = j_hist_sdcp_copy->GetXaxis()->GetBinCenter(bin);
      double xwidth = j_hist_sdcp_copy->GetXaxis()->GetBinWidth(bin);
      if((tsum/integral)<contlevel1) {
        j_hist_sdcp_copy->SetBinContent(bin,-1.0);
        j_hist_sdcp_1sig->SetBinContent(bin,0.);
        j_hist_sdcp_2sig->SetBinContent(bin,0.);
        j_hist_sdcp_3sig->SetBinContent(bin,0.);
        if(xval<j_sdcp_1sig_low && xval<j_bf) j_sdcp_1sig_low = xval - xwidth/2.;
        if(xval>j_sdcp_1sig_up && xval>j_bf) j_sdcp_1sig_up = xval + xwidth/2.;
      }
      if((tsum/integral)<contlevel2  && (tsum / integral > contlevel1) ) {
        j_hist_sdcp_copy->SetBinContent(bin,-5.0);
        j_hist_sdcp_2sig->SetBinContent(bin,0.);
        j_hist_sdcp_3sig->SetBinContent(bin,0.);
        if(xval<j_sdcp_2sig_low && xval<j_bf) j_sdcp_2sig_low = xval - xwidth/2.;
        if(xval>j_sdcp_2sig_up && xval>j_bf) j_sdcp_2sig_up = xval + xwidth/2.;
      }
      if((tsum/integral)<contlevel4  && (tsum / integral > contlevel1) ) {
        j_hist_sdcp_copy->SetBinContent(bin,-9.0);
        j_hist_sdcp_3sig->SetBinContent(bin,0.);
        if(xval<j_sdcp_3sig_low && xval<j_bf) j_sdcp_3sig_low = xval - xwidth/2.;
        if(xval>j_sdcp_3sig_up && xval>j_bf) j_sdcp_3sig_up = xval + xwidth/2.;
      }
      tsum+=tmax;
    }

    j_hist_1sig->SetLineStyle(9);
    j_hist_sdcp_1sig->SetLineStyle(9);
    j_hist_2sig->SetLineStyle(7);
    j_hist_sdcp_2sig->SetLineStyle(7);
    j_hist_3sig->SetLineStyle(2);
    j_hist_sdcp_3sig->SetLineStyle(2);

    auto ldash = std::make_unique<TH1D>("ldash", "solid",  10, -0.04, 0.04);
    auto sdash = std::make_unique<TH1D>("sdash", "dashed", 10, -0.04, 0.04);
    auto fdash = std::make_unique<TH1D>("fdash", "fdashed",10, -0.04, 0.04);
    ldash->SetLineColor(kBlack);
    sdash->SetLineColor(kBlack);
    fdash->SetLineColor(kBlack);
    ldash->SetLineWidth(2);
    sdash->SetLineWidth(2);
    fdash->SetLineWidth(2);
    ldash->SetLineStyle(9);
    sdash->SetLineStyle(7);
    fdash->SetLineStyle(2);

    double vertUp = 0.5 * j_hist->GetMaximum();
    auto jline_1sig_low = std::make_unique<TLine>(j_1sig_low, 0., j_1sig_low, vertUp);
    auto jline_2sig_low = std::make_unique<TLine>(j_2sig_low, 0., j_2sig_low, vertUp);
    auto jline_3sig_low = std::make_unique<TLine>(j_3sig_low, 0., j_3sig_low, vertUp);

    auto jline_1sig_up = std::make_unique<TLine>(j_1sig_up, 0., j_1sig_up,vertUp);
    auto jline_2sig_up = std::make_unique<TLine>(j_2sig_up, 0., j_2sig_up,vertUp);
    auto jline_3sig_up = std::make_unique<TLine>(j_3sig_up, 0., j_3sig_up,vertUp);

    auto jline_sdcp_1sig_low = std::make_unique<TLine>(j_sdcp_1sig_low, 0., j_sdcp_1sig_low, vertUp);
    auto jline_sdcp_2sig_low = std::make_unique<TLine>(j_sdcp_2sig_low, 0., j_sdcp_2sig_low, vertUp);
    auto jline_sdcp_3sig_low = std::make_unique<TLine>(j_sdcp_3sig_low, 0., j_sdcp_3sig_low, vertUp);

    auto jline_sdcp_1sig_up = std::make_unique<TLine>(j_sdcp_1sig_up, 0., j_sdcp_1sig_up, vertUp);
    auto jline_sdcp_2sig_up = std::make_unique<TLine>(j_sdcp_2sig_up, 0., j_sdcp_2sig_up, vertUp);
    auto jline_sdcp_3sig_up = std::make_unique<TLine>(j_sdcp_3sig_up, 0., j_sdcp_3sig_up, vertUp);

    double arrowLength = 0.003;
    double arrowHeight = vertUp;

    auto MakeArrow = [&](double x, Color_t color, Width_t width) -> std::unique_ptr<TArrow> {
      auto arrow = std::make_unique<TArrow>(x, arrowHeight, x - arrowLength, arrowHeight, 0.02, ">");
      arrow->SetLineColor(color);
      arrow->SetLineWidth(width);
      return arrow;
    };

    auto j_arrow_1sig_up = MakeArrow(j_1sig_up, j_hist_1sig->GetLineColor(), j_hist_1sig->GetLineWidth());
    auto j_arrow_2sig_up = MakeArrow(j_2sig_up, j_hist_2sig->GetLineColor(), j_hist_2sig->GetLineWidth());
    auto j_arrow_3sig_up = MakeArrow(j_3sig_up, j_hist_3sig->GetLineColor(), j_hist_3sig->GetLineWidth());

    auto j_sdcp_arrow_1sig_up = MakeArrow(j_sdcp_1sig_up, j_hist_sdcp_1sig->GetLineColor(), j_hist_sdcp_1sig->GetLineWidth());
    auto j_sdcp_arrow_2sig_up = MakeArrow(j_sdcp_2sig_up, j_hist_sdcp_2sig->GetLineColor(), j_hist_sdcp_2sig->GetLineWidth());
    auto j_sdcp_arrow_3sig_up = MakeArrow(j_sdcp_3sig_up, j_hist_sdcp_3sig->GetLineColor(), j_hist_sdcp_3sig->GetLineWidth());

    MACH3LOG_DEBUG("j_1sig_low = {:.4f}, j_2sig_low = {:.4f}, j_3sig_low = {:.4f}", j_1sig_low, j_2sig_low, j_3sig_low);
    MACH3LOG_DEBUG("j_1sig_up = {:.4f}, j_2sig_up = {:.4f}, j_3sig_up = {:.4f}", j_1sig_up, j_2sig_up, j_3sig_up);

    auto CopyLineStyle = [](const TH1D* src, TLine* dst) {
      dst->SetLineColor(src->GetLineColor());
      dst->SetLineStyle(src->GetLineStyle());
      dst->SetLineWidth(src->GetLineWidth());
    };

    CopyLineStyle(j_hist_1sig.get(), jline_1sig_low.get());
    CopyLineStyle(j_hist_1sig.get(), jline_1sig_up.get());
    CopyLineStyle(j_hist_2sig.get(), jline_2sig_low.get());
    CopyLineStyle(j_hist_2sig.get(), jline_2sig_up.get());
    CopyLineStyle(j_hist_3sig.get(), jline_3sig_low.get());
    CopyLineStyle(j_hist_3sig.get(), jline_3sig_up.get());

    CopyLineStyle(j_hist_sdcp_1sig.get(), jline_sdcp_1sig_low.get());
    CopyLineStyle(j_hist_sdcp_1sig.get(), jline_sdcp_1sig_up.get());
    CopyLineStyle(j_hist_sdcp_2sig.get(), jline_sdcp_2sig_low.get());
    CopyLineStyle(j_hist_sdcp_2sig.get(), jline_sdcp_2sig_up.get());
    CopyLineStyle(j_hist_sdcp_3sig.get(), jline_sdcp_3sig_low.get());
    CopyLineStyle(j_hist_sdcp_3sig.get(), jline_sdcp_3sig_up.get());

    auto leg = std::make_unique<TLegend>(0.45, 0.60, 0.75, 0.90);
    leg->SetTextSize(0.05);
    leg->SetFillStyle(0);
    leg->SetNColumns(1);
    leg->SetTextFont(132);
    leg->SetBorderSize(0);

    leg->AddEntry(j_hist.get(), "Prior flat in #delta_{CP}", "l");
    leg->AddEntry(j_hist_sdcp.get(), "Prior flat in sin#delta_{CP}", "l");
    leg->AddEntry(ldash.get(), "68% CI", "l");
    leg->AddEntry(sdash.get(), "90% CI", "l");
    leg->AddEntry(fdash.get(), "99% CI", "l");

    j_hist->GetYaxis()->SetRangeUser(0., j_hist->GetMaximum()*1.15);
    j_hist->Draw("h");
    j_hist_sdcp->Draw("same h");

    jline_sdcp_1sig_up->Draw("same");
    jline_sdcp_2sig_up->Draw("same");
    jline_sdcp_3sig_up->Draw("same");
    jline_1sig_up->Draw("same");
    jline_2sig_up->Draw("same");
    jline_3sig_up->Draw("same");

    j_arrow_1sig_up->Draw();
    j_arrow_2sig_up->Draw();
    j_arrow_3sig_up->Draw();
    j_sdcp_arrow_1sig_up->Draw();
    j_sdcp_arrow_2sig_up->Draw();
    j_sdcp_arrow_3sig_up->Draw();
    leg->Draw("same");

    auto ttext = std::make_unique<TText>();
    ttext->SetNDC(); // Use normalized device coordinates
    ttext->SetTextSize(0.03); // Adjust size as needed
    ttext->SetTextAlign(13); // Align left-top

    if (hierarchy == 1) ttext->DrawText(0.15, 0.85, "Normal Ordering");
    else if (hierarchy == 0) ttext->DrawText(0.15, 0.85, "Both Orderings");
    else if (hierarchy == -1) ttext->DrawText(0.15, 0.85, "Inverted Ordering");

    gPad->RedrawAxis();
    Posterior->Update();
    gPad->Update();

    Posterior->Print(CanvasName);

    if(hierarchy == 1)  Posterior->Write("jarl1D_NH_comp");
    else if(hierarchy == 0)  Posterior->Write("jarl1D_both_comp");
    else if(hierarchy == -1)  Posterior->Write("jarl1D_IH_comp");
  }

  gErrorIgnoreLevel = originalErrorLevel;
}
