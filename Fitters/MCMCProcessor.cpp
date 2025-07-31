#include "MCMCProcessor.h"

_MaCh3_Safe_Include_Start_ //{
#include "TChain.h"
#include "TF1.h"
_MaCh3_Safe_Include_End_ //}

//Only if GPU is enabled
#ifdef MaCh3_CUDA
#include "Fitters/gpuMCMCProcessorUtils.cuh"
#endif

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"

// ****************************
MCMCProcessor::MCMCProcessor(const std::string &InputFile) :
  Chain(nullptr), StepCut(""), MadePostfit(false) {
// ****************************
  MCMCFile = InputFile;

  SetMaCh3LoggerFormat();
  MaCh3Utils::MaCh3Welcome();
  MACH3LOG_INFO("Making post-fit processor for: {}", MCMCFile);

  ParStep = nullptr;
  StepNumber = nullptr;
    
  Posterior = nullptr;
  hviolin = nullptr;
  hviolin_prior = nullptr;

  OutputFile = nullptr;
  
  ParamSums = nullptr;
  BatchedAverages = nullptr;
  SampleValues = nullptr;
  SystValues = nullptr;
  AccProbValues = nullptr;
  AccProbBatchedAverages = nullptr;

  //KS:Hardcoded should be a way to get it via config or something
  plotRelativeToPrior = false;
  printToPDF = false;
  plotBinValue = false;
  PlotFlatPrior = true;
  CacheMCMC = false;
  ApplySmoothing = true;
  FancyPlotNames = true;
  doDiagMCMC = false;

  // KS: ROOT can compile FFT code but it will crash during run time. Turn off FFT dynamically
#ifdef MaCh3_FFT
  useFFTAutoCorrelation = true;
#else
  useFFTAutoCorrelation = false;
#endif
  OutputSuffix = "_Process";
  Post2DPlotThreshold = 1.e-5;

  nDraw = 0;
  nEntries = 0;
  UpperCut = M3::_BAD_INT_;
  nSteps = 0;
  nBatches = 0;
  AutoCorrLag = 0;
  nSysts = 0;
  nSamples = 0;
  
  nBins = 70;
  DrawRange = 1.5;
  
  Posterior1DCut = "";
  //KS:Those keep basic information for ParameterEnum
  ParamNames.resize(kNParameterEnum);
  ParamCentral.resize(kNParameterEnum);
  ParamNom.resize(kNParameterEnum);
  ParamErrors.resize(kNParameterEnum);
  ParamFlat.resize(kNParameterEnum);
  ParamTypeStartPos.resize(kNParameterEnum);
  nParam.resize(kNParameterEnum);
  CovPos.resize(kNParameterEnum);
  CovNamePos.resize(kNParameterEnum);
  CovConfig.resize(kNParameterEnum);

  for(int i = 0; i < kNParameterEnum; i++)
  {
    ParamTypeStartPos[i] = 0;
    nParam[i] = 0;
  }
  //Only if GPU is enabled
  #ifdef MaCh3_CUDA
   ParStep_cpu = nullptr;
   NumeratorSum_cpu = nullptr;
   ParamSums_cpu = nullptr;
   DenomSum_cpu = nullptr;

   ParStep_gpu = nullptr;
   NumeratorSum_gpu = nullptr;
   ParamSums_gpu = nullptr;
   DenomSum_gpu = nullptr;
  #endif
}

// ****************************
// The destructor
MCMCProcessor::~MCMCProcessor() {
// ****************************
  // Close the pdf file
  MACH3LOG_INFO("Closing pdf in MCMCProcessor: {}", CanvasName.Data());
  CanvasName += "]";
  if(printToPDF) Posterior->Print(CanvasName);

  delete Covariance;
  delete Correlation;
  delete Central_Value;
  delete Means;
  delete Errors;
  delete Means_Gauss;
  delete Errors_Gauss;
  delete Means_HPD;
  delete Errors_HPD;
  delete Errors_HPD_Positive;
  delete Errors_HPD_Negative;

  for (int i = 0; i < nDraw; ++i)
  {
    if(hpost[i] != nullptr) delete hpost[i];
  }
  if(CacheMCMC)
  {
    for (int i = 0; i < nDraw; ++i) 
    {
      for (int j = 0; j < nDraw; ++j) 
      {
        delete hpost2D[i][j];
      }
      delete[] ParStep[i];
    }
    delete[] ParStep;
  }
  if(StepNumber != nullptr) delete[] StepNumber;

  if(OutputFile != nullptr) OutputFile->Close();
  if(OutputFile != nullptr) delete OutputFile;
  delete Chain;
}

// ***************
void MCMCProcessor::Initialise(){
// ***************
  // Scan the ROOT file for useful branches
  ScanInput();

  // Setup the output
  SetupOutput();
}

// ***************
void MCMCProcessor::GetPostfit(TVectorD *&Central_PDF, TVectorD *&Errors_PDF, TVectorD *&Central_G, TVectorD *&Errors_G, TVectorD *&Peak_Values) {
// ***************
  // Make the post fit
  MakePostfit();

  // We now have the private members
  Central_PDF = Means;
  Errors_PDF = Errors;
  Central_G = Means_Gauss;
  Errors_G = Errors_Gauss;
  Peak_Values = Means_HPD;
}

// ***************
// Get post-fits for the ParameterEnum type, e.g. xsec params, ND params or flux params etc
void MCMCProcessor::GetPostfit_Ind(TVectorD *&PDF_Central, TVectorD *&PDF_Errors, TVectorD *&Peak_Values, ParameterEnum kParam) {
// ***************
  // Make the post fit
  MakePostfit();

  // Loop over the loaded param types
  const int ParamTypeSize = int(ParamType.size());
  int ParamNumber = 0;
  for (int i = 0; i < ParamTypeSize; ++i) {
    if (ParamType[i] != kParam) continue;
    (*PDF_Central)(ParamNumber) = (*Means)(i);
    (*PDF_Errors)(ParamNumber) = (*Errors)(i);
    (*Peak_Values)(ParamNumber) = (*Means_HPD)(i);
    ++ParamNumber;
  }
}

// ***************
void MCMCProcessor::GetCovariance(TMatrixDSym *&Cov, TMatrixDSym *&Corr) {
// ***************
  if (CacheMCMC) MakeCovariance_MP();
  else MakeCovariance();
  Cov = static_cast<TMatrixDSym*>(Covariance->Clone());
  Corr = static_cast<TMatrixDSym*>(Correlation->Clone());
}

// ***************
void MCMCProcessor::MakeOutputFile() {
// ***************
  //KS: ROOT hates me... but we can create several instances of MCMC Processor, each with own TCanvas ROOT is mad and will delete if there is more than one canvas with the same name, so we add random number to avoid issue
  auto rand = std::make_unique<TRandom3>(0);
  const int uniform = int(rand->Uniform(0, 10000));
  // Open a TCanvas to write the posterior onto
  Posterior = std::make_unique<TCanvas>(("Posterior" + std::to_string(uniform)).c_str(), ("Posterior" + std::to_string(uniform)).c_str(), 0, 0, 1024, 1024);
  //KS: No idea why but ROOT changed treatment of violin in R6. If you have non uniform binning this will results in very hard to see violin plots.
  TCandle::SetScaledViolin(false);

  Posterior->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Posterior->SetTickx();
  Posterior->SetTicky();

  Posterior->SetBottomMargin(0.1);
  Posterior->SetTopMargin(0.05);
  Posterior->SetRightMargin(0.03);
  Posterior->SetLeftMargin(0.15);

  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;
  
  // Output file to write to
  OutputName = MCMCFile + OutputSuffix +".root";

  // Output file
  OutputFile = new TFile(OutputName.c_str(), "recreate");
  OutputFile->cd();
}

// ****************************
//CW: Function to make the post-fit
void MCMCProcessor::MakePostfit() {
// ****************************
  // Check if we've already made post-fit
  if (MadePostfit == true) return;
  MadePostfit = true;

  // Check if the output file is ready
  if (OutputFile == nullptr) MakeOutputFile();
  
  MACH3LOG_INFO("MCMCProcessor is making post-fit plots...");

  int originalErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  // Directory for posteriors
  TDirectory *PostDir = OutputFile->mkdir("Post");
  TDirectory *PostHistDir = OutputFile->mkdir("Post_1d_hists");
  
  // nDraw is number of draws we want to do
  for (int i = 0; i < nDraw; ++i)
  {
    if (i % (nDraw/5) == 0) {
      MaCh3Utils::PrintProgressBar(i, nDraw);
    }
    OutputFile->cd();
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(i, Prior, PriorError, Title);

    // This holds the posterior density
    const double maxi = Chain->GetMaximum(BranchNames[i]);
    const double mini = Chain->GetMinimum(BranchNames[i]);
    // This holds the posterior density
    hpost[i] = new TH1D(BranchNames[i], BranchNames[i], nBins, mini, maxi);
    hpost[i]->SetMinimum(0);
    hpost[i]->GetYaxis()->SetTitle("Steps");
    hpost[i]->GetYaxis()->SetNoExponent(false);

    //KS: Apply additional Cuts, like mass ordering
    std::string CutPosterior1D = "";
    if(Posterior1DCut != "")
    {
      CutPosterior1D = StepCut +" && " + Posterior1DCut;
    }
    else CutPosterior1D = StepCut;

    // Project BranchNames[i] onto hpost, applying stepcut
    Chain->Project(BranchNames[i], BranchNames[i], CutPosterior1D.c_str());

    if(ApplySmoothing) hpost[i]->Smooth();

    (*Central_Value)(i) = Prior;

    double Mean, Err, Err_p, Err_m;
    GetArithmetic(hpost[i], Mean, Err);
    (*Means)(i) = Mean;
    (*Errors)(i) = Err;

    GetGaussian(hpost[i], Gauss.get(), Mean, Err);
    (*Means_Gauss)(i) = Mean;
    (*Errors_Gauss)(i) = Err;

    GetHPD(hpost[i], Mean, Err, Err_p, Err_m);
    (*Means_HPD)(i) = Mean;
    (*Errors_HPD)(i) = Err;
    (*Errors_HPD_Positive)(i) = Err_p;
    (*Errors_HPD_Negative)(i) = Err_m;

    // Write the results from the projection into the TVectors and TMatrices
    (*Covariance)(i,i) = (*Errors)(i)*(*Errors)(i);
    (*Correlation)(i,i) = 1.0;

    //KS: This need to be before SetMaximum(), this way plot is nicer as line end at the maximum
    auto hpd = std::make_unique<TLine>((*Means_HPD)(i), hpost[i]->GetMinimum(), (*Means_HPD)(i), hpost[i]->GetMaximum());
    SetTLineStyle(hpd.get(), kBlack, 2, kSolid);
    
    hpost[i]->SetLineWidth(2);
    hpost[i]->SetLineColor(kBlue-1);
    hpost[i]->SetMaximum(hpost[i]->GetMaximum()*DrawRange);
    hpost[i]->SetTitle(Title);
    hpost[i]->GetXaxis()->SetTitle(hpost[i]->GetTitle());
    
    // Now make the TLine for the Asimov
    auto Asimov = std::make_unique<TLine>(Prior, hpost[i]->GetMinimum(), Prior, hpost[i]->GetMaximum());
    SetTLineStyle(Asimov.get(), kRed-3, 2, kDashed);

    auto leg = std::make_unique<TLegend>(0.12, 0.6, 0.6, 0.97);
    SetLegendStyle(leg.get(), 0.04);
    leg->AddEntry(hpost[i], Form("#splitline{PDF}{#mu = %.2f, #sigma = %.2f}", hpost[i]->GetMean(), hpost[i]->GetRMS()), "l");
    leg->AddEntry(Gauss.get(), Form("#splitline{Gauss}{#mu = %.2f, #sigma = %.2f}", Gauss->GetParameter(1), Gauss->GetParameter(2)), "l");
    leg->AddEntry(hpd.get(), Form("#splitline{HPD}{#mu = %.2f, #sigma = %.2f (+%.2f-%.2f)}", (*Means_HPD)(i), (*Errors_HPD)(i), (*Errors_HPD_Positive)(i), (*Errors_HPD_Negative)(i)), "l");
    leg->AddEntry(Asimov.get(), Form("#splitline{Prior}{x = %.2f , #sigma = %.2f}", Prior, PriorError), "l");

    //CW: Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
    if (hpost[i]->GetMaximum() == hpost[i]->Integral()*DrawRange) 
    {
      MACH3LOG_WARN("Found fixed parameter: {} ({}), moving on", Title, i);
      IamVaried[i] = false;
      //KS:Set mean and error to prior for fixed parameters, it looks much better when fixed parameter has mean on prior rather than on 0 with 0 error.
      (*Means_HPD)(i)  = Prior;
      (*Errors_HPD)(i) = PriorError;
      continue;
    }

    // Store that this parameter is indeed being varied
    IamVaried[i] = true;

    // Write to file
    Posterior->SetName(Title);
    Posterior->SetTitle(Title);

    // Draw onto the TCanvas
    hpost[i]->Draw();
    hpd->Draw("same");
    Asimov->Draw("same");
    leg->Draw("same");  
    
    if(printToPDF) Posterior->Print(CanvasName);
        
    // cd into params directory in root file
    PostDir->cd();
    Posterior->Write();
    
    hpost[i]->SetName(Title);
    hpost[i]->SetTitle(Title);
    PostHistDir->cd();
    hpost[i]->Write();
  } // end the for loop over nDraw

  OutputFile->cd();
  TTree *SettingsBranch = new TTree("Settings", "Settings");
  int CrossSectionParameters = nParam[kXSecPar];
  SettingsBranch->Branch("CrossSectionParameters", &CrossSectionParameters);
  int CrossSectionParametersStartingPos = ParamTypeStartPos[kXSecPar];
  SettingsBranch->Branch("CrossSectionParametersStartingPos", &CrossSectionParametersStartingPos);
  int FluxParameters = GetGroup("Flux");
  SettingsBranch->Branch("FluxParameters", &FluxParameters);
  
  int NDParameters = nParam[kNDPar];
  SettingsBranch->Branch("NDParameters", &NDParameters);
  int NDParametersStartingPos = ParamTypeStartPos[kNDPar];
  SettingsBranch->Branch("NDParametersStartingPos", &NDParametersStartingPos);

  int FDParameters = nParam[kFDDetPar];
  SettingsBranch->Branch("FDParameters", &FDParameters);
  int FDParametersStartingPos = ParamTypeStartPos[kFDDetPar];
  SettingsBranch->Branch("FDParametersStartingPos", &FDParametersStartingPos);
  
  SettingsBranch->Branch("NDSamplesBins", &NDSamplesBins);
  SettingsBranch->Branch("NDSamplesNames", &NDSamplesNames);

  SettingsBranch->Fill();
  SettingsBranch->Write();
  delete SettingsBranch;

  TDirectory *Names = OutputFile->mkdir("Names");
  Names->cd();
  for (std::vector<TString>::iterator it = BranchNames.begin(); it != BranchNames.end(); ++it) {
    TObjString((*it)).Write();
  }
  Names->Close();
  delete Names;

  OutputFile->cd();
  Central_Value->Write("Central_Value");
  Means->Write("PDF_Means");
  Errors->Write("PDF_Error");
  Means_Gauss->Write("Gauss_Means");
  Errors_Gauss->Write("Gauss_Errors");
  Means_HPD->Write("Means_HPD");
  Errors_HPD->Write("Errors_HPD");
  Errors_HPD_Positive->Write("Errors_HPD_Positive");
  Errors_HPD_Negative->Write("Errors_HPD_Negative");

  PostDir->Close();
  delete PostDir;
  PostHistDir->Close();
  delete PostHistDir;

  // restore original warning setting
  gErrorIgnoreLevel = originalErrorLevel;
} // Have now written the postfit projections

// *******************
//CW: Draw the postfit
void MCMCProcessor::DrawPostfit() {
// *******************
  if (OutputFile == nullptr) MakeOutputFile();

  // Make the prefit plot
  std::unique_ptr<TH1D> prefit = MakePrefit();

  prefit->GetXaxis()->SetTitle("");
  // cd into the output file
  OutputFile->cd();
 
  std::string CutPosterior1D = "";
  if(Posterior1DCut != "")
  {
    CutPosterior1D = StepCut +" && " + Posterior1DCut;
  }
  else CutPosterior1D = StepCut;

  // Make a TH1D of the central values and the errors
  TH1D *paramPlot = new TH1D("paramPlot", "paramPlot", nDraw, 0, nDraw);
  paramPlot->SetName("mach3params");
  paramPlot->SetTitle(CutPosterior1D.c_str());
  paramPlot->SetFillStyle(3001);
  paramPlot->SetFillColor(kBlue-1);
  paramPlot->SetMarkerColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerStyle(20);
  paramPlot->SetLineColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerSize(prefit->GetMarkerSize());
  paramPlot->GetXaxis()->SetTitle("");

  // Same but with Gaussian output
  TH1D *paramPlot_Gauss = static_cast<TH1D*>(paramPlot->Clone());
  paramPlot_Gauss->SetMarkerColor(kOrange-5);
  paramPlot_Gauss->SetMarkerStyle(23);
  paramPlot_Gauss->SetLineWidth(2);
  paramPlot_Gauss->SetMarkerSize((prefit->GetMarkerSize())*0.75);
  paramPlot_Gauss->SetFillColor(paramPlot_Gauss->GetMarkerColor());
  paramPlot_Gauss->SetFillStyle(3244);
  paramPlot_Gauss->SetLineColor(paramPlot_Gauss->GetMarkerColor());
  paramPlot_Gauss->GetXaxis()->SetTitle("");

  // Same but with Gaussian output
  TH1D *paramPlot_HPD = static_cast<TH1D*>(paramPlot->Clone());
  paramPlot_HPD->SetMarkerColor(kBlack);
  paramPlot_HPD->SetMarkerStyle(25);
  paramPlot_HPD->SetLineWidth(2);
  paramPlot_HPD->SetMarkerSize((prefit->GetMarkerSize())*0.5);
  paramPlot_HPD->SetFillColor(0);
  paramPlot_HPD->SetFillStyle(0);
  paramPlot_HPD->SetLineColor(paramPlot_HPD->GetMarkerColor());
  paramPlot_HPD->GetXaxis()->SetTitle("");

  // Set labels and data
  for (int i = 0; i < nDraw; ++i)
  {
    //Those keep which parameter type we run currently and realtive number  
    int ParamEnu = ParamType[i];
    int ParamNo = i - ParamTypeStartPos[ParameterEnum(ParamEnu)];

    //KS: Slightly hacky way to get relative to prior or nominal as this is convention we use
    //This only applies for xsec for other systematic types doesn't matter
    double CentralValueTemp = 0;
    double Central, Central_gauss, Central_HPD;
    double Err, Err_Gauss, Err_HPD;
    
    if(plotRelativeToPrior)
    {
      CentralValueTemp = ParamCentral[ParamEnu][ParamNo];
      // Normalise the prior relative the nominal/prior, just the way we get our fit results in MaCh3
      if ( CentralValueTemp != 0) 
      {
        Central = (*Means)(i) / CentralValueTemp;
        Err = (*Errors)(i) / CentralValueTemp;
            
        Central_gauss = (*Means_Gauss)(i) / CentralValueTemp;
        Err_Gauss = (*Errors_Gauss)(i) / CentralValueTemp;
        
        Central_HPD = (*Means_HPD)(i) / CentralValueTemp;
        Err_HPD = (*Errors_HPD)(i) / CentralValueTemp;
      } 
      else {
        Central = 1+(*Means)(i);
        Err = (*Errors)(i);
            
        Central_gauss = 1+(*Means_Gauss)(i);
        Err_Gauss = (*Errors_Gauss)(i);
        
        Central_HPD = 1+(*Means_HPD)(i) ;
        Err_HPD = (*Errors_HPD)(i);
      }
    }
    //KS: Just get value of each parameter without dividing by prior
    else
    {
      Central = (*Means)(i);
      Err = (*Errors)(i);
          
      Central_gauss = (*Means_Gauss)(i);
      Err_Gauss = (*Errors_Gauss)(i);
      
      Central_HPD = (*Means_HPD)(i) ;
      Err_HPD = (*Errors_HPD)(i);
    }
    
    paramPlot->SetBinContent(i+1, Central);
    paramPlot->SetBinError(i+1, Err);

    paramPlot_Gauss->SetBinContent(i+1, Central_gauss);
    paramPlot_Gauss->SetBinError(i+1, Err_Gauss);

    paramPlot_HPD->SetBinContent(i+1, Central_HPD);
    paramPlot_HPD->SetBinError(i+1, Err_HPD);
    
    paramPlot->GetXaxis()->SetBinLabel(i+1, prefit->GetXaxis()->GetBinLabel(i+1));
    paramPlot_Gauss->GetXaxis()->SetBinLabel(i+1, prefit->GetXaxis()->GetBinLabel(i+1));
    paramPlot_HPD->GetXaxis()->SetBinLabel(i+1, prefit->GetXaxis()->GetBinLabel(i+1));
  }
  prefit->GetXaxis()->LabelsOption("v");
  paramPlot->GetXaxis()->LabelsOption("v");\
  paramPlot_Gauss->GetXaxis()->LabelsOption("v");
  paramPlot_HPD->GetXaxis()->LabelsOption("v");

  // Make a TLegend
  auto CompLeg = std::make_unique<TLegend>(0.33, 0.73, 0.76, 0.95);
  CompLeg->AddEntry(prefit.get(), "Prefit", "fp");
  CompLeg->AddEntry(paramPlot, "Postfit PDF", "fp");
  CompLeg->AddEntry(paramPlot_Gauss, "Postfit Gauss", "fp");
  CompLeg->AddEntry(paramPlot_HPD, "Postfit HPD", "lfep");
  CompLeg->SetFillColor(0);
  CompLeg->SetFillStyle(0);
  CompLeg->SetLineWidth(0);
  CompLeg->SetLineStyle(0);
  CompLeg->SetBorderSize(0);

  const std::vector<double> Margins = GetMargins(Posterior);
  Posterior->SetBottomMargin(0.2);

  OutputFile->cd();

  //KS: Plot Xsec and Flux
  /// @todo this need revision
  if (nParam[kXSecPar] > 0)
  {
    const int Start = ParamTypeStartPos[kXSecPar];
    const int nFlux = GetGroup("Flux");
    // Plot the xsec parameters (0 to ~nXsec-nFlux) nXsec == xsec + flux, quite confusing I know
    // Have already looked through the branches earlier
    if(plotRelativeToPrior)  prefit->GetYaxis()->SetTitle("Variation rel. prior"); 
    else prefit->GetYaxis()->SetTitle("Parameter Value");
    prefit->GetYaxis()->SetRangeUser(-2.5, 2.5);

    prefit->GetXaxis()->SetRangeUser(Start, Start + nParam[kXSecPar]-nFlux);
    paramPlot->GetXaxis()->SetRangeUser(Start, Start + nParam[kXSecPar]-nFlux);
    paramPlot_Gauss->GetXaxis()->SetRangeUser(Start, Start+ nParam[kXSecPar]-nFlux);
    paramPlot_HPD->GetXaxis()->SetRangeUser(Start, Start + nParam[kXSecPar]-nFlux);

    // Write the individual ones
    prefit->Write("param_xsec_prefit");
    paramPlot->Write("param_xsec");
    paramPlot_Gauss->Write("param_xsec_gaus");
    paramPlot_HPD->Write("param_xsec_HPD");

    // And the combined
    prefit->Draw("e2");
    paramPlot->Draw("e2, same");
    paramPlot_Gauss->Draw("e2, same");
    paramPlot_HPD->Draw("e1, same");
    CompLeg->Draw("same");
    Posterior->Write("param_xsec_canv");
    if(printToPDF) Posterior->Print(CanvasName);
    Posterior->Clear();
  
    OutputFile->cd();
    // Plot the flux parameters (nXSec to nxsec+nflux) if enabled
    // Have already looked through the branches earlier
    prefit->GetYaxis()->SetRangeUser(0.7, 1.3);
    paramPlot->GetYaxis()->SetRangeUser(0.7, 1.3);
    paramPlot_Gauss->GetYaxis()->SetRangeUser(0.7, 1.3);
    paramPlot_HPD->GetYaxis()->SetRangeUser(0.7, 1.3);
    
    prefit->GetXaxis()->SetRangeUser(Start + nParam[kXSecPar]-nFlux, Start + nParam[kXSecPar]);
    paramPlot->GetXaxis()->SetRangeUser(Start + nParam[kXSecPar]-nFlux, Start + nParam[kXSecPar]);
    paramPlot_Gauss->GetXaxis()->SetRangeUser(Start + nParam[kXSecPar]-nFlux, Start + nParam[kXSecPar]);
    paramPlot_HPD->GetXaxis()->SetRangeUser(Start + nParam[kXSecPar]-nFlux, Start + nParam[kXSecPar]);

    prefit->Write("param_flux_prefit");
    paramPlot->Write("param_flux");
    paramPlot_Gauss->Write("param_flux_gaus");
    paramPlot_HPD->Write("param_flux_HPD");

    prefit->Draw("e2");
    paramPlot->Draw("e2, same");
    paramPlot_Gauss->Draw("e1, same");
    paramPlot_HPD->Draw("e1, same");
    CompLeg->Draw("same");
    Posterior->Write("param_flux_canv");
    if(printToPDF) Posterior->Print(CanvasName);
    Posterior->Clear();
  }
  if(nParam[kNDPar] > 0)
  {
    int Start = ParamTypeStartPos[kNDPar];
    int NDbinCounter = Start;
    //KS: Make prefit postfit for each ND sample, having all of them at the same plot is unreadable
    for(unsigned int i = 0; i < NDSamplesNames.size(); i++ )
    {
      std::string NDname = NDSamplesNames[i];
      NDbinCounter += NDSamplesBins[i];
      OutputFile->cd();
      prefit->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      prefit->GetYaxis()->SetRangeUser(0.6, 1.4);
      prefit->GetXaxis()->SetRangeUser(Start, NDbinCounter);

      paramPlot->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      paramPlot->GetYaxis()->SetRangeUser(0.6, 1.4);
      paramPlot->GetXaxis()->SetRangeUser(Start, NDbinCounter);
      paramPlot->SetTitle(CutPosterior1D.c_str());

      paramPlot_Gauss->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      paramPlot_Gauss->GetYaxis()->SetRangeUser(0.6, 1.4);
      paramPlot_Gauss->GetXaxis()->SetRangeUser(Start, NDbinCounter);
      paramPlot_Gauss->SetTitle(CutPosterior1D.c_str());

      paramPlot_HPD->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      paramPlot_HPD->GetYaxis()->SetRangeUser(0.6, 1.4);
      paramPlot_HPD->GetXaxis()->SetRangeUser(Start, NDbinCounter);
      paramPlot_HPD->SetTitle(CutPosterior1D.c_str());

      prefit->Write(("param_"+NDname+"_prefit").c_str());
      paramPlot->Write(("param_"+NDname).c_str());
      paramPlot_Gauss->Write(("param_"+NDname+"_gaus").c_str());
      paramPlot_HPD->Write(("param_"+NDname+"_HPD").c_str());

      prefit->Draw("e2");
      paramPlot->Draw("e2, same");
      paramPlot_Gauss->Draw("e1, same");
      paramPlot_HPD->Draw("e1, same");
      CompLeg->Draw("same");
      Posterior->Write(("param_"+NDname+"_canv").c_str());
      if(printToPDF) Posterior->Print(CanvasName);
      Posterior->Clear();
      Start += NDSamplesBins[i];
    }
  }
  delete paramPlot;
  delete paramPlot_Gauss;
  delete paramPlot_HPD;

  //KS: Return Margin to default one
  SetMargins(Posterior, Margins);
}

// *********************
// Make fancy Credible Intervals plots
void MCMCProcessor::MakeCredibleIntervals(const std::vector<double>& CredibleIntervals,
                                          const std::vector<Color_t>& CredibleIntervalsColours,
                                          const bool CredibleInSigmas) {
// *********************
  if(hpost[0] == nullptr) MakePostfit();

  MACH3LOG_INFO("Making Credible Intervals ");
  const double LeftMargin = Posterior->GetLeftMargin();
  Posterior->SetLeftMargin(0.15);

  // KS: Sanity check of size and ordering is correct
  CheckCredibleIntervalsOrder(CredibleIntervals, CredibleIntervalsColours);
  const int nCredible = int(CredibleIntervals.size());
  std::vector<std::unique_ptr<TH1D>> hpost_copy(nDraw);
  std::vector<std::vector<std::unique_ptr<TH1D>>> hpost_cl(nDraw);

  //KS: Copy all histograms to be thread safe
  for (int i = 0; i < nDraw; ++i)
  {
    hpost_copy[i] = M3::Clone<TH1D>(hpost[i], Form("hpost_copy_%i", i));
    hpost_cl[i].resize(nCredible);
    for (int j = 0; j < nCredible; ++j)
    {
      hpost_cl[i][j] = M3::Clone<TH1D>(hpost[i], Form("hpost_copy_%i_CL_%f", i, CredibleIntervals[j]));

      //KS: Reset to get rid to TF1 otherwise we run into segfault :(
      hpost_cl[i][j]->Reset("");
      hpost_cl[i][j]->Fill(0.0, 0.0);
    }
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nDraw; ++i)
  {
    /// Scale the histograms so it shows the posterior probability
    hpost_copy[i]->Scale(1. / hpost_copy[i]->Integral());
    for (int j = 0; j < nCredible; ++j)
    {
      // Scale the histograms before gettindg credible intervals
      hpost_cl[i][j]->Scale(1. / hpost_cl[i][j]->Integral());
      GetCredibleIntervalSig(hpost_copy[i], hpost_cl[i][j], CredibleInSigmas, CredibleIntervals[j]);

      hpost_cl[i][j]->SetFillColor(CredibleIntervalsColours[j]);
      hpost_cl[i][j]->SetLineWidth(1);
    }
    hpost_copy[i]->GetYaxis()->SetTitleOffset(1.8);
    hpost_copy[i]->SetLineWidth(1);
    hpost_copy[i]->SetMaximum(hpost_copy[i]->GetMaximum()*1.2);
    hpost_copy[i]->SetLineWidth(2);
    hpost_copy[i]->SetLineColor(kBlack);
    hpost_copy[i]->GetYaxis()->SetTitle("Posterior Probability");
  }

  OutputFile->cd();
  TDirectory *CredibleDir = OutputFile->mkdir("Credible");

  for (int i = 0; i < nDraw; ++i)
  {
    if(!IamVaried[i]) continue;

    // Now make the TLine for the Asimov
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(i, Prior, PriorError, Title);

    auto Asimov = std::make_unique<TLine>(Prior, hpost_copy[i]->GetMinimum(), Prior, hpost_copy[i]->GetMaximum());
    SetTLineStyle(Asimov.get(), kRed-3, 2, kDashed);

    auto legend = std::make_unique<TLegend>(0.20, 0.7, 0.4, 0.92);
    SetLegendStyle(legend.get(), 0.03);
    hpost_copy[i]->Draw("HIST");

    for (int j = 0; j < nCredible; ++j)
        hpost_cl[i][j]->Draw("HIST SAME");
    for (int j = nCredible-1; j >= 0; --j)
    {
      if(CredibleInSigmas)
        legend->AddEntry(hpost_cl[i][j].get(), Form("%.0f#sigma Credible Interval", CredibleIntervals[j]), "f");
      else
        legend->AddEntry(hpost_cl[i][j].get(), Form("%.0f%% Credible Interval", CredibleIntervals[j]*100), "f");
    }
    legend->AddEntry(Asimov.get(), Form("#splitline{Prior}{x = %.2f , #sigma = %.2f}", Prior, PriorError), "l");
    legend->Draw("SAME");
    Asimov->Draw("SAME");

    // Write to file
    Posterior->SetName(hpost[i]->GetName());
    Posterior->SetTitle(hpost[i]->GetTitle());

    if(printToPDF) Posterior->Print(CanvasName);
    // cd into directory in root file
    CredibleDir->cd();
    Posterior->Write();
  }
  CredibleDir->Close();
  delete CredibleDir;

  OutputFile->cd();

  //Set back to normal
  Posterior->SetLeftMargin(LeftMargin);
}

// *********************
// Make fancy violin plots
void MCMCProcessor::MakeViolin() {
// *********************
  //KS: Make sure we have steps
  if(!CacheMCMC) CacheSteps();

  MACH3LOG_INFO("Producing Violin Plot");

  //KS: Find min and max to make histogram in range
  double maxi_y = Chain->GetMaximum(BranchNames[0]);
  double mini_y = Chain->GetMinimum(BranchNames[0]);
  for (int i = 1; i < nDraw; ++i)
  {
    const double max_val = Chain->GetMaximum(BranchNames[i]);
    const double min_val = Chain->GetMinimum(BranchNames[i]);
  
    maxi_y = std::max(maxi_y, max_val);
    mini_y = std::min(mini_y, min_val);
  }

  const int vBins = (maxi_y-mini_y)*25;
  hviolin = std::make_unique<TH2D>("hviolin", "hviolin", nDraw, 0, nDraw, vBins, mini_y, maxi_y);
  hviolin->SetDirectory(nullptr);
  //KS: Prior has larger errors so we increase range and number of bins
  constexpr int PriorFactor = 4;
  hviolin_prior = std::make_unique<TH2D>("hviolin_prior", "hviolin_prior", nDraw, 0, nDraw, PriorFactor*vBins, PriorFactor*mini_y, PriorFactor*maxi_y);
  hviolin_prior->SetDirectory(nullptr);

  auto rand = std::make_unique<TRandom3>(0);
  std::vector<double> PriorVec(nDraw);
  std::vector<double> PriorErrorVec(nDraw);
  std::vector<bool> PriorFlatVec(nDraw);

  for (int x = 0; x < nDraw; ++x)
  {
    TString Title;
    double Prior, PriorError;

    GetNthParameter(x, Prior, PriorError, Title);
    //Set fancy labels
    hviolin->GetXaxis()->SetBinLabel(x+1, Title);
    hviolin_prior->GetXaxis()->SetBinLabel(x+1, Title);
    PriorVec[x] = Prior;
    PriorErrorVec[x] = PriorError;

    ParameterEnum ParType = ParamType[x];
    int ParamTemp = x - ParamTypeStartPos[ParType];
    PriorFlatVec[x] = ParamFlat[ParType][ParamTemp];
  }

  TStopwatch clock;
  clock.Start();

  // nDraw is number of draws we want to do
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int x = 0; x < nDraw; ++x)
  {
    //KS: Consider another treatment for fixed params
    //if (IamVaried[x] == false) continue;
    for (int k = 0; k < nEntries; ++k)
    {
      //KS: Burn in cut
      if(StepNumber[k] < BurnInCut) continue;

      //Only used for Suboptimatlity
      if(StepNumber[k] > UpperCut) continue;

      //KS: We know exactly which x bin we will end up, find y bin. This allow to avoid coslty Fill() and enable multithreading becasue I am master of faster
      const double y = hviolin->GetYaxis()->FindBin(ParStep[x][k]);
      hviolin->SetBinContent(x+1, y,  hviolin->GetBinContent(x+1, y)+1);
    }

    //KS: If we set option to not plot flat prior and param has flat prior then we skip this step
    if(!(!PlotFlatPrior && PriorFlatVec[x]))
    {
      for (int k = 0; k < nEntries; ++k)
      {
        const double Entry = rand->Gaus(PriorVec[x], PriorErrorVec[x]);
        const double y = hviolin_prior->GetYaxis()->FindBin(Entry);
        hviolin_prior->SetBinContent(x+1, y,  hviolin_prior->GetBinContent(x+1, y)+1);
      }
    }
  } // end the for loop over nDraw
  clock.Stop();
  MACH3LOG_INFO("Making Violin plot took {:.2f}s to finish for {} steps", clock.RealTime(), nEntries);

  //KS: Tells how many parameters in one canvas we want
  constexpr int IntervalsSize = 10;
  const int NIntervals = nDraw/IntervalsSize;

  hviolin->GetYaxis()->SetTitle("Parameter Value");
  hviolin->GetXaxis()->SetTitle();
  hviolin->GetXaxis()->LabelsOption("v");
  
  hviolin_prior->GetYaxis()->SetTitle("Parameter Value");
  hviolin_prior->GetXaxis()->SetTitle();
  hviolin_prior->GetXaxis()->LabelsOption("v");

  hviolin_prior->SetLineColor(kRed);
  hviolin_prior->SetMarkerColor(kRed);
  hviolin_prior->SetFillColorAlpha(kRed, 0.35);
  hviolin_prior->SetMarkerStyle(20);
  hviolin_prior->SetMarkerSize(0.5);

  // These control violin width, if you use larger then 1 they will most likely overlay, so be cautious
  hviolin_prior->SetBarWidth(1.0);
  hviolin_prior->SetBarOffset(0);

  hviolin->SetLineColor(kBlue);
  hviolin->SetMarkerColor(kBlue);
  hviolin->SetFillColorAlpha(kBlue, 0.35);
  hviolin->SetMarkerStyle(20);
  hviolin->SetMarkerSize(1.0);
  
  const double BottomMargin = Posterior->GetBottomMargin();
  Posterior->SetBottomMargin(0.2);
    
  OutputFile->cd();
  hviolin->Write("param_violin");
  hviolin_prior->Write("param_violin_prior");
  //KS: This is mostly for example plots, we have full file in the ROOT file so can do much better plot later
  hviolin->GetYaxis()->SetRangeUser(-1, +2);
  hviolin_prior->GetYaxis()->SetRangeUser(-1, +2);
  for (int i = 0; i < NIntervals+1; ++i)
  {
    hviolin->GetXaxis()->SetRangeUser(i*IntervalsSize, i*IntervalsSize+IntervalsSize);
    hviolin_prior->GetXaxis()->SetRangeUser(i*IntervalsSize, i*IntervalsSize+IntervalsSize);
    if(i == NIntervals+1)
    {
      hviolin->GetXaxis()->SetRangeUser(i*IntervalsSize, nDraw); 
      hviolin_prior->GetXaxis()->SetRangeUser(i*IntervalsSize, nDraw);
    }
    //KS: ROOT6 has some additional options, consider updating it. more https://root.cern/doc/master/classTHistPainter.html#HP140b
    hviolin_prior->Draw("violinX(03100300)");
    hviolin->Draw("violinX(03100300) SAME");
    if(printToPDF) Posterior->Print(CanvasName);
  }
  //KS: Return Margin to default one
  Posterior->SetBottomMargin(BottomMargin);
}

// *********************
// Make the post-fit covariance matrix in all dimensions
void MCMCProcessor::MakeCovariance() {
// *********************
  if (OutputFile == nullptr) MakeOutputFile();

  bool HaveMadeDiagonal = false;
  MACH3LOG_INFO("Making post-fit covariances...");
  // Check that the diagonal entries have been filled
  // i.e. MakePostfit() has been called
  for (int i = 0; i < nDraw; ++i) {
    if ((*Covariance)(i,i) == M3::_BAD_DOUBLE_) {
      HaveMadeDiagonal = false;
      MACH3LOG_INFO("Have not run diagonal elements in covariance, will do so now by calling MakePostfit()");
      break;
    } else {
      HaveMadeDiagonal = true;
    }
  }

  if (HaveMadeDiagonal == false) {
    MakePostfit();
  }
  gStyle->SetPalette(55);
  // Now we are sure we have the diagonal elements, let's make the off-diagonals
  for (int i = 0; i < nDraw; ++i)
  {
    if (i % (nDraw/5) == 0)
      MaCh3Utils::PrintProgressBar(i, nDraw);

    TString Title_i = "";
    double Prior_i, PriorError;

    GetNthParameter(i, Prior_i, PriorError, Title_i);
    
    const double min_i = Chain->GetMinimum(BranchNames[i]);
    const double max_i = Chain->GetMaximum(BranchNames[i]);

    // Loop over the other parameters to get the correlations
    for (int j = 0; j <= i; ++j) {
      // Skip the diagonal elements which we've already done above
      if (j == i) continue;

      // If this parameter isn't varied
      if (IamVaried[j] == false) {
        (*Covariance)(i,j) = 0.0;
        (*Covariance)(j,i) = (*Covariance)(i,j);
        (*Correlation)(i,j) = 0.0;
        (*Correlation)(j,i) = (*Correlation)(i,j);
        continue;
      }

      TString Title_j = "";
      double Prior_j, PriorError_j;
      GetNthParameter(j, Prior_j, PriorError_j, Title_j);

      OutputFile->cd();

      // The draw which we want to perform
      TString DrawMe = BranchNames[j]+":"+BranchNames[i];

      const double max_j = Chain->GetMaximum(BranchNames[j]);
      const double min_j = Chain->GetMinimum(BranchNames[j]);

      // TH2F to hold the Correlation 
      std::unique_ptr<TH2D> hpost_2D = std::make_unique<TH2D>(DrawMe, DrawMe, nBins, min_i, max_i, nBins, min_j, max_j);
      hpost_2D->SetMinimum(0);
      hpost_2D->GetXaxis()->SetTitle(Title_i);
      hpost_2D->GetYaxis()->SetTitle(Title_j);
      hpost_2D->GetZaxis()->SetTitle("Steps");

      // The draw command we want, i.e. draw param j vs param i
      Chain->Project(DrawMe, DrawMe, StepCut.c_str());
      
      if(ApplySmoothing) hpost_2D->Smooth();
      // Get the Covariance for these two parameters
      (*Covariance)(i,j) = hpost_2D->GetCovariance();
      (*Covariance)(j,i) = (*Covariance)(i,j);

      (*Correlation)(i,j) = hpost_2D->GetCorrelationFactor();
      (*Correlation)(j,i) = (*Correlation)(i,j);

      if(printToPDF)
      {
        //KS: Skip Flux Params
        if(ParamType[i] == kXSecPar && ParamType[j] == kXSecPar)
        {
          if(std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold)
          {
            Posterior->cd();
            hpost_2D->Draw("colz");
            Posterior->SetName(hpost_2D->GetName());
            Posterior->SetTitle(hpost_2D->GetTitle());
            Posterior->Print(CanvasName);
          }
        }
      }
      // Write it to root file
      //OutputFile->cd();
      //if( std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold ) hpost_2D->Write();
    } // End j loop
  } // End i loop
  OutputFile->cd();
  Covariance->Write("Covariance");
  Correlation->Write("Correlation");
}

// ***************
//KS: Cache all steps to allow multithreading, hit RAM quite a bit
void MCMCProcessor::CacheSteps() {
// ***************
  if(CacheMCMC == true) return;
  
  CacheMCMC = true;
  
  if(ParStep != nullptr)
  {
    MACH3LOG_ERROR("It look like ParStep was already filled ");
    MACH3LOG_ERROR("Even though it is used for MakeCovariance_MP and for DiagMCMC ");
    MACH3LOG_ERROR("it has different structure in both for cache hits, sorry ");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  MACH3LOG_INFO("Caching input tree...");
  MACH3LOG_INFO("Allocating {:.2f} MB", double(sizeof(double)*nDraw*nEntries)/1.E6);
  TStopwatch clock;
  clock.Start();
  
  ParStep = new double*[nDraw];
  StepNumber = new int[nEntries];
  
  hpost2D.resize(nDraw);
  for (int i = 0; i < nDraw; ++i) 
  {
    ParStep[i] = new double[nEntries];
    hpost2D[i].resize(nDraw);
    for (int j = 0; j < nEntries; ++j)
    {
      ParStep[i][j] = -999.99;
      //KS: Set this only once
      if(i == 0) StepNumber[j] = -999.99;
    }
  }

  // Set all the branches to off
  Chain->SetBranchStatus("*", false);
  int stepBranch = 0;
  double* ParValBranch = new double[nEntries]();
  // Turn on the branches which we want for parameters
  for (int i = 0; i < nDraw; ++i)
  {
    Chain->SetBranchStatus(BranchNames[i].Data(), true);
    Chain->SetBranchAddress(BranchNames[i].Data(), &ParValBranch[i]);
  }
  Chain->SetBranchStatus("step", true);
  Chain->SetBranchAddress("step", &stepBranch);
  const Long64_t countwidth = nEntries/10;

  // Loop over the entries
  //KS: This is really a bottleneck right now, thus revisit with ROOT6 https://pep-root6.github.io/docs/analysis/parallell/root.html
  for (Long64_t j = 0; j < nEntries; ++j) 
  {
    if (j % countwidth == 0) {
        MaCh3Utils::PrintProgressBar(j, nEntries);
        MaCh3Utils::EstimateDataTransferRate(Chain, j);
    } else {
      Chain->GetEntry(j);
    }
    StepNumber[j] = stepBranch;
    // Set the branch addresses for params
    for (int i = 0; i < nDraw; ++i) 
    {
      ParStep[i][j] = ParValBranch[i];
    }
  }
  delete[] ParValBranch;

  // Set all the branches to on
  Chain->SetBranchStatus("*", true);
  
  // KS: Set temporary branch address to allow min/max, otherwise ROOT can segfaults
  double tempVal = 0.0;
  std::vector<double> Min_Chain(nDraw);
  std::vector<double> Max_Chain(nDraw);
  for (int i = 0; i < nDraw; ++i)
  {
    Chain->SetBranchAddress(BranchNames[i].Data(), &tempVal);
    Min_Chain[i] = Chain->GetMinimum(BranchNames[i]);
    Max_Chain[i] = Chain->GetMaximum(BranchNames[i]);
  }

  // Cache max and min in chain for covariance matrix
  for (int i = 0; i < nDraw; ++i)
  {
    TString Title_i = "";
    double Prior_i, PriorError_i;
    GetNthParameter(i, Prior_i, PriorError_i, Title_i);

    for (int j = 0; j <= i; ++j)
    {
      // TH2D to hold the Correlation
      hpost2D[i][j] = new TH2D(Form("hpost2D_%i_%i",i,j), Form("hpost2D_%i_%i",i,j), nBins, Min_Chain[i], Max_Chain[i], nBins, Min_Chain[j], Max_Chain[j]);
      TString Title_j = "";
      double Prior_j, PriorError_j;
      GetNthParameter(j, Prior_j, PriorError_j, Title_j);

      hpost2D[i][j]->SetMinimum(0);
      hpost2D[i][j]->GetXaxis()->SetTitle(Title_i);
      hpost2D[i][j]->GetYaxis()->SetTitle(Title_j);
      hpost2D[i][j]->GetZaxis()->SetTitle("Steps");
    }
  }
  clock.Stop();
  MACH3LOG_INFO("Caching steps took {:.2f}s to finish for {} steps", clock.RealTime(), nEntries );
}

// *********************
// Make the post-fit covariance matrix in all dimensions
void MCMCProcessor::MakeCovariance_MP(bool Mute) {
// *********************
  if (OutputFile == nullptr) MakeOutputFile();
    
  if(!CacheMCMC) CacheSteps();
  
  bool HaveMadeDiagonal = false;    
  // Check that the diagonal entries have been filled
  // i.e. MakePostfit() has been called
  for (int i = 0; i < nDraw; ++i) {
    if ((*Covariance)(i,i) == M3::_BAD_DOUBLE_) {
      HaveMadeDiagonal = false;
      MACH3LOG_WARN("Have not run diagonal elements in covariance, will do so now by calling MakePostfit()");
      break;
    } else {
      HaveMadeDiagonal = true;
    }
  }
    
  if (HaveMadeDiagonal == false) MakePostfit();
  if(!Mute) MACH3LOG_INFO("Calculating covariance matrix");
  TStopwatch clock;
  if(!Mute) clock.Start();

  gStyle->SetPalette(55);
  // Now we are sure we have the diagonal elements, let's make the off-diagonals
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nDraw; ++i)
  {    
    for (int j = 0; j <= i; ++j)
    {
      // Skip the diagonal elements which we've already done above
      if (j == i) continue;
      
      // If this parameter isn't varied
      if (IamVaried[j] == false) {
        (*Covariance)(i,j) = 0.0;
        (*Covariance)(j,i) = (*Covariance)(i,j);
        (*Correlation)(i,j) = 0.0;
        (*Correlation)(j,i) = (*Correlation)(i,j);
        continue;
      }
      hpost2D[i][j]->SetMinimum(0);
          
      for (int k = 0; k < nEntries; ++k)
      {
        //KS: Burn in cut
        if(StepNumber[k] < BurnInCut) continue;

        //KS: Fill histogram with cached steps
        hpost2D[i][j]->Fill(ParStep[i][k], ParStep[j][k]);
      }
      if(ApplySmoothing) hpost2D[i][j]->Smooth();
      
      // Get the Covariance for these two parameters
      (*Covariance)(i,j) = hpost2D[i][j]->GetCovariance();
      (*Covariance)(j,i) = (*Covariance)(i,j);

      //KS: Since we already have covariance consider calculating correlation using it, right now we effectively calculate covariance twice
      //https://root.cern.ch/doc/master/TH2_8cxx_source.html#l01099
      (*Correlation)(i,j) = hpost2D[i][j]->GetCorrelationFactor();
      (*Correlation)(j,i) = (*Correlation)(i,j);
    }// End j loop
  }// End i loop

  if(!Mute) {
    clock.Stop();
    MACH3LOG_INFO("Making Covariance took {:.2f}s to finish for {} steps", clock.RealTime(), nEntries);
  }
  OutputFile->cd();
  if(printToPDF)
  {
    Posterior->cd();
    for (int i = 0; i < nDraw; ++i)
    {    
      for (int j = 0; j <= i; ++j)
      {
        // Skip the diagonal elements which we've already done above
        if (j == i) continue;
        if (IamVaried[j] == false) continue;

        if(ParamType[i] == kXSecPar && ParamType[j] == kXSecPar)
        {
          if(std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold)
          {
            hpost2D[i][j]->Draw("colz");
            Posterior->SetName(hpost2D[i][j]->GetName());
            Posterior->SetTitle(hpost2D[i][j]->GetTitle());
            Posterior->Print(CanvasName);
          }
        }
        //if( std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold) hpost2D[i][j]->Write();
      }// End j loop
    }// End i loop
  } //end if pdf
  if(!Mute) {
    Covariance->Write("Covariance");
    Correlation->Write("Correlation");
  }
}

// *********************
// Based on @cite roberts2009adaptive
// all credits for finding and studying it goes to Henry
void MCMCProcessor::MakeSubOptimality(const int NIntervals) {
// *********************

  //Save burn in cut, at the end of the loop we will return to default values
  const int DefaultUpperCut = UpperCut;
  const int DefaultBurnInCut = BurnInCut;
  bool defaultPrintToPDF = printToPDF;
  BurnInCut = 0;
  UpperCut = 0;
  printToPDF = false;

  //Set via config in future
  int MaxStep = nSteps;
  int MinStep = 0;
  const int IntervalsSize = nSteps/NIntervals;

  MACH3LOG_INFO("Making Suboptimality");
  TStopwatch clock;
  clock.Start();

  std::unique_ptr<TH1D> SubOptimality = std::make_unique<TH1D>("Suboptimality", "Suboptimality", NIntervals, MinStep, MaxStep);
  SubOptimality->GetXaxis()->SetTitle("Step");
  SubOptimality->GetYaxis()->SetTitle("Suboptimality");
  SubOptimality->SetLineWidth(2);
  SubOptimality->SetLineColor(kBlue);

  for(int i = 0; i < NIntervals; ++i)
  {
    //Reset our cov matrix
    ResetHistograms();

    //Set threshold for calculating new matrix
    UpperCut = i*IntervalsSize;
    //Calculate cov matrix
    MakeCovariance_MP(true);

    //Calculate eigen values
    TMatrixDSymEigen eigen(*Covariance);
    TVectorD eigen_values;
    eigen_values.ResizeTo(eigen.GetEigenValues());
    eigen_values = eigen.GetEigenValues();

    //KS: Converting from ROOT to vector as to make using other libraires (Eigen) easier in future
    std::vector<double> EigenValues(eigen_values.GetNrows());
    for(unsigned int j = 0; j < EigenValues.size(); j++)
    {
      EigenValues[j] = eigen_values(j);
    }
    const double SubOptimalityValue = GetSubOptimality(EigenValues, nDraw);
    SubOptimality->SetBinContent(i+1, SubOptimalityValue);
  }
  clock.Stop();
  MACH3LOG_INFO("Making Suboptimality took {:.2f}s to finish for {} steps", clock.RealTime(), nEntries);

  UpperCut = DefaultUpperCut;
  BurnInCut = DefaultBurnInCut;
  printToPDF = defaultPrintToPDF;

  SubOptimality->Draw("l");
  Posterior->SetName(SubOptimality->GetName());
  Posterior->SetTitle(SubOptimality->GetTitle());

  if(printToPDF) Posterior->Print(CanvasName);
  // Write it to root file
  OutputFile->cd();
  Posterior->Write();
}

// *********************
// Make the covariance plots
void MCMCProcessor::DrawCovariance() {
// *********************
  const double RightMargin  = Posterior->GetRightMargin();
  Posterior->SetRightMargin(0.15);

  // The Covariance matrix from the fit
  std::unique_ptr<TH2D> hCov = std::make_unique<TH2D>("hCov", "hCov", nDraw, 0, nDraw, nDraw, 0, nDraw);
  hCov->GetZaxis()->SetTitle("Covariance");
  // The Covariance matrix square root, with correct sign
  std::unique_ptr<TH2D> hCovSq = std::make_unique<TH2D>("hCovSq", "hCovSq", nDraw, 0, nDraw, nDraw, 0, nDraw);
  hCovSq->GetZaxis()->SetTitle("Covariance");
  // The Correlation
  std::unique_ptr<TH2D> hCorr = std::make_unique<TH2D>("hCorr", "hCorr", nDraw, 0, nDraw, nDraw, 0, nDraw);
  hCorr->GetZaxis()->SetTitle("Correlation");
  hCorr->SetMinimum(-1);
  hCorr->SetMaximum(1);
  hCov->GetXaxis()->SetLabelSize(0.015);
  hCov->GetYaxis()->SetLabelSize(0.015);
  hCovSq->GetXaxis()->SetLabelSize(0.015);
  hCovSq->GetYaxis()->SetLabelSize(0.015);
  hCorr->GetXaxis()->SetLabelSize(0.015);
  hCorr->GetYaxis()->SetLabelSize(0.015);

  // Loop over the Covariance matrix entries
  for (int i = 0; i < nDraw; ++i)
  {
    TString titlex = "";
    double nom, err;
    GetNthParameter(i, nom, err, titlex);
    
    hCov->GetXaxis()->SetBinLabel(i+1, titlex);
    hCovSq->GetXaxis()->SetBinLabel(i+1, titlex);
    hCorr->GetXaxis()->SetBinLabel(i+1, titlex);

    for (int j = 0; j < nDraw; ++j)
    {
      // The value of the Covariance
      const double cov = (*Covariance)(i,j);
      const double corr = (*Correlation)(i,j);

      hCov->SetBinContent(i+1, j+1, cov);
      hCovSq->SetBinContent(i+1, j+1, ((cov > 0) - (cov < 0))*std::sqrt(std::fabs(cov)));
      hCorr->SetBinContent(i+1, j+1, corr);

      TString titley = "";
      double nom_j, err_j;
      GetNthParameter(j, nom_j, err_j, titley);

      hCov->GetYaxis()->SetBinLabel(j+1, titley);
      hCovSq->GetYaxis()->SetBinLabel(j+1, titley);
      hCorr->GetYaxis()->SetBinLabel(j+1, titley);
    }
  }

  // Take away the stat box
  gStyle->SetOptStat(0);
  if(plotBinValue)gStyle->SetPaintTextFormat("4.1f"); //Precision of value in matrix element
  // Make pretty Correlation colors (red to blue)
  const int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  // cd into the correlation directory
  OutputFile->cd();

  Posterior->cd();
  Posterior->Clear();
  if(plotBinValue) hCov->Draw("colz text");
  else hCov->Draw("colz");
  if(printToPDF) Posterior->Print(CanvasName);

  Posterior->cd();
  Posterior->Clear();
  if(plotBinValue) hCorr->Draw("colz text");
  else hCorr->Draw("colz");
  if(printToPDF) Posterior->Print(CanvasName);

  hCov->Write("Covariance_plot");
  hCovSq->Write("Covariance_sq_plot");
  hCorr->Write("Correlation_plot");
  
  //Back to normal
  Posterior->SetRightMargin(RightMargin);
  DrawCorrelations1D();
}

// *********************
//KS: Make the 1D projections of Correlations inspired by Henry's slides (page 28) https://www.t2k.org/asg/oagroup/meeting/2023/2023-07-10-oa-pre-meeting/MaCh3FDUpdate
void MCMCProcessor::DrawCorrelations1D() {
// *********************
  //KS: Store it as we go back to them at the end
  const std::vector<double> Margins = GetMargins(Posterior);
  const int OptTitle = gStyle->GetOptTitle();

  Posterior->SetTopMargin(0.1);
  Posterior->SetBottomMargin(0.2);
  gStyle->SetOptTitle(1);

  constexpr int Nhists = 3;
  //KS: Highest value is just meant bo be sliglhy higher than 1 to catch >,
  constexpr double Thresholds[Nhists+1] = {0, 0.25, 0.5, 1.0001};
  constexpr Color_t CorrColours[Nhists] = {kRed-10, kRed-6,  kRed};

  //KS: This store necessary entries for stripped covariance which store only "meaningful correlations
  std::vector<std::vector<double>> CorrOfInterest;
  CorrOfInterest.resize(nDraw);
  std::vector<std::vector<std::string>> NameCorrOfInterest;
  NameCorrOfInterest.resize(nDraw);

  std::vector<std::vector<std::unique_ptr<TH1D>>> Corr1DHist(nDraw);
  //KS: Initialising ROOT objects is never safe in MP loop
  for(int i = 0; i < nDraw; ++i)
  {
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(i, Prior, PriorError, Title);

    Corr1DHist[i].resize(Nhists);
    for(int j = 0; j < Nhists; ++j)
    {
      Corr1DHist[i][j] = std::make_unique<TH1D>(Form("Corr1DHist_%i_%i", i, j), Form("Corr1DHist_%i_%i", i, j), nDraw, 0, nDraw);
      Corr1DHist[i][j]->SetTitle(Form("%s",Title.Data()));
      Corr1DHist[i][j]->GetYaxis()->SetTitle("Correlation");
      Corr1DHist[i][j]->SetFillColor(CorrColours[j]);
      Corr1DHist[i][j]->SetLineColor(kBlack);

      for (int k = 0; k < nDraw; ++k)
      {
        TString Title_y = "";
        double Prior_y = 1.0;
        double PriorError_y = 1.0;
        GetNthParameter(k, Prior_y, PriorError_y, Title_y);
        Corr1DHist[i][j]->GetXaxis()->SetBinLabel(k+1, Title_y.Data());
      }
    }
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int i = 0; i < nDraw; ++i)
  {
    for(int j = 0; j < nDraw; ++j)
    {
      for(int k = 0; k < Nhists; ++k)
      {
        const double TempEntry = std::fabs((*Correlation)(i,j));
        if(Thresholds[k+1] > TempEntry && TempEntry >= Thresholds[k])
        {
          Corr1DHist[i][k]->SetBinContent(j+1, (*Correlation)(i,j));
        }
      }
      if(std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold && i != j)
      {
        CorrOfInterest[i].push_back((*Correlation)(i,j));
        NameCorrOfInterest[i].push_back(Corr1DHist[i][0]->GetXaxis()->GetBinLabel(j+1));
      }
    }
  }

  TDirectory *CorrDir = OutputFile->mkdir("Corr1D");
  CorrDir->cd();

  for(int i = 0; i < nDraw; i++)
  {
    if (IamVaried[i] == false) continue;

    Corr1DHist[i][0]->SetMaximum(+1.);
    Corr1DHist[i][0]->SetMinimum(-1.);
    Corr1DHist[i][0]->Draw();
    for(int k = 1; k < Nhists; k++) {
      Corr1DHist[i][k]->Draw("SAME");
    }

    auto leg = std::make_unique<TLegend>(0.3, 0.75, 0.6, 0.90);
    SetLegendStyle(leg.get(), 0.02);
    for(int k = 0; k < Nhists; k++) {
      leg->AddEntry(Corr1DHist[i][k].get(), Form("%.2f > |Corr| >= %.2f", Thresholds[k+1], Thresholds[k]), "f");
    }
    leg->Draw("SAME");

    Posterior->Write(Corr1DHist[i][0]->GetTitle());
    if(printToPDF) Posterior->Print(CanvasName);
  }

  //KS: Plot only meaningful correlations
  for(int i = 0; i < nDraw; i++)
  {
    const int size = int(CorrOfInterest[i].size());

    if(size == 0) continue;
    auto Corr1DHist_Reduced = std::make_unique<TH1D>("Corr1DHist_Reduced", "Corr1DHist_Reduced", size, 0, size);
    Corr1DHist_Reduced->SetTitle(Corr1DHist[i][0]->GetTitle());
    Corr1DHist_Reduced->GetYaxis()->SetTitle("Correlation");
    Corr1DHist_Reduced->SetFillColor(kBlue);
    Corr1DHist_Reduced->SetLineColor(kBlue);

    for (int j = 0; j < size; ++j)
    {
      Corr1DHist_Reduced->GetXaxis()->SetBinLabel(j+1, NameCorrOfInterest[i][j].c_str());
      Corr1DHist_Reduced->SetBinContent(j+1, CorrOfInterest[i][j]);
    }
    Corr1DHist_Reduced->GetXaxis()->LabelsOption("v");

    Corr1DHist_Reduced->SetMaximum(+1.);
    Corr1DHist_Reduced->SetMinimum(-1.);
    Corr1DHist_Reduced->Draw();

    Posterior->Write(Form("%s_Red", Corr1DHist_Reduced->GetTitle()));
    if(printToPDF) Posterior->Print(CanvasName);
  }

  CorrDir->Close();
  delete CorrDir;
  OutputFile->cd();

  SetMargins(Posterior, Margins);
  gStyle->SetOptTitle(OptTitle);
}

// *********************
// Make fancy Credible Intervals plots
void MCMCProcessor::MakeCredibleRegions(const std::vector<double>& CredibleRegions,
                                        const std::vector<Style_t>& CredibleRegionStyle,
                                        const std::vector<Color_t>& CredibleRegionColor,
                                        const bool CredibleInSigmas, 
					const bool Draw2DPosterior,
					const bool DrawBestFit) {
// *********************
  if(hpost2D.size() == 0) MakeCovariance_MP();
  MACH3LOG_INFO("Making Credible Regions");

  CheckCredibleRegionsOrder(CredibleRegions, CredibleRegionStyle, CredibleRegionColor);
  const int nCredible = int(CredibleRegions.size());

  std::vector<std::vector<std::unique_ptr<TH2D>>> hpost_2D_copy(nDraw);
  std::vector<std::vector<std::vector<std::unique_ptr<TH2D>>>> hpost_2D_cl(nDraw);
  //KS: Copy all histograms to be thread safe
  for (int i = 0; i < nDraw; ++i)
  {
    hpost_2D_copy[i].resize(nDraw);
    hpost_2D_cl[i].resize(nDraw);
    for (int j = 0; j <= i; ++j)
    {
      hpost_2D_copy[i][j] = M3::Clone<TH2D>(hpost2D[i][j], Form("hpost_copy_%i_%i", i, j));
      hpost_2D_cl[i][j].resize(nCredible);
      for (int k = 0; k < nCredible; ++k)
      {
        hpost_2D_cl[i][j][k] = M3::Clone<TH2D>(hpost2D[i][j], Form("hpost_copy_%i_%i_CL_%f", i, j, CredibleRegions[k]));
      }
    }
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  //Calculate credible histogram
  for (int i = 0; i < nDraw; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      for (int k = 0; k < nCredible; ++k)
      {
        GetCredibleRegionSig(hpost_2D_cl[i][j][k], CredibleInSigmas, CredibleRegions[k]);
        hpost_2D_cl[i][j][k]->SetLineColor(CredibleRegionColor[k]);
        hpost_2D_cl[i][j][k]->SetLineWidth(2);
        hpost_2D_cl[i][j][k]->SetLineStyle(CredibleRegionStyle[k]);
      }
    }
  }

  gStyle->SetPalette(51);
  for (int i = 0; i < nDraw; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      // Skip the diagonal elements which we've already done above
      if (j == i) continue;
      if (IamVaried[j] == false) continue;

      auto legend = std::make_unique<TLegend>(0.20, 0.7, 0.4, 0.92);
      legend->SetTextColor(kRed);
      SetLegendStyle(legend.get(), 0.03);

      //Get Best point
      auto bestfitM = std::make_unique<TGraph>(1);
      const int MaxBin = hpost_2D_copy[i][j]->GetMaximumBin();
      int Mbx, Mby, Mbz;
      hpost_2D_copy[i][j]->GetBinXYZ(MaxBin, Mbx, Mby, Mbz);
      const double Mx = hpost_2D_copy[i][j]->GetXaxis()->GetBinCenter(Mbx);
      const double My = hpost_2D_copy[i][j]->GetYaxis()->GetBinCenter(Mby);

      bestfitM->SetPoint(0, Mx, My);
      bestfitM->SetMarkerStyle(22);
      bestfitM->SetMarkerSize(1);
      bestfitM->SetMarkerColor(kMagenta);
    
      //Plot default 2D posterior

      if(Draw2DPosterior){
      hpost_2D_copy[i][j]->Draw("COLZ");
      }
      else{
      hpost_2D_copy[i][j]->Draw("AXIS");
      }

      //Now credible regions
      for (int k = 0; k < nCredible; ++k)
        hpost_2D_cl[i][j][k]->Draw("CONT3 SAME");
      for (int k = nCredible-1; k >= 0; --k)
      {
        if(CredibleInSigmas)
          legend->AddEntry(hpost_2D_cl[i][j][k].get(), Form("%.0f#sigma Credible Interval", CredibleRegions[k]), "l");
        else
          legend->AddEntry(hpost_2D_cl[i][j][k].get(), Form("%.0f%% Credible Region", CredibleRegions[k]*100), "l");
      }
      legend->Draw("SAME");
  
    if(DrawBestFit){
      legend->AddEntry(bestfitM.get(),"Best Fit","p");
      bestfitM->Draw("SAME.P");
     }

      // Write to file
      Posterior->SetName(hpost2D[i][j]->GetName());
      Posterior->SetTitle(hpost2D[i][j]->GetTitle());

      //KS: Print only regions with correlation greater than specified value, by default 0.2. This is done to avoid dumping thousands of plots
      if(printToPDF && std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold) Posterior->Print(CanvasName);
      // Write it to root file
      //OutputFile->cd();
      //if( std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold ) Posterior->Write();
    }
  }

  OutputFile->cd();
}

// *********************
// Make fancy triangle plot for selected parameters
void MCMCProcessor::MakeTrianglePlot(const std::vector<std::string>& ParNames,
                                     // 1D
                                     const std::vector<double>& CredibleIntervals,
                                     const std::vector<Color_t>& CredibleIntervalsColours,
                                     //2D
                                     const std::vector<double>& CredibleRegions,
                                     const std::vector<Style_t>& CredibleRegionStyle,
                                     const std::vector<Color_t>& CredibleRegionColor,
                                     // Other
                                     const bool CredibleInSigmas) {
// *********************
  if(hpost2D.size() == 0) MakeCovariance_MP();
  MACH3LOG_INFO("Making Triangle Plot");

  const int nParamPlot = int(ParNames.size());
  std::vector<int> ParamNumber;
  for(int j = 0; j < nParamPlot; ++j)
  {
    int ParamNo = GetParamIndexFromName(ParNames[j]);
    if(ParamNo == M3::_BAD_INT_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Will not plot Triangle plot", ParNames[j]);
      return;
    }
    ParamNumber.push_back(ParamNo);
  }

  //KS: Store it as we go back to them at the end
  const std::vector<double> Margins = GetMargins(Posterior);
  Posterior->SetTopMargin(0.001);
  Posterior->SetBottomMargin(0.001);
  Posterior->SetLeftMargin(0.001);
  Posterior->SetRightMargin(0.001);

  // KS: We later format hist several times so make one unfired lambda
  auto FormatHistogram = [](auto& hist) {
    hist->GetXaxis()->SetTitle("");
    hist->GetYaxis()->SetTitle("");
    hist->SetTitle("");

    hist->GetXaxis()->SetLabelSize(0.1);
    hist->GetYaxis()->SetLabelSize(0.1);

    hist->GetXaxis()->SetNdivisions(4);
    hist->GetYaxis()->SetNdivisions(4);
  };

  Posterior->cd();
  Posterior->Clear();
  Posterior->Update();

  //KS: We sort to have parameters from highest to lowest, this is related to how we make 2D projections in MakeCovariance_MP
  std::sort(ParamNumber.begin(), ParamNumber.end(),  std::greater<int>());

  //KS: Calculate how many pads/plots we need
  int Npad = 0;
  for(int j = 1; j < nParamPlot+1; j++) Npad += j;
  Posterior->cd();
  // KS: Sanity check of size and ordering is correct
  CheckCredibleIntervalsOrder(CredibleIntervals, CredibleIntervalsColours);
  CheckCredibleRegionsOrder(CredibleRegions, CredibleRegionStyle, CredibleRegionColor);

  const int nCredibleIntervals = int(CredibleIntervals.size());
  const int nCredibleRegions = int(CredibleRegions.size());

  //KS: Initialise Tpad histograms etc we will need
  std::vector<TPad*> TrianglePad(Npad);
  //KS: 1D copy of posterior, we need it as we modify them
  std::vector<std::unique_ptr<TH1D>> hpost_copy(nParamPlot);
  std::vector<std::vector<std::unique_ptr<TH1D>>> hpost_cl(nParamPlot);
  std::vector<std::unique_ptr<TText>> TriangleText(nParamPlot * 2);
  std::vector<std::unique_ptr<TH2D>> hpost_2D_copy(Npad-nParamPlot);
  std::vector<std::vector<std::unique_ptr<TH2D>>> hpost_2D_cl(Npad-nParamPlot);
  gStyle->SetPalette(51);

  //KS: Super convoluted way of calculating ranges for our pads, trust me it works...
  std::vector<double> X_Min(nParamPlot);
  std::vector<double> X_Max(nParamPlot);
  X_Min[0] = 0.10;
  double xScale = (0.95 - (X_Min[0]+0.05))/nParamPlot;
  //KS: 0.05 is because we need additional offset for labels
  X_Max[0] = X_Min[0]+xScale+0.05;
  for(int i = 1; i < nParamPlot; i++)
  {
    X_Min[i] = X_Max[i-1];
    X_Max[i] = X_Min[i]+xScale;
  }
  std::vector<double> Y_Min(nParamPlot);
  std::vector<double> Y_Max(nParamPlot);
  Y_Max[0] = 0.95;
  //KS: 0.10 is becasue we need additional offset for labels
  double yScale = std::fabs(0.10 - (Y_Max[0]))/nParamPlot;
  Y_Min[0] = Y_Max[0]-yScale;
  for(int i = 1; i < nParamPlot; i++)
  {
    Y_Max[i] = Y_Min[i-1];
    Y_Min[i] = Y_Max[i]-yScale;
  }

  //KS: We store as numbering of isn't straightforward
  int counterPad = 0, counterText = 0, counterPost = 0, counter2DPost = 0;
  //KS: We start from top of the plot, might be confusing but works very well
  for(int y = 0; y < nParamPlot; y++)
  {
    //KS: start from left and go right, depending on y
    for(int x = 0; x <= y; x++)
    {
      //KS: Need to go to canvas every time to have our pads in the same canvas, not pads in the pads
      Posterior->cd();
      TrianglePad[counterPad] = new TPad(Form("TPad_%i", counterPad), Form("TPad_%i", counterPad),
                                         X_Min[x], Y_Min[y], X_Max[x], Y_Max[y]);

      TrianglePad[counterPad]->SetTopMargin(0);
      TrianglePad[counterPad]->SetRightMargin(0);

      TrianglePad[counterPad]->SetGrid();
      TrianglePad[counterPad]->SetFrameBorderMode(0);
      TrianglePad[counterPad]->SetBorderMode(0);
      TrianglePad[counterPad]->SetBorderSize(0);

      //KS: Corresponds to bottom part of the plot, need margins for labels
      TrianglePad[counterPad]->SetBottomMargin(y == (nParamPlot - 1) ? 0.1 : 0);
      //KS: Corresponds to left part, need margins for labels
      TrianglePad[counterPad]->SetLeftMargin(x == 0 ? 0.15 : 0);

      TrianglePad[counterPad]->Draw();
      TrianglePad[counterPad]->cd();

      //KS:if diagonal plot main posterior
      if(x == y)
      {
        hpost_copy[counterPost] = M3::Clone<TH1D>(hpost[ParamNumber[x]], Form("hpost_copy_%i", ParamNumber[x]));
        hpost_cl[counterPost].resize(nCredibleIntervals);
        /// Scale the histograms so it shows the posterior probability
        hpost_copy[counterPost]->Scale(1. / hpost_copy[counterPost]->Integral());
        for (int j = 0; j < nCredibleIntervals; ++j)
        {
          hpost_cl[counterPost][j] = M3::Clone<TH1D>(hpost[ParamNumber[x]], Form("hpost_copy_%i_CL_%f", ParamNumber[x], CredibleIntervals[j]));
          //KS: Reset to get rid to TF1 otherwise we run into segfault :(
          hpost_cl[counterPost][j]->Reset("");
          hpost_cl[counterPost][j]->Fill(0.0, 0.0);

          // Scale the histograms before gettindg credible intervals
          hpost_cl[counterPost][j]->Scale(1. / hpost_cl[counterPost][j]->Integral());
          GetCredibleIntervalSig(hpost_copy[counterPost], hpost_cl[counterPost][j], CredibleInSigmas, CredibleIntervals[j]);

          hpost_cl[counterPost][j]->SetFillColor(CredibleIntervalsColours[j]);
          hpost_cl[counterPost][j]->SetLineWidth(1);
        }

        hpost_copy[counterPost]->SetMaximum(hpost_copy[counterPost]->GetMaximum()*1.2);
        hpost_copy[counterPost]->SetLineWidth(2);
        hpost_copy[counterPost]->SetLineColor(kBlack);

        //KS: Don't want any titles
        FormatHistogram(hpost_copy[counterPost]);

        hpost_copy[counterPost]->Draw("HIST");
        for (int j = 0; j < nCredibleIntervals; ++j){
          hpost_cl[counterPost][j]->Draw("HIST SAME");
        }
        counterPost++;
      }
      //KS: Here we plot 2D credible regions
      else
      {
        hpost_2D_copy[counter2DPost] = M3::Clone<TH2D>(hpost2D[ParamNumber[x]][ParamNumber[y]],
                                                       Form("hpost_copy_%i_%i", ParamNumber[x], ParamNumber[y]));
        hpost_2D_cl[counter2DPost].resize(nCredibleRegions);
        //KS: Now copy for every credible region
        for (int k = 0; k < nCredibleRegions; ++k)
        {
          hpost_2D_cl[counter2DPost][k] = M3::Clone<TH2D>(hpost2D[ParamNumber[x]][ParamNumber[y]],
                                                          Form("hpost_copy_%i_%i_CL_%f", ParamNumber[x], ParamNumber[y], CredibleRegions[k]));
          GetCredibleRegionSig(hpost_2D_cl[counter2DPost][k], CredibleInSigmas, CredibleRegions[k]);

          hpost_2D_cl[counter2DPost][k]->SetLineColor(CredibleRegionColor[k]);
          hpost_2D_cl[counter2DPost][k]->SetLineWidth(2);
          hpost_2D_cl[counter2DPost][k]->SetLineStyle(CredibleRegionStyle[k]);
        }
        //KS: Don't want any titles
        FormatHistogram(hpost_2D_copy[counter2DPost]);

        hpost_2D_copy[counter2DPost]->Draw("COL");
        //Now credible regions
        for (int k = 0; k < nCredibleRegions; ++k){
          hpost_2D_cl[counter2DPost][k]->Draw("CONT3 SAME");
        }
        counter2DPost++;
      }
      //KS: Corresponds to bottom part of the plot
      if(y == (nParamPlot-1))
      {
        Posterior->cd();
        TriangleText[counterText] = std::make_unique<TText>(X_Min[x]+ (X_Max[x]-X_Min[x])/4, 0.04, hpost[ParamNumber[x]]->GetTitle());
        //KS: Unfortunately for many plots or long names this can go out of bounds :(
        TriangleText[counterText]->SetTextSize(0.015);
        TriangleText[counterText]->SetNDC(true);
        TriangleText[counterText]->Draw();
        counterText++;
      }
      //KS: Corresponds to left part
      if(x == 0)
      {
        Posterior->cd();
        TriangleText[counterText] = std::make_unique<TText>(0.04, Y_Min[y] + (Y_Max[y]-Y_Min[y])/4, hpost[ParamNumber[y]]->GetTitle());
        //KS: Rotate as this is y axis
        TriangleText[counterText]->SetTextAngle(90);
        //KS: Unfortunately for many plots or long names this can go out of bounds :(
        TriangleText[counterText]->SetTextSize(0.015);
        TriangleText[counterText]->SetNDC(true);
        TriangleText[counterText]->Draw();
        counterText++;
      }
      Posterior->Update();
      counterPad++;
    }
  }

  Posterior->cd();
  auto legend = std::make_unique<TLegend>(0.60, 0.7, 0.9, 0.9);
  SetLegendStyle(legend.get(), 0.03);
  //KS: Legend is shared so just take first histograms
  for (int j = nCredibleIntervals-1; j >= 0; --j)
  {
    if(CredibleInSigmas)
      legend->AddEntry(hpost_cl[0][j].get(), Form("%.0f#sigma Credible Interval", CredibleIntervals[j]), "f");
    else
      legend->AddEntry(hpost_cl[0][j].get(), Form("%.0f%% Credible Interval", CredibleRegions[j]*100), "f");
  }
  for (int k = nCredibleRegions-1; k >= 0; --k)
  {
    if(CredibleInSigmas)
      legend->AddEntry(hpost_2D_cl[0][k].get(), Form("%.0f#sigma Credible Region", CredibleRegions[k]), "l");
    else
      legend->AddEntry(hpost_2D_cl[0][k].get(), Form("%.0f%% Credible Region", CredibleRegions[k]*100), "l");
  }
  legend->Draw("SAME");
  Posterior->Update();

  // Write to file
  Posterior->SetName("TrianglePlot");
  Posterior->SetTitle("TrianglePlot");

  if(printToPDF) Posterior->Print(CanvasName);
  // Write it to root file
  OutputFile->cd();
  Posterior->Write();

  //KS: Remove allocated structures
  for(int i = 0; i < Npad; i++) delete TrianglePad[i];

  //KS: Restore margin
  SetMargins(Posterior, Margins);
}

// **************************
// Scan the input trees
void MCMCProcessor::ScanInput() {
// **************************
  // KS: This can reduce time necessary for caching even by half
  #ifdef MULTITHREAD
  //ROOT::EnableImplicitMT();
  #endif

  // Open the Chain
  Chain = new TChain("posteriors","posteriors");
  Chain->Add(MCMCFile.c_str());

  nEntries = int(Chain->GetEntries());
  
  //Only is suboptimality we might want to change it, therefore set it high enough so it doesn't affect other functionality
  UpperCut = nEntries+1;

  // Get the list of branches
  TObjArray* brlis = Chain->GetListOfBranches();

  // Get the number of branches
  nBranches = brlis->GetEntries();

  BranchNames.reserve(nBranches);
  ParamType.reserve(nBranches);

  // Read the input Covariances
  ReadInputCov();

  // Set all the branches to off
  Chain->SetBranchStatus("*", false);

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  for (int i = 0; i < nBranches; i++)
  {
    // Get the TBranch and its name
    TBranch* br = static_cast<TBranch*>(brlis->At(i));
    if(!br){
      MACH3LOG_ERROR("Invalid branch at position {}", i);
      throw MaCh3Exception(__FILE__,__LINE__);
    }
    TString bname = br->GetName();

    //KS: Exclude parameter types
    bool rejected = false;
    for(unsigned int ik = 0; ik < ExcludedTypes.size(); ++ik )
    {
      if(bname.BeginsWith(ExcludedTypes[ik]))
      {
        rejected = true;
        break;
      }
    }
    if(rejected) continue;

    // Turn on the branches which we want for parameters
    Chain->SetBranchStatus(bname.Data(), true);

    if (bname.BeginsWith("ndd_"))
    {
      BranchNames.push_back(bname);
      ParamType.push_back(kNDPar);
      nParam[kNDPar]++;
    }
    else if (bname.BeginsWith("skd_joint_"))
    {
      BranchNames.push_back(bname);
      ParamType.push_back(kFDDetPar);
      nParam[kFDDetPar]++;
    }

    //KS: as a bonus get LogL systematic
    if (bname.BeginsWith("LogL_sample_")) {
      SampleName_v.push_back(bname);
      nSamples++;
    }
    else if (bname.BeginsWith("LogL_systematic_")) {
      SystName_v.push_back(bname);
      nSysts++;
    }
  }
  nDraw = int(BranchNames.size());

  // Read the input Covariances
  ReadInputCovLegacy();
  
  // Check order of parameter types
  ScanParameterOrder();

  IamVaried.resize(nDraw, true);

  // Print useful Info
  PrintInfo();

  nSteps = Chain->GetMaximum("step");
  // Set the step cut to be 20%
  int cut = nSteps/5;
  SetStepCut(cut);

  // Basically allow loading oscillation parameters
  LoadAdditionalInfo();
}

// ****************************
// Set up the output files and canvases
void MCMCProcessor::SetupOutput() {
// ****************************
  // Make sure we can read files located anywhere and strip the .root ending
  MCMCFile = MCMCFile.substr(0, MCMCFile.find(".root"));

  // Check if the output file is ready
  if (OutputFile == nullptr) MakeOutputFile();
  
  CanvasName = MCMCFile + OutputSuffix + ".pdf[";
  if(printToPDF) Posterior->Print(CanvasName);

  // Once the pdf file is open no longer need to bracket
  CanvasName.ReplaceAll("[","");

  // We fit with this Gaussian
  Gauss = std::make_unique<TF1>("Gauss", "[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])", -5, 5);
  Gauss->SetLineWidth(2);
  Gauss->SetLineColor(kOrange-5);

  // Declare the TVectors
  Covariance = new TMatrixDSym(nDraw);
  Correlation = new TMatrixDSym(nDraw);
  Central_Value = new TVectorD(nDraw);
  Means = new TVectorD(nDraw);
  Errors = new TVectorD(nDraw);
  Means_Gauss = new TVectorD(nDraw);
  Errors_Gauss = new TVectorD(nDraw);
  Means_HPD    = new TVectorD(nDraw);
  Errors_HPD   = new TVectorD(nDraw);
  Errors_HPD_Positive = new TVectorD(nDraw);
  Errors_HPD_Negative = new TVectorD(nDraw);

  // Initialise to something silly
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nDraw; ++i)
  {
    (*Central_Value)(i) = M3::_BAD_DOUBLE_;
    (*Means)(i) = M3::_BAD_DOUBLE_;
    (*Errors)(i) = M3::_BAD_DOUBLE_;
    (*Means_Gauss)(i) = M3::_BAD_DOUBLE_;
    (*Errors_Gauss)(i) = M3::_BAD_DOUBLE_;
    (*Means_HPD)(i) = M3::_BAD_DOUBLE_;
    (*Errors_HPD)(i) = M3::_BAD_DOUBLE_;
    (*Errors_HPD_Positive)(i) = M3::_BAD_DOUBLE_;
    (*Errors_HPD_Negative)(i) = M3::_BAD_DOUBLE_;
    for (int j = 0; j < nDraw; ++j) {
      (*Covariance)(i, j) = M3::_BAD_DOUBLE_;
      (*Correlation)(i, j) = M3::_BAD_DOUBLE_;
    }
  }
  hpost.resize(nDraw);
}

// ****************************
// Check order of parameter types
void MCMCProcessor::ScanParameterOrder() {
// *****************************
  for(int i = 0; i < kNParameterEnum; i++)
  {
    for(unsigned int j = 0; j < ParamType.size(); j++)
    {
      if(ParamType[j] == ParameterEnum(i)) 
      {
        //KS: When we find that i-th parameter types start at j, save and move to the next parameter.
        ParamTypeStartPos[i] = j;
        break;
      }
    }
  }
}
    
// *****************************
// Make the prefit plots
std::unique_ptr<TH1D> MCMCProcessor::MakePrefit() {
// *****************************
  if (OutputFile == nullptr) MakeOutputFile();

  auto PreFitPlot = std::make_unique<TH1D>("Prefit", "Prefit", nDraw, 0, nDraw);
  PreFitPlot->SetDirectory(nullptr);
  for (int i = 0; i < PreFitPlot->GetNbinsX() + 1; ++i) {
    PreFitPlot->SetBinContent(i+1, 0);
    PreFitPlot->SetBinError(i+1, 0);
  }

  //KS: Slightly hacky way to get relative to prior or nominal as this is convention we use,
  //Only applies for xsec, for other systematic it make no difference
  double CentralValueTemp, Central, Error;

  // Set labels and data
  for (int i = 0; i < nDraw; ++i)
  {
    //Those keep which parameter type we run currently and relative number
    int ParamEnum = ParamType[i];
    int ParamNo = i - ParamTypeStartPos[ParameterEnum(ParamEnum)];
    CentralValueTemp = ParamCentral[ParamEnum][ParamNo];
    if(plotRelativeToPrior) 
    {
      // Normalise the prior relative the nominal/prior, just the way we get our fit results in MaCh3
      if ( CentralValueTemp != 0) {
        Central = ParamCentral[ParamEnum][ParamNo] / CentralValueTemp;
        Error = ParamErrors[ParamEnum][ParamNo]/CentralValueTemp;
      } else {
        Central = CentralValueTemp + 1.0;
        Error = ParamErrors[ParamEnum][ParamNo];
      }
    }
    else
    {
      Central = CentralValueTemp;
      Error = ParamErrors[ParamEnum][ParamNo];
    }
    //KS: If plotting error for param with flat prior is turned off and given param really has flat prior set error to 0
    if(!PlotFlatPrior && ParamFlat[ParamEnum][ParamNo]) {
      Error = 0.;
    }
    PreFitPlot->SetBinContent(i+1, Central);
    PreFitPlot->SetBinError(i+1, Error);
    PreFitPlot->GetXaxis()->SetBinLabel(i+1, ParamNames[ParamEnum][ParamNo]);
  }
  PreFitPlot->SetDirectory(nullptr);

  PreFitPlot->SetFillStyle(1001);
  PreFitPlot->SetFillColor(kRed-3);
  PreFitPlot->SetMarkerStyle(21);
  PreFitPlot->SetMarkerSize(2.4);
  PreFitPlot->SetMarkerColor(kWhite);
  PreFitPlot->SetLineColor(PreFitPlot->GetFillColor());
  PreFitPlot->GetXaxis()->LabelsOption("v");

  return PreFitPlot;
}

// **************************
//CW: Read the input Covariance matrix entries
// Get stuff like parameter input errors, names, and so on
void MCMCProcessor::ReadInputCov() {
// **************************
  FindInputFiles();
  if(CovPos[kXSecPar].back() != "none") ReadModelFile();
}

// **************************
//CW: Read the input Covariance matrix entries
// Get stuff like parameter input errors, names, and so on
void MCMCProcessor::ReadInputCovLegacy() {
// **************************
  FindInputFilesLegacy();
  if(nParam[kNDPar] > 0)    ReadNDFile();
  if(nParam[kFDDetPar] > 0) ReadFDFile();
}

// **************************
// Read the output MCMC file and find what inputs were used
void MCMCProcessor::FindInputFiles() {
// **************************
  // Now read the MCMC file
  TFile *TempFile = new TFile(MCMCFile.c_str(), "open");
  TDirectory* CovarianceFolder = TempFile->Get<TDirectory>("CovarianceFolder");

  // Get the settings for the MCMC
  TMacro *Config = TempFile->Get<TMacro>("MaCh3_Config");

  if (Config == nullptr) {
    MACH3LOG_ERROR("Didn't find MaCh3_Config tree in MCMC file! {}", MCMCFile);
    TempFile->ls();
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  //KS:Most inputs are in ${MACH3}/inputs/blarb.root
  if (std::getenv("MACH3") != nullptr) {
    MACH3LOG_INFO("Found MACH3 environment variable: {}", std::getenv("MACH3"));
  }

  MACH3LOG_INFO("Loading YAML config from MCMC chain");

  YAML::Node Settings = TMacroToYAML(*Config);

  bool InputNotFound = false;
  //CW: Get the xsec Covariance matrix
  CovPos[kXSecPar] = GetFromManager<std::vector<std::string>>(Settings["General"]["Systematics"]["XsecCovFile"], {"none"});
  if(CovPos[kXSecPar].back() == "none")
  {
    MACH3LOG_WARN("Couldn't find XsecCov branch in output");
    InputNotFound = true;
  }

  TMacro *XsecConfig = M3::GetConfigMacroFromChain(CovarianceFolder);
  if (XsecConfig == nullptr) {
    MACH3LOG_WARN("Didn't find Config_xsec_cov tree in MCMC file! {}", MCMCFile);
  } else {
    CovConfig[kXSecPar] = TMacroToYAML(*XsecConfig);
  }
  if(InputNotFound) MaCh3Utils::PrintConfig(Settings);

  if (const char * mach3_env = std::getenv("MACH3"))
  {
    for(size_t i = 0; i < CovPos[kXSecPar].size(); i++)
      CovPos[kXSecPar][i].insert(0, std::string(mach3_env)+"/");
  }

  // Delete the TTrees and the input file handle since we've now got the settings we need
  delete Config;
  delete XsecConfig;

  // Delete the MCMCFile pointer we're reading
  CovarianceFolder->Close();
  delete CovarianceFolder;
  TempFile->Close();
  delete TempFile;
}

// **************************
// Read the output MCMC file and find what inputs were used
void MCMCProcessor::FindInputFilesLegacy() {
// **************************
  // Now read the MCMC file
  TFile *TempFile = new TFile(MCMCFile.c_str(), "open");

  // Get the settings for the MCMC
  TMacro *Config = TempFile->Get<TMacro>("MaCh3_Config");

  if (Config == nullptr) {
    MACH3LOG_ERROR("Didn't find MaCh3_Config tree in MCMC file! {}", MCMCFile);
    TempFile->ls();
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  YAML::Node Settings = TMacroToYAML(*Config);

  //CW: And the ND Covariance matrix
  CovPos[kNDPar].push_back(GetFromManager<std::string>(Settings["General"]["Systematics"]["NDCovFile"], "none"));

  if(CovPos[kNDPar].back() == "none") {
    MACH3LOG_WARN("Couldn't find NDCov (legacy) branch in output");
  } else{
    //If the FD Cov is not none, then you need the name of the covariance object to grab
    CovNamePos[kNDPar].push_back(GetFromManager<std::string>(Settings["General"]["Systematics"]["NDCovName"], "none"));
    MACH3LOG_INFO("Given NDCovFile {} and NDCovName {}", CovPos[kNDPar].back(), CovNames[kNDPar].back());
  }

  //CW: And the FD Covariance matrix
  CovPos[kFDDetPar].push_back(GetFromManager<std::string>(Settings["General"]["Systematics"]["FDCovFile"], "none"));

  if(CovPos[kFDDetPar].back() == "none") {
    MACH3LOG_WARN("Couldn't find FDCov (legacy) branch in output");
  } else {
    //If the FD Cov is not none, then you need the name of the covariance object to grab
    CovNamePos[kFDDetPar].push_back(GetFromManager<std::string>(Settings["General"]["Systematics"]["FDCovName"], "none"));
    MACH3LOG_INFO("Given FDCovFile {} and FDCovName {}", CovPos[kFDPar].back(), CovNames[kFDPar].back());
  }

  if (const char * mach3_env = std::getenv("MACH3"))
  {
    for(size_t i = 0; i < CovPos[kNDPar].size(); i++)
      CovPos[kNDPar][i].insert(0, std::string(mach3_env)+"/");

    for(size_t i = 0; i < CovPos[kFDDetPar].size(); i++)
      CovPos[kFDDetPar][i].insert(0, std::string(mach3_env)+"/");
  }
  TempFile->Close();
  delete TempFile;
}

// ***************
// Read the model file and get the input central values and errors
void MCMCProcessor::ReadModelFile() {
// ***************
  YAML::Node XSecFile = CovConfig[kXSecPar];

  auto systematics = XSecFile["Systematics"];
  int paramIndex  = 0;
  for (auto it = systematics.begin(); it != systematics.end(); ++it, ++paramIndex )
  {
    auto const &param = *it;
    // Push back the name
    std::string TempString = (param["Systematic"]["Names"]["FancyName"].as<std::string>());

    bool rejected = false;
    for (unsigned int ik = 0; ik < ExcludedNames.size(); ++ik)
    {
      if (TempString.rfind(ExcludedNames.at(ik), 0) == 0)
      {
        rejected = true;
        break;
      }
    }
    if(rejected) continue;

    ParamNames[kXSecPar].push_back(TempString);
    ParamCentral[kXSecPar].push_back(param["Systematic"]["ParameterValues"]["PreFitValue"].as<double>());
    ParamNom[kXSecPar].push_back(param["Systematic"]["ParameterValues"]["Generated"].as<double>());
    ParamErrors[kXSecPar].push_back(param["Systematic"]["Error"].as<double>() );
    ParamFlat[kXSecPar].push_back(GetFromManager<bool>(param["Systematic"]["FlatPrior"], false));

    ParameterGroup.push_back(param["Systematic"]["ParameterGroup"].as<std::string>());

    nParam[kXSecPar]++;
    ParamType.push_back(kXSecPar);
    // Params from osc group have branch name equal to fancy name while all others are basically xsec_0 for example
    if(ParameterGroup.back() == "Osc") {
      BranchNames.push_back(ParamNames[kXSecPar].back());
    } else {
      BranchNames.push_back("xsec_" + std::to_string(paramIndex));
    }

    // Check that the branch exists before setting address
    if (!Chain->GetBranch(BranchNames.back())) {
      MACH3LOG_ERROR("Couldn't find branch '{}'", BranchNames.back());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
}

// ***************
// Read the ND cov file and get the input central values and errors
void MCMCProcessor::ReadNDFile() {
// ***************
  // Do the same for the ND280
  TFile *NDdetFile = new TFile(CovPos[kNDPar].back().c_str(), "open");
  if (NDdetFile->IsZombie()) {
    MACH3LOG_ERROR("Couldn't find NDdetFile {}", CovPos[kNDPar].back());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  NDdetFile->cd();

  TMatrixDSym *NDdetMatrix = NDdetFile->Get<TMatrixDSym>(CovNamePos[kNDPar].back().c_str());
  TVectorD *NDdetNominal = NDdetFile->Get<TVectorD>("det_weights");
  TDirectory *BinningDirectory = NDdetFile->Get<TDirectory>("Binning");

  for (int i = 0; i < NDdetNominal->GetNrows(); ++i)
  {
    ParamNom[kNDPar].push_back( (*NDdetNominal)(i) );
    ParamCentral[kNDPar].push_back( (*NDdetNominal)(i) );

    ParamErrors[kNDPar].push_back( std::sqrt((*NDdetMatrix)(i,i)) );
    ParamNames[kNDPar].push_back( Form("ND Det %i", i) );
    //KS: Currently we can only set it via config, change it in future
    ParamFlat[kNDPar].push_back( false );
  }

  TIter next(BinningDirectory->GetListOfKeys());
  TKey *key = nullptr;
  // Loop through all entries
  while ((key = static_cast<TKey*>(next())))
  {
    std::string name = std::string(key->GetName());
    TH2Poly* RefPoly = BinningDirectory->Get<TH2Poly>((name).c_str());
    int size = RefPoly->GetNumberOfBins();
    NDSamplesBins.push_back(size);
    NDSamplesNames.push_back(RefPoly->GetTitle());
  }

  NDdetFile->Close();
  delete NDdetFile;
}

// ***************
// Read the FD cov file and get the input central values and errors
void MCMCProcessor::ReadFDFile() {
// ***************
  // Do the same for the FD
  TFile *FDdetFile = new TFile(CovPos[kFDDetPar].back().c_str(), "open");
  if (FDdetFile->IsZombie()) {
    MACH3LOG_ERROR("Couldn't find FDdetFile {}", CovPos[kFDDetPar].back());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  FDdetFile->cd();

  TMatrixD *FDdetMatrix = FDdetFile->Get<TMatrixD>(CovNamePos[kFDDetPar].back().c_str());

  for (int i = 0; i < FDdetMatrix->GetNrows(); ++i)
  {
    //KS: FD parameters start at 1. in contrary to ND280
    ParamNom[kFDDetPar].push_back(1.);
    ParamCentral[kFDDetPar].push_back(1.);

    ParamErrors[kFDDetPar].push_back( std::sqrt((*FDdetMatrix)(i,i)) );
    ParamNames[kFDDetPar].push_back( Form("FD Det %i", i) );

    //KS: Currently we can only set it via config, change it in future
    ParamFlat[kFDDetPar].push_back( false );
  }
  //KS: The last parameter is p scale
  //ETA: we need to be careful here, this is only true for SK in the T2K beam analysis...
  if(FancyPlotNames) ParamNames[kFDDetPar].back() = "Momentum Scale";

  FDdetFile->Close();
  delete FDdetFile;
  delete FDdetMatrix;
}

// ***************
// Make the step cut from a string
void MCMCProcessor::SetStepCut(const std::string& Cuts) {
// ***************
  StepCut = Cuts;
  BurnInCut = std::stoi( Cuts );
}

// ***************
// Make the step cut from an int
void MCMCProcessor::SetStepCut(const int Cuts) {
// ***************
  std::stringstream TempStream;
  TempStream << "step > " << Cuts;
  StepCut = TempStream.str();
  BurnInCut = Cuts;
}

// ***************
// Pass central value
void MCMCProcessor::GetNthParameter(const int param, double &Prior, double &PriorError, TString &Title) const {
// **************************
  ParameterEnum ParType = ParamType[param];
  int ParamNo = M3::_BAD_INT_;
  ParamNo = param - ParamTypeStartPos[ParType];

  Prior = ParamCentral[ParType][ParamNo];
  PriorError = ParamErrors[ParType][ParamNo];
  Title = ParamNames[ParType][ParamNo];
}

// ***************
// Find Param Index based on name
int MCMCProcessor::GetParamIndexFromName(const std::string& Name){
// **************************
  int ParamNo = M3::_BAD_INT_;
  for (int i = 0; i < nDraw; ++i)
  {
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(i, Prior, PriorError, Title);

    if(Name == Title)
    {
      ParamNo = i;
      break;
    }
  }
  return ParamNo;
}

// **************************************************
// Helper function to reset histograms
void MCMCProcessor::ResetHistograms() {
// **************************************************
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nDraw; ++i) 
  {
    for (int j = 0; j <= i; ++j)
    {
      // TH2D to hold the Correlation 
      hpost2D[i][j]->Reset("");
      hpost2D[i][j]->Fill(0.0, 0.0, 0.0);
    }
  }
}

// **************************
// KS: Get Super Fancy Polar Plot
void MCMCProcessor::GetPolarPlot(const std::vector<std::string>& ParNames){
// **************************
  if(hpost[0] == nullptr) MakePostfit();

  std::vector<double> Margins = GetMargins(Posterior);

  Posterior->SetTopMargin(0.1);
  Posterior->SetBottomMargin(0.1);
  Posterior->SetLeftMargin(0.1);
  Posterior->SetRightMargin(0.1);
  Posterior->Update();

  MACH3LOG_INFO("Calculating Polar Plot");
  TDirectory *PolarDir = OutputFile->mkdir("PolarDir");
  PolarDir->cd();

  for(unsigned int k = 0; k < ParNames.size(); ++k)
  {
    //KS: First we need to find parameter number based on name
    int ParamNo = GetParamIndexFromName(ParNames[k]);
    if(ParamNo == M3::_BAD_INT_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Will not calculate Polar Plot", ParNames[k]);
      continue;
    }

    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(ParamNo, Prior, PriorError, Title);

    std::vector<double> x_val(nBins);
    std::vector<double> y_val(nBins);

    double xmin = 0;
    double xmax = 2*TMath::Pi();

    double Integral = hpost[ParamNo]->Integral();
    for (Int_t ipt = 0; ipt < nBins; ipt++)
    {
      x_val[ipt] = ipt*(xmax-xmin)/nBins+xmin;
      y_val[ipt] = hpost[ParamNo]->GetBinContent(ipt+1)/Integral;
    }

    auto PolarGraph = std::make_unique<TGraphPolar>(nBins, x_val.data(), y_val.data());
    PolarGraph->SetLineWidth(2);
    PolarGraph->SetFillStyle(3001);
    PolarGraph->SetLineColor(kRed);
    PolarGraph->SetFillColor(kRed);
    PolarGraph->Draw("AFL");

    auto Text = std::make_unique<TText>(0.6, 0.1, Title);
    Text->SetTextSize(0.04);
    Text->SetNDC(true);
    Text->Draw("");

    Posterior->Print(CanvasName);
    Posterior->Write(Title);
  } //End loop over parameters

  PolarDir->Close();
  delete PolarDir;

  OutputFile->cd();

  SetMargins(Posterior, Margins);
}

// **************************
// Get Bayes Factor for particular parameter
void MCMCProcessor::GetBayesFactor(const std::vector<std::string>& ParNames,
                                   const std::vector<std::vector<double>>& Model1Bounds,
                                   const std::vector<std::vector<double>>& Model2Bounds,
                                   const std::vector<std::vector<std::string>>& ModelNames){
// **************************
  if(hpost[0] == nullptr) MakePostfit();

  MACH3LOG_INFO("Calculating Bayes Factor");
  if((ParNames.size() != Model1Bounds.size()) || (Model2Bounds.size() != Model1Bounds.size())  || (Model2Bounds.size() != ModelNames.size()))
  {
    MACH3LOG_ERROR("Size doesn't match");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  for(unsigned int k = 0; k < ParNames.size(); ++k)
  {
    //KS: First we need to find parameter number based on name
    int ParamNo = GetParamIndexFromName(ParNames[k]);
    if(ParamNo == M3::_BAD_INT_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Will not calculate Bayes Factor", ParNames[k]);
      continue;
    }

    const double M1_min = Model1Bounds[k][0];
    const double M2_min = Model2Bounds[k][0];
    const double M1_max = Model1Bounds[k][1];
    const double M2_max = Model2Bounds[k][1];

    long double IntegralMode1 = hpost[ParamNo]->Integral(hpost[ParamNo]->FindFixBin(M1_min), hpost[ParamNo]->FindFixBin(M1_max));
    long double IntegralMode2 = hpost[ParamNo]->Integral(hpost[ParamNo]->FindFixBin(M2_min), hpost[ParamNo]->FindFixBin(M2_max));

    double BayesFactor = 0.;
    std::string Name = "";
    //KS: Calc Bayes Factor
    //If M1 is more likely
    if(IntegralMode1 >= IntegralMode2)
    {
      BayesFactor = IntegralMode1/IntegralMode2;
      Name = "\\mathfrak{B}(" + ModelNames[k][0]+ "/" + ModelNames[k][1] + ") = " + std::to_string(BayesFactor);
    }
    else //If M2 is more likely
    {
      BayesFactor = IntegralMode2/IntegralMode1;
      Name = "\\mathfrak{B}(" + ModelNames[k][1]+ "/" + ModelNames[k][0] + ") = " + std::to_string(BayesFactor);
    }
    std::string JeffreysScale = GetJeffreysScale(BayesFactor);
    std::string DunneKabothScale = GetDunneKaboth(BayesFactor);

    MACH3LOG_INFO("{} for {}", Name, ParNames[k]);
    MACH3LOG_INFO("Following Jeffreys Scale = {}", JeffreysScale);
    MACH3LOG_INFO("Following Dunne-Kaboth Scale = {}", DunneKabothScale);
    MACH3LOG_INFO("");
  }
}

// **************************
// KS: Get Savage Dickey point hypothesis test
void MCMCProcessor::GetSavageDickey(const std::vector<std::string>& ParNames,
                                    const std::vector<double>& EvaluationPoint,
                                    const std::vector<std::vector<double>>& Bounds){
// **************************
  if((ParNames.size() != EvaluationPoint.size()) || (Bounds.size() != EvaluationPoint.size()))
  {
    MACH3LOG_ERROR("Size doesn't match");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  
  if(hpost[0] == nullptr) MakePostfit();

  MACH3LOG_INFO("Calculating Savage Dickey");
  TDirectory *SavageDickeyDir = OutputFile->mkdir("SavageDickey");
  SavageDickeyDir->cd();
  
  for(unsigned int k = 0; k < ParNames.size(); ++k)
  {
    //KS: First we need to find parameter number based on name
    int ParamNo = GetParamIndexFromName(ParNames[k]);
    if(ParamNo == M3::_BAD_INT_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Will not calculate SavageDickey", ParNames[k]);
      continue;
    }
    
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    bool FlatPrior = false;
    GetNthParameter(ParamNo, Prior, PriorError, Title);
    
    ParameterEnum ParType = ParamType[ParamNo];
    int ParamTemp = ParamNo - ParamTypeStartPos[ParType];
    FlatPrior = ParamFlat[ParType][ParamTemp];
    
    auto PosteriorHist = M3::Clone<TH1D>(hpost[ParamNo], std::string(Title));
    RemoveFitter(PosteriorHist.get(), "Gauss");
            
    std::unique_ptr<TH1D> PriorHist;
    //KS: If flat prior we need to have well defined bounds otherwise Prior distribution will not make sense
    if(FlatPrior)
    {
      int NBins = PosteriorHist->GetNbinsX();
      if(Bounds[k][0] > Bounds[k][1])
      {
        MACH3LOG_ERROR("Lower bound is higher than upper bound");
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
      PriorHist = std::make_unique<TH1D>("PriorHist", Title, NBins, Bounds[k][0], Bounds[k][1]);
      
      double FlatProb = ( Bounds[k][1] - Bounds[k][0]) / NBins;
      for (int g = 0; g < NBins + 1; ++g) 
      {
        PriorHist->SetBinContent(g+1, FlatProb);
      }
    }
    else //KS: Otherwise throw from Gaussian
    {
      PriorHist = M3::Clone<TH1D>(PosteriorHist.get(), "Prior");
      PriorHist->Reset("");
      PriorHist->Fill(0.0, 0.0);
      
      auto rand = std::make_unique<TRandom3>(0);
      //KS: Throw nice gaussian, just need big number to have smooth distribution
      for(int g = 0; g < 1000000; ++g)
      {
        PriorHist->Fill(rand->Gaus(Prior, PriorError));
      }
    }
    SavageDickeyPlot(PriorHist, PosteriorHist, std::string(Title), EvaluationPoint[k]);
  } //End loop over parameters

  SavageDickeyDir->Close();
  delete SavageDickeyDir;

  OutputFile->cd();
}

// **************************
// KS: Get Savage Dickey point hypothesis test
void MCMCProcessor::SavageDickeyPlot(std::unique_ptr<TH1D>& PriorHist,
                                     std::unique_ptr<TH1D>& PosteriorHist,
                                     const std::string& Title,
                                     const double EvaluationPoint) const {
// **************************
  // Area normalise the distributions
  PriorHist->Scale(1./PriorHist->Integral(), "width");
  PosteriorHist->Scale(1./PosteriorHist->Integral(), "width");

  PriorHist->SetLineColor(kRed);
  PriorHist->SetMarkerColor(kRed);
  PriorHist->SetFillColorAlpha(kRed, 0.35);
  PriorHist->SetFillStyle(1001);
  PriorHist->GetXaxis()->SetTitle(Title.c_str());
  PriorHist->GetYaxis()->SetTitle("Posterior Probability");
  PriorHist->SetMaximum(PosteriorHist->GetMaximum()*1.5);
  PriorHist->GetYaxis()->SetLabelOffset(999);
  PriorHist->GetYaxis()->SetLabelSize(0);
  PriorHist->SetLineWidth(2);
  PriorHist->SetLineStyle(kSolid);

  PosteriorHist->SetLineColor(kBlue);
  PosteriorHist->SetMarkerColor(kBlue);
  PosteriorHist->SetFillColorAlpha(kBlue, 0.35);
  PosteriorHist->SetFillStyle(1001);

  PriorHist->Draw("hist");
  PosteriorHist->Draw("hist same");

  double ProbPrior = PriorHist->GetBinContent(PriorHist->FindBin(EvaluationPoint));
  //KS: In case we go so far away that prior is 0, set this to small value to avoid dividing by 0
  if(ProbPrior < 0) ProbPrior = 0.00001;
  double ProbPosterior = PosteriorHist->GetBinContent(PosteriorHist->FindBin(EvaluationPoint));
  double SavageDickey = ProbPosterior/ProbPrior;

  std::string DunneKabothScale = GetDunneKaboth(SavageDickey);
  //Get Best point
  std::unique_ptr<TGraph> PostPoint(new TGraph(1));
  PostPoint->SetPoint(0, EvaluationPoint, ProbPosterior);
  PostPoint->SetMarkerStyle(20);
  PostPoint->SetMarkerSize(1);
  PostPoint->Draw("P same");

  std::unique_ptr<TGraph> PriorPoint(new TGraph(1));
  PriorPoint->SetPoint(0, EvaluationPoint, ProbPrior);
  PriorPoint->SetMarkerStyle(20);
  PriorPoint->SetMarkerSize(1);
  PriorPoint->Draw("P same");

  auto legend = std::make_unique<TLegend>(0.12, 0.6, 0.6, 0.97);
  SetLegendStyle(legend.get(), 0.04);
  legend->AddEntry(PriorHist.get(), "Prior", "l");
  legend->AddEntry(PosteriorHist.get(), "Posterior", "l");
  legend->AddEntry(PostPoint.get(), Form("SavageDickey = %.2f, (%s)", SavageDickey, DunneKabothScale.c_str()),"");
  legend->Draw("same");

  Posterior->Print(CanvasName);
  Posterior->Write(Title.c_str());
}

// **************************
// KS: Reweight prior of MCMC chain to another
void MCMCProcessor::ReweightPrior(const std::vector<std::string>& Names,
                                  const std::vector<double>& NewCentral,
                                  const std::vector<double>& NewError) {
// **************************
  MACH3LOG_INFO("Reweighting Prior");

  if( (Names.size() != NewCentral.size()) || (NewCentral.size() != NewError.size()))
  {
    MACH3LOG_ERROR("Size of passed vectors doesn't match in ReweightPrior");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  std::vector<int> Param;
  std::vector<double> OldCentral;
  std::vector<double> OldError;
  std::vector<bool> FlatPrior;

  //KS: First we need to find parameter number based on name
  for(unsigned int k = 0; k < Names.size(); ++k)
  {
    //KS: First we need to find parameter number based on name
    int ParamNo = GetParamIndexFromName(Names[k]);
    if(ParamNo == M3::_BAD_INT_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Can't reweight Prior", Names[k]);
      return;
    }

    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(ParamNo, Prior, PriorError, Title);

    Param.push_back(ParamNo);
    OldCentral.push_back(Prior);
    OldError.push_back(PriorError);

    ParameterEnum ParType = ParamType[ParamNo];
    int ParamTemp = ParamNo - ParamTypeStartPos[ParType];

    FlatPrior.push_back(ParamFlat[ParType][ParamTemp]);
  }
  std::vector<double> ParameterPos(Names.size());

  std::string InputFile = MCMCFile+".root";
  std::string OutputFilename = MCMCFile + "_reweighted.root";

  //KS: Simply create copy of file and add there new branch
  int ret = system(("cp " + InputFile + " " + OutputFilename).c_str());
  if (ret != 0)
    MACH3LOG_WARN("Error: system call to copy file failed with code {}", ret);

  TFile *OutputChain = new TFile(OutputFilename.c_str(), "UPDATE");
  OutputChain->cd();
  TTree *post = OutputChain->Get<TTree>("posteriors");

  double Weight = 1.;

  post->SetBranchStatus("*",false);
  // Set the branch addresses for params
  for (unsigned int j = 0; j < Names.size(); ++j) {
    post->SetBranchStatus(BranchNames[Param[j]].Data(), true);
    post->SetBranchAddress(BranchNames[Param[j]].Data(), &ParameterPos[j]);
  }
  TBranch *bpt = post->Branch("Weight", &Weight, "Weight/D");
  post->SetBranchStatus("Weight", true);

  for (int i = 0; i < nEntries; ++i)
  {
    post->GetEntry(i);
    Weight = 1.;

    //KS: Calculate reweight weight. Weights are multiplicative so we can do several reweights at once. FIXME Big limitation is that code only works for uncorrelated parameters :(
    for (unsigned int j = 0; j < Names.size(); ++j)
    {
      double new_chi = (ParameterPos[j] - NewCentral[j])/NewError[j];
      double new_prior = std::exp(-0.5 * new_chi * new_chi);

      double old_chi = -1;
      double old_prior = -1;
      if(FlatPrior[j]) {
        old_prior = 1.0;
      } else {
        old_chi = (ParameterPos[j] - OldCentral[j])/OldError[j];
        old_prior = std::exp(-0.5 * old_chi * old_chi);
      }
      Weight *= new_prior/old_prior;
    }
    bpt->Fill();
  }
  post->SetBranchStatus("*",true);
  OutputChain->cd();
  post->Write("posteriors", TObject::kOverwrite);

  // KS: Save reweight metadeta
  std::ostringstream yaml_stream;
  yaml_stream << "Weight:\n";
  for (size_t k = 0; k < Names.size(); ++k) {
    yaml_stream << "    " << Names[k] << ": [" << NewCentral[k] << ", " << NewError[k] << "]\n";
  }
  std::string yaml_string = yaml_stream.str();
  YAML::Node root = STRINGtoYAML(yaml_string);
  TMacro ConfigSave = YAMLtoTMacro(root, "Reweight_Config");
  ConfigSave.Write();

  OutputChain->Close();
  delete OutputChain;

  OutputFile->cd();
}

// **************************
// Diagnose the MCMC
void MCMCProcessor::ParameterEvolution(const std::vector<std::string>& Names,
                                       const std::vector<int>& NIntervals) {
// **************************
  MACH3LOG_INFO("Parameter Evolution gif");

  //KS: First we need to find parameter number based on name
  for(unsigned int k = 0; k < Names.size(); ++k)
  {
    //KS: First we need to find parameter number based on name
    int ParamNo = GetParamIndexFromName(Names[k]);
    if(ParamNo == M3::_BAD_INT_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Can't reweight Prior", Names[k]);
      continue;
    }

    const int IntervalsSize = nSteps/NIntervals[k];
    // ROOT won't overwrite gifs so we need to delete the file if it's there already
    int ret = system(fmt::format("rm {}.gif",Names[k]).c_str());
    if (ret != 0){
      MACH3LOG_WARN("Error: system call to delete {} failed with code {}", Names[k], ret);
    }

    // This holds the posterior density
    const double maxi = Chain->GetMaximum(BranchNames[ParamNo]);
    const double mini = Chain->GetMinimum(BranchNames[ParamNo]);

    int Counter = 0;
    for(int i = NIntervals[k]-1; i >= 0; --i)
    {
      // This holds the posterior density
      TH1D* EvePlot = new TH1D(BranchNames[ParamNo], BranchNames[ParamNo], nBins, mini, maxi);
      EvePlot->SetMinimum(0);
      EvePlot->GetYaxis()->SetTitle("PDF");
      EvePlot->GetYaxis()->SetNoExponent(false);

      //KS: Apply additional Cuts, like mass ordering
      std::string CutPosterior1D = "step > " + std::to_string(i*IntervalsSize+IntervalsSize);

      std::string TextTitle = "Steps = 0 - "+std::to_string(Counter*IntervalsSize+IntervalsSize);
      // Project BranchNames[ParamNo] onto hpost, applying stepcut
      Chain->Project(BranchNames[ParamNo], BranchNames[ParamNo], CutPosterior1D.c_str());

      EvePlot->SetLineWidth(2);
      EvePlot->SetLineColor(kBlue-1);
      EvePlot->SetTitle(Names[k].c_str());
      EvePlot->GetXaxis()->SetTitle(EvePlot->GetTitle());
      EvePlot->GetYaxis()->SetLabelOffset(1000);
      if(ApplySmoothing) EvePlot->Smooth();

      EvePlot->Scale(1. / EvePlot->Integral());
      EvePlot->Draw("HIST");

      TText text(0.3, 0.8, TextTitle.c_str());
      text.SetTextFont (43);
      text.SetTextSize (40);
      text.SetNDC(true);
      text.Draw("SAME");

      if(i == 0) Posterior->Print((Names[k] + ".gif++20").c_str()); // produces infinite loop animated GIF
      else Posterior->Print((Names[k] + ".gif+20").c_str()); // add picture to .gif

      delete EvePlot;
      Counter++;
    }
  }
}

// **************************
// Diagnose the MCMC
void MCMCProcessor::DiagMCMC() {
// **************************
  // Prepare branches etc for DiagMCMC
  PrepareDiagMCMC();

  // Draw the simple trace matrices
  ParamTraces();

  // Get the batched means
  BatchedMeans();

  // Draw the auto-correlations
  if (useFFTAutoCorrelation) {
    AutoCorrelation_FFT();
  } else {
    AutoCorrelation();
  }

  // Calculate Power Spectrum for each param
  PowerSpectrumAnalysis();

  // Get Geweke Z score helping select burn-in
  GewekeDiagnostic();

  // Draw acceptance Probability
  AcceptanceProbabilities();
}

// **************************
//CW: Prepare branches etc. for DiagMCMC
void MCMCProcessor::PrepareDiagMCMC() {
// **************************
  doDiagMCMC = true;
    
  if(ParStep != nullptr) {
    MACH3LOG_ERROR("It look like ParStep was already filled ");
    MACH3LOG_ERROR("Even though it is used for MakeCovariance_MP and for DiagMCMC");
    MACH3LOG_ERROR("it has different structure in both for cache hits, sorry ");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  if(nBatches == 0) {
    MACH3LOG_ERROR("nBatches is equal to 0");
    MACH3LOG_ERROR("please use SetnBatches to set other value fore example 20");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
    
  // Initialise ParStep
  ParStep = new double*[nDraw]();
  for (int j = 0; j < nDraw; ++j) {
    ParStep[j] = new double[nEntries]();
    for (int i = 0; i < nEntries; ++i) {
      ParStep[j][i] = -999.99;
    }
  }

  SampleValues = new double*[nEntries]();
  SystValues = new double*[nEntries]();
  AccProbValues = new double[nEntries]();
  StepNumber = new int[nEntries]();
  for (int i = 0; i < nEntries; ++i) {
    SampleValues[i] = new double[nSamples]();
    SystValues[i] = new double[nSysts]();

    for (int j = 0; j < nSamples; ++j) {
      SampleValues[i][j] = -999.99;
    }
    for (int j = 0; j < nSysts; ++j) {
      SystValues[i][j] = -999.99;
    }
    AccProbValues[i] = -999.99;
    StepNumber[i] = -999.99;
  }

  // Initialise the sums
  ParamSums = new double[nDraw]();
  for (int i = 0; i < nDraw; ++i) {
    ParamSums[i] = 0.0;
  }
  MACH3LOG_INFO("Reading input tree...");
  TStopwatch clock;
  clock.Start();

  // Set all the branches to off
  Chain->SetBranchStatus("*", false);

  // 10 entries output
  const int countwidth = nEntries/10;

  // Can also do the batched means here to minimize excessive loops
  // The length of each batch
  const int BatchLength = nEntries/nBatches+1;
  BatchedAverages = new double*[nBatches]();
  AccProbBatchedAverages = new double[nBatches]();
  for (int i = 0; i < nBatches; ++i) {
    BatchedAverages[i] = new double[nDraw];
    AccProbBatchedAverages[i] = 0;
    for (int j = 0; j < nDraw; ++j) {
      BatchedAverages[i][j] = 0.0;
    }
  }
  std::vector<double> ParStepBranch(nDraw);
  std::vector<double> SampleValuesBranch(nSamples);
  std::vector<double> SystValuesBranch(nSysts);
  int StepNumberBranch = 0;
  double AccProbValuesBranch = 0;
  // Set the branch addresses for params
  for (int j = 0; j < nDraw; ++j) {
    Chain->SetBranchStatus(BranchNames[j].Data(), true);
    Chain->SetBranchAddress(BranchNames[j].Data(), &ParStepBranch[j]);
  }
  // Set the branch addresses for samples
  for (int j = 0; j < nSamples; ++j) {
    Chain->SetBranchStatus(SampleName_v[j].Data(), true);
    Chain->SetBranchAddress(SampleName_v[j].Data(), &SampleValuesBranch[j]);
  }
  // Set the branch addresses for systematics
  for (int j = 0; j < nSysts; ++j) {
    Chain->SetBranchStatus(SystName_v[j].Data(), true);
    Chain->SetBranchAddress(SystName_v[j].Data(), &SystValuesBranch[j]);
  }
  // Only needed for Geweke right now
  Chain->SetBranchStatus("step", true);
  Chain->SetBranchAddress("step", &StepNumberBranch);
  // Turn on the branches which we want for acc prob
  Chain->SetBranchStatus("accProb", true);
  Chain->SetBranchAddress("accProb", &AccProbValuesBranch);

  // Loop over the entries
  //KS: This is really a bottleneck right now, thus revisit with ROOT6 https://pep-root6.github.io/docs/analysis/parallell/root.html
  for (int i = 0; i < nEntries; ++i) {
    // Fill up the arrays
    Chain->GetEntry(i);

    if (i % countwidth == 0)
      MaCh3Utils::PrintProgressBar(i, nEntries);

    // Set the branch addresses for params
    for (int j = 0; j < nDraw; ++j) {
      ParStep[j][i] = ParStepBranch[j];
    }
    // Set the branch addresses for samples
    for (int j = 0; j < nSamples; ++j) {
      SampleValues[i][j] = SampleValuesBranch[j];
    }
    // Set the branch addresses for systematics
    for (int j = 0; j < nSysts; ++j) {
      SystValues[i][j] = SystValuesBranch[j];
    }
      
    // Set the branch addresses for Acceptance Probability
    AccProbValues[i] = AccProbValuesBranch;
    StepNumber[i] = StepNumberBranch;

    // Find which batch the event belongs in
    int BatchNumber = -1;
    // I'm so lazy! But it's OK, the major overhead here is GetEntry: saved by ROOT!
    for (int j = 0; j < nBatches; ++j) {
      if (i < (j+1)*BatchLength) {
        BatchNumber = j;
        break;
      }
    }
    // Fill up the sum for each j param
    for (int j = 0; j < nDraw; ++j) {
      ParamSums[j] += ParStep[j][i];
      BatchedAverages[BatchNumber][j] += ParStep[j][i];
    }
    
    //KS: Could easily add this to above loop but I accProb is different beast so better keep it like this
    AccProbBatchedAverages[BatchNumber] += AccProbValues[i];
  }
  clock.Stop();
  MACH3LOG_INFO("Took {:.2f}s to finish caching statistic for Diag MCMC with {} steps", clock.RealTime(), nEntries);

  // Make the sums into average
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nDraw; ++i) {
    ParamSums[i] /= double(nEntries);
    for (int j = 0; j < nBatches; ++j) {
      // Divide by the total number of events in the batch
      BatchedAverages[j][i] /= BatchLength;
      if(i == 0) AccProbBatchedAverages[j] /= BatchLength; //KS: we have only one accProb, keep it like this for now
    }
  }

  // And make our sweet output file
  if (OutputFile == nullptr) MakeOutputFile();
}

// *****************
//CW: Draw trace plots of the parameters i.e. parameter vs step
void MCMCProcessor::ParamTraces() {
// *****************
  if (ParStep == nullptr) PrepareDiagMCMC();
  MACH3LOG_INFO("Making trace plots...");
  // Make the TH1Ds
  std::vector<TH1D*> TraceParamPlots(nDraw);
  std::vector<TH1D*> TraceSamplePlots(nSamples);
  std::vector<TH1D*> TraceSystsPlots(nSysts);

  // Set the titles and limits for TH2Ds
  for (int j = 0; j < nDraw; ++j) {
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    
    GetNthParameter(j, Prior, PriorError, Title);
    std::string HistName = Form("%s_%s_Trace", Title.Data(), BranchNames[j].Data());
    TraceParamPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nEntries, 0, nEntries);
    TraceParamPlots[j]->GetXaxis()->SetTitle("Step");
    TraceParamPlots[j]->GetYaxis()->SetTitle("Parameter Variation");
  }

  for (int j = 0; j < nSamples; ++j) {
    std::string HistName = SampleName_v[j].Data();
    TraceSamplePlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nEntries, 0, nEntries);
    TraceSamplePlots[j]->GetXaxis()->SetTitle("Step");
    TraceSamplePlots[j]->GetYaxis()->SetTitle("Sample -logL");
  }

  for (int j = 0; j < nSysts; ++j) {
    std::string HistName = SystName_v[j].Data();
    TraceSystsPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nEntries, 0, nEntries);
    TraceSystsPlots[j]->GetXaxis()->SetTitle("Step");
    TraceSystsPlots[j]->GetYaxis()->SetTitle("Systematic -logL");
  }

  // Have now made the empty TH1Ds, now for writing content to them!
  // Loop over the number of parameters to draw their traces
  // Each histogram
#ifdef MULTITHREAD
  MACH3LOG_INFO("Using multi-threading...");
  #pragma omp parallel for
#endif
  for (int i = 0; i < nEntries; ++i) {
    // Set bin content for the ith bin to the parameter values
    for (int j = 0; j < nDraw; ++j) {
      TraceParamPlots[j]->SetBinContent(i, ParStep[j][i]);
    }
    for (int j = 0; j < nSamples; ++j) {
      TraceSamplePlots[j]->SetBinContent(i, SampleValues[i][j]);
    }
    for (int j = 0; j < nSysts; ++j) {
      TraceSystsPlots[j]->SetBinContent(i, SystValues[i][j]);
    }
  }

  // Write the output and delete the TH2Ds
  TDirectory *TraceDir = OutputFile->mkdir("Trace");
  TraceDir->cd();
  for (int j = 0; j < nDraw; ++j) {
    // Fit a linear function to the traces
    TF1 *Fitter = new TF1("Fitter","[0]", nEntries/2, nEntries);
    Fitter->SetLineColor(kRed);
    TraceParamPlots[j]->Fit("Fitter","Rq");
    TraceParamPlots[j]->Write();
    delete Fitter;
    delete TraceParamPlots[j];
  }

  TDirectory *LLDir = OutputFile->mkdir("LogL");
  LLDir->cd();
  for (int j = 0; j < nSamples; ++j) {
    TraceSamplePlots[j]->Write();
    delete TraceSamplePlots[j];
    delete[] SampleValues[j];
  }
  delete[] SampleValues;

  for (int j = 0; j < nSysts; ++j) {
    TraceSystsPlots[j]->Write();
    delete TraceSystsPlots[j];
    delete SystValues[j];
  }
  delete[] SystValues;

  TraceDir->Close();
  delete TraceDir;

  OutputFile->cd();
}

// *********************************
// MJR: Calculate autocorrelations using the FFT algorithm.
//      Fast, even on CPU, and get all lags for free.
void MCMCProcessor::AutoCorrelation_FFT() {
// *********************************
  if (ParStep == nullptr) PrepareDiagMCMC();

  TStopwatch clock;
  clock.Start();
  const int nLags = AutoCorrLag;
  MACH3LOG_INFO("Making auto-correlations for nLags = {}", nLags);

  // Prep outputs
  OutputFile->cd();
  TDirectory* AutoCorrDir = OutputFile->mkdir("Auto_corr");
  std::vector<TH1D*> LagKPlots(nDraw);
  std::vector<std::vector<double>> LagL(nDraw);

  // Arrays needed to perform FFT using ROOT
  std::vector<double> ACFFT(nEntries, 0.0);  // Main autocorrelation array
  std::vector<double> ParVals(nEntries, 0.0);  // Param values for full chain
  std::vector<double> ParValsFFTR(nEntries, 0.0);  // FFT Real part
  std::vector<double> ParValsFFTI(nEntries, 0.0);  // FFT Imaginary part
  std::vector<double> ParValsFFTSquare(nEntries, 0.0);  // FFT Absolute square
  std::vector<double> ParValsComplex(nEntries, 0.0);  // Input Imaginary values (0)

  // Create forward and reverse FFT objects. I don't love using ROOT here,
  // but it works so I can't complain
  TVirtualFFT* fftf = TVirtualFFT::FFT(1, &nEntries, "C2CFORWARD");
  TVirtualFFT* fftb = TVirtualFFT::FFT(1, &nEntries, "C2CBACKWARD");

  // Loop over all pars and calculate the full autocorrelation function using FFT
  for (int j = 0; j < nDraw; ++j) {
    // Initialize
    LagL[j].resize(nLags);
    for (int i = 0; i < nEntries; ++i) {
      ParVals[i] = ParStep[j][i]-ParamSums[j]; // Subtract the mean to make it numerically tractable
      ParValsComplex[i] = 0.; // Reset dummy array
    }

    // Transform
    fftf->SetPointsComplex(ParVals.data(), ParValsComplex.data());
    fftf->Transform();
    fftf->GetPointsComplex(ParValsFFTR.data(), ParValsFFTI.data());

    // Square the results to get the power spectrum
    for (int i = 0; i < nEntries; ++i) {
      ParValsFFTSquare[i] = ParValsFFTR[i]*ParValsFFTR[i] + ParValsFFTI[i]*ParValsFFTI[i];
    }

    // Transforming back gives the autocovariance
    fftb->SetPointsComplex(ParValsFFTSquare.data(), ParValsComplex.data());
    fftb->Transform();
    fftb->GetPointsComplex(ACFFT.data(), ParValsComplex.data());

    // Divide by norm to get autocorrelation
    double normAC = ACFFT[0];
    for (int i = 0; i < nEntries; ++i) {
      ACFFT[i] /= normAC;
    }

    // Get plotting info
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(j, Prior, PriorError, Title);
    std::string HistName = Form("%s_%s_Lag", Title.Data(), BranchNames[j].Data());

    // Initialize Lag plot
    LagKPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nLags, 0.0, nLags);
    LagKPlots[j]->GetXaxis()->SetTitle("Lag");
    LagKPlots[j]->GetYaxis()->SetTitle("Auto-correlation function");

    // Fill plot
    for (int k = 0; k < nLags; ++k) {
      LagL[j][k] = ACFFT[k];
      LagKPlots[j]->SetBinContent(k, ACFFT[k]);
    }

    // Write and clean up
    AutoCorrDir->cd();
    LagKPlots[j]->Write();
    delete LagKPlots[j];
  }

  //KS: This is different diagnostic however it relies on calculated Lag, thus we call it before we delete LagKPlots
  CalculateESS(nLags, LagL);

  AutoCorrDir->Close();
  delete AutoCorrDir;

  OutputFile->cd();

  clock.Stop();
  MACH3LOG_INFO("Making auto-correlations took {:.2f}s", clock.RealTime());
}

// *********************************
//KS: Calculate autocorrelations supports both OpenMP and CUDA :)
void MCMCProcessor::AutoCorrelation() {
// *********************************
  if (ParStep == nullptr) PrepareDiagMCMC();

  TStopwatch clock;
  clock.Start();
  const int nLags = AutoCorrLag;
  MACH3LOG_INFO("Making auto-correlations for nLags = {}", nLags);

  // The sum of (Y-Ymean)^2 over all steps for each parameter
  std::vector<std::vector<double>> DenomSum(nDraw);
  std::vector<std::vector<double>> NumeratorSum(nDraw);
  std::vector<std::vector<double>> LagL(nDraw);
  for (int j = 0; j < nDraw; ++j) {
    DenomSum[j].resize(nLags);
    NumeratorSum[j].resize(nLags);
    LagL[j].resize(nLags);
  }
  std::vector<TH1D*> LagKPlots(nDraw);
  // Loop over the parameters of interest
  for (int j = 0; j < nDraw; ++j)
  {
    // Loop over each lag
    for (int k = 0; k < nLags; ++k) {
      NumeratorSum[j][k] = 0.0;
      DenomSum[j][k] = 0.0;
      LagL[j][k] = 0.0;
    }

    // Make TH1Ds for each parameter which hold the lag
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    
    GetNthParameter(j, Prior, PriorError, Title);
    std::string HistName = Form("%s_%s_Lag", Title.Data(), BranchNames[j].Data());
    LagKPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nLags, 0.0, nLags);
    LagKPlots[j]->GetXaxis()->SetTitle("Lag");
    LagKPlots[j]->GetYaxis()->SetTitle("Auto-correlation function");
  }
//KS: If CUDA is not enabled do calculations on CPU
#ifndef MaCh3_CUDA
  // Loop over the lags
  //CW: Each lag is independent so might as well multi-thread them!
  #ifdef MULTITHREAD
  MACH3LOG_INFO("Using multi-threading...");
  #pragma omp parallel for collapse(2)
  #endif      // Loop over the number of parameters
  for (int j = 0; j < nDraw; ++j) {
    for (int k = 0; k < nLags; ++k) {
      // Loop over the number of entries
      for (int i = 0; i < nEntries; ++i) {
        const double Diff = ParStep[j][i]-ParamSums[j];

        // Only sum the numerator up to i = N-k
        if (i < nEntries-k) {
          const double LagTerm = ParStep[j][i+k]-ParamSums[j];
          const double Product = Diff*LagTerm;
          NumeratorSum[j][k] += Product;
        }
        // Square the difference to form the denominator
        const double Denom = Diff*Diff;
        DenomSum[j][k] += Denom;
      }
    }
  }
#else //NOW GPU specific code
  MACH3LOG_INFO("Using GPU");
  //KS: This allocates memory and copy data from CPU to GPU
  PrepareGPU_AutoCorr(nLags);

  //KS: This runs the main kernel and copy results back to CPU
  RunGPU_AutoCorr(
    ParStep_gpu,
    ParamSums_gpu,
    NumeratorSum_gpu,
    DenomSum_gpu,
    NumeratorSum_cpu,
    DenomSum_cpu);

  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2)
  #endif
  //KS: Now that that we received data from GPU convert it to CPU-like format
  for (int j = 0; j < nDraw; ++j)
  {
    for (int k = 0; k < nLags; ++k)
    {
      const int temp_index = j*nLags+k;
      NumeratorSum[j][k] = NumeratorSum_cpu[temp_index];
      DenomSum[j][k] = DenomSum_cpu[temp_index];
    }
  }
  //delete auxiliary variables
  delete[] NumeratorSum_cpu;
  delete[] DenomSum_cpu;
  delete[] ParStep_cpu;
  delete[] ParamSums_cpu;

  //KS: Delete stuff at GPU as well
  CleanupGPU_AutoCorr(
      ParStep_gpu,
      NumeratorSum_gpu,
      ParamSums_gpu,
      DenomSum_gpu);

//KS: End of GPU specific code
#endif

  OutputFile->cd();
  TDirectory *AutoCorrDir = OutputFile->mkdir("Auto_corr");
  // Now fill the LagK auto-correlation plots
  for (int j = 0; j < nDraw; ++j) {
    for (int k = 0; k < nLags; ++k) {
      LagL[j][k] = NumeratorSum[j][k]/DenomSum[j][k];
      LagKPlots[j]->SetBinContent(k, NumeratorSum[j][k]/DenomSum[j][k]);
    }
    AutoCorrDir->cd();
    LagKPlots[j]->Write();
    delete LagKPlots[j];
  }

  //KS: This is different diagnostic however it relies on calculated Lag, thus we call it before we delete LagKPlots
  CalculateESS(nLags, LagL);

  delete[] ParamSums;

  AutoCorrDir->Close();
  delete AutoCorrDir;

  OutputFile->cd();

  clock.Stop();
  MACH3LOG_INFO("Making auto-correlations took {:.2f}s", clock.RealTime());
}

#ifdef MaCh3_CUDA
// **************************
//KS: Allocates memory and copy data from CPU to GPU
void MCMCProcessor::PrepareGPU_AutoCorr(const int nLags) {
// **************************
  //KS: Create temporary arrays that will communicate with GPU code
  ParStep_cpu = new float[nDraw*nEntries];
  NumeratorSum_cpu = new float[nDraw*nLags];
  DenomSum_cpu = new float[nDraw*nLags];
  ParamSums_cpu = new float[nDraw];

  #ifdef MULTITHREAD
  //KS: Open parallel region
  #pragma omp parallel
  {
  #endif
    //KS: Operations are independent thus we are using nowait close
    #ifdef MULTITHREAD
    #pragma omp for nowait
    #endif
    for (int i = 0; i < nDraw; ++i)
    {
      //KS: We basically need this to convert from double to float for GPU
      ParamSums_cpu[i]  = ParamSums[i];
    }

    #ifdef MULTITHREAD
    #pragma omp for collapse(2) nowait
    #endif
    for (int j = 0; j < nDraw; ++j)
    {
      for (int k = 0; k < nLags; ++k)
      {
        const int temp = j*nLags+k;
        NumeratorSum_cpu[temp] = 0.0;
        DenomSum_cpu[temp] = 0.0;
      }
    }

    #ifdef MULTITHREAD
    #pragma omp for collapse(2)
    #endif
    for (int j = 0; j < nDraw; ++j)
    {
      for (int i = 0; i < nEntries; ++i)
      {
        const int temp = j*nEntries+i;
        ParStep_cpu[temp] = ParStep[j][i];
      }
    }
  #ifdef MULTITHREAD
  //KS: End parallel region
  }
  #endif

  //KS: First allocate memory on GPU
  InitGPU_AutoCorr(&ParStep_gpu,
                   &NumeratorSum_gpu,
                   &ParamSums_gpu,
                   &DenomSum_gpu,

                   nEntries,
                   nDraw,
                   nLags);


  //KS: Now copy from CPU to GPU
  CopyToGPU_AutoCorr(ParStep_cpu,
                     NumeratorSum_cpu,
                     ParamSums_cpu,
                     DenomSum_cpu,

                     ParStep_gpu,
                     NumeratorSum_gpu,
                     ParamSums_gpu,
                     DenomSum_gpu);
}
#endif


// **************************
// KS: calc Effective Sample Size Following @cite StanManual
// Furthermore we calculate Sampling efficiency following @cite hanson2008mcmc
// Rule of thumb is to have efficiency above 25%
void MCMCProcessor::CalculateESS(const int nLags, const std::vector<std::vector<double>>& LagL) {
// **************************
  if(LagL.size() == 0)
  {
    MACH3LOG_ERROR("Size of LagL is 0");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  MACH3LOG_INFO("Making ESS plots...");
  TVectorD* EffectiveSampleSize = new TVectorD(nDraw);
  TVectorD* SamplingEfficiency = new TVectorD(nDraw);
  std::vector<double> TempDenominator(nDraw);

  constexpr int Nhists = 5;
  constexpr double Thresholds[Nhists + 1] = {1, 0.02, 0.005, 0.001, 0.0001, 0.0};
  constexpr Color_t ESSColours[Nhists] = {kGreen, kGreen + 2, kYellow, kOrange, kRed};

  //KS: This histogram is inspired by the following: @cite gabry2024visual
  std::vector<std::unique_ptr<TH1D>> EffectiveSampleSizeHist(Nhists);
  for(int i = 0; i < Nhists; ++i)
  {
    EffectiveSampleSizeHist[i] =
      std::make_unique<TH1D>(Form("EffectiveSampleSizeHist_%i", i), Form("EffectiveSampleSizeHist_%i", i), nDraw, 0, nDraw);
    EffectiveSampleSizeHist[i]->GetYaxis()->SetTitle("N_{eff}/N");
    EffectiveSampleSizeHist[i]->SetFillColor(ESSColours[i]);
    EffectiveSampleSizeHist[i]->SetLineColor(ESSColours[i]);
    EffectiveSampleSizeHist[i]->Sumw2();
    for (int j = 0; j < nDraw; ++j)
    {
      TString Title = "";
      double Prior = 1.0, PriorError = 1.0;
      GetNthParameter(j, Prior, PriorError, Title);
      EffectiveSampleSizeHist[i]->GetXaxis()->SetBinLabel(j+1, Title.Data());
    }
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  //KS: Calculate ESS and MCMC efficiency for each parameter
  for (int j = 0; j < nDraw; ++j)
  {
    (*EffectiveSampleSize)(j) = M3::_BAD_DOUBLE_;
    (*SamplingEfficiency)(j) = M3::_BAD_DOUBLE_;
    TempDenominator[j] = 0.;
    //KS: Firs sum over all Calculated autocorrelations
    for (int k = 0; k < nLags; ++k)
    {
      TempDenominator[j] += LagL[j][k];
    }
    TempDenominator[j] = 1+2*TempDenominator[j];
    (*EffectiveSampleSize)(j) = double(nEntries)/TempDenominator[j];
    // 100 because we convert to percentage
    (*SamplingEfficiency)(j) = 100 * 1/TempDenominator[j];

    for(int i = 0; i < Nhists; ++i)
    {
      EffectiveSampleSizeHist[i]->SetBinContent(j+1, 0);
      EffectiveSampleSizeHist[i]->SetBinError(j+1, 0);

      const double TempEntry = std::fabs((*EffectiveSampleSize)(j)) / double(nEntries);
      if(Thresholds[i] >= TempEntry && TempEntry > Thresholds[i+1])
      {
        if( std::isnan((*EffectiveSampleSize)(j)) ) continue;
        EffectiveSampleSizeHist[i]->SetBinContent(j+1, TempEntry);
      }
    }
  }

  //KS Write to the output tree
  //Save to file
  OutputFile->cd();
  EffectiveSampleSize->Write("EffectiveSampleSize");
  SamplingEfficiency->Write("SamplingEfficiency");

  EffectiveSampleSizeHist[0]->SetTitle("Effective Sample Size");
  EffectiveSampleSizeHist[0]->Draw();
  for(int i = 1; i < Nhists; ++i)
  {
    EffectiveSampleSizeHist[i]->Draw("SAME");
  }

  auto leg = std::make_unique<TLegend>(0.2, 0.7, 0.6, 0.95);
  SetLegendStyle(leg.get(), 0.03);
  for(int i = 0; i < Nhists; ++i)
  {
    leg->AddEntry(EffectiveSampleSizeHist[i].get(), Form("%.4f >= N_{eff}/N > %.4f", Thresholds[i], Thresholds[i+1]), "f");
  }  leg->Draw("SAME");

  Posterior->Write("EffectiveSampleSizeCanvas");

  //Delete all variables
  delete EffectiveSampleSize;
  delete SamplingEfficiency;
}

// **************************
//CW: Batched means, literally read from an array and chuck into TH1D
void MCMCProcessor::BatchedMeans() {
// **************************
  if (BatchedAverages == nullptr) PrepareDiagMCMC();
  MACH3LOG_INFO("Making BatchedMeans plots...");
  
  std::vector<TH1D*> BatchedParamPlots(nDraw);
  for (int j = 0; j < nDraw; ++j) {
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    
    GetNthParameter(j, Prior, PriorError, Title);
    
    std::string HistName = Form("%s_%s_batch", Title.Data(), BranchNames[j].Data());
    BatchedParamPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nBatches, 0, nBatches);
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int j = 0; j < nDraw; ++j) {
    for (int i = 0; i < nBatches; ++i) {
      BatchedParamPlots[j]->SetBinContent(i+1, BatchedAverages[i][j]);
      const int BatchRangeLow = double(i)*double(nEntries)/double(nBatches);
      const int BatchRangeHigh = double(i+1)*double(nEntries)/double(nBatches);
      std::stringstream ss;
      ss << BatchRangeLow << " - " << BatchRangeHigh;
      BatchedParamPlots[j]->GetXaxis()->SetBinLabel(i+1, ss.str().c_str());
    }
  }

  TDirectory *BatchDir = OutputFile->mkdir("Batched_means");
  BatchDir->cd();
  for (int j = 0; j < nDraw; ++j) {
    TF1 *Fitter = new TF1("Fitter","[0]", 0, nBatches);
    Fitter->SetLineColor(kRed);
    BatchedParamPlots[j]->Fit("Fitter","Rq");
    BatchedParamPlots[j]->Write();
    delete Fitter;
    delete BatchedParamPlots[j];
  }

  //KS: Get the batched means variance estimation and variable indicating if number of batches is sensible
  // We do this before deleting BatchedAverages
  BatchedAnalysis();

  for (int i = 0; i < nBatches; ++i) {
    delete BatchedAverages[i];
  }
  delete[] BatchedAverages;

  BatchDir->Close();
  delete BatchDir;

  OutputFile->cd();
}

// **************************
// Get the batched means variance estimation and variable indicating if number of batches is sensible
void MCMCProcessor::BatchedAnalysis() {
// **************************
  if(BatchedAverages == nullptr)
  {
    MACH3LOG_ERROR("BatchedAverages haven't been initialises or have been deleted something is wrong");
    MACH3LOG_ERROR("I need it and refuse to go further");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Calculate variance estimator using batched means following @cite chakraborty2019estimating see Eq. 1.2
  TVectorD* BatchedVariance = new TVectorD(nDraw);
  //KS: The hypothesis is rejected if C > z α for a given confidence level α. If the batch means do not pass the test, Correlated is reported for the half-width on the statistical reports following @cite rossetti2024batch alternatively for more old-school see Alexopoulos and Seila 1998 section 3.4.3
  TVectorD* C_Test_Statistics = new TVectorD(nDraw);
 
  std::vector<double> OverallBatchMean(nDraw);
  std::vector<double> C_Rho_Nominator(nDraw);
  std::vector<double> C_Rho_Denominator(nDraw);
  std::vector<double> C_Nominator(nDraw);
  std::vector<double> C_Denominator(nDraw);
  const int BatchLength = nEntries/nBatches+1;
//KS: Start parallel region
#ifdef MULTITHREAD
#pragma omp parallel
{
#endif
  #ifdef MULTITHREAD
  #pragma omp for
  #endif
  //KS: First calculate mean of batched means for each param and Initialise everything to 0
  for (int j = 0; j < nDraw; ++j)
  {
    OverallBatchMean[j] = 0.0;
    C_Rho_Nominator[j] = 0.0;
    C_Rho_Denominator[j] = 0.0;
    C_Nominator[j] = 0.0;
    C_Denominator[j] = 0.0;

    (*BatchedVariance)(j) = 0.0;
    (*C_Test_Statistics)(j) = 0.0;
    for (int i = 0; i < nBatches; ++i)
    {
      OverallBatchMean[j] += BatchedAverages[i][j];
    }
    OverallBatchMean[j] /= nBatches;
  }

  #ifdef MULTITHREAD
  #pragma omp for nowait
  #endif
  //KS: next loop is completely independent thus nowait clause
  for (int j = 0; j < nDraw; ++j)
  {
    for (int i = 0; i < nBatches; ++i)
    {
      (*BatchedVariance)(j) += (OverallBatchMean[j] - BatchedAverages[i][j])*(OverallBatchMean[j] - BatchedAverages[i][j]);
    }
    (*BatchedVariance)(j) = (BatchLength/(nBatches-1))* (*BatchedVariance)(j);
  }
  
  //KS: Now we focus on C test statistic, again use nowait as next is calculation is independent
  #ifdef MULTITHREAD
  #pragma omp for nowait
  #endif
  for (int j = 0; j < nDraw; ++j)
  {
    C_Nominator[j] = (OverallBatchMean[j] - BatchedAverages[0][j])*(OverallBatchMean[j] - BatchedAverages[0][j]) + 
                      (OverallBatchMean[j] - BatchedAverages[nBatches-1][j])*(OverallBatchMean[j] - BatchedAverages[nBatches-1][j]);
    for (int i = 0; i < nBatches; ++i)
    {
      C_Denominator[j] += (OverallBatchMean[j] - BatchedAverages[i][j])*(OverallBatchMean[j] - BatchedAverages[i][j]);
    }
    C_Denominator[j] = 2*C_Denominator[j];
  }
  
  //KS: We still calculate C and for this we need rho wee need autocorrelations between batches
  #ifdef MULTITHREAD
  #pragma omp for
  #endif
  for (int j = 0; j < nDraw; ++j)
  {
    for (int i = 0; i < nBatches-1; ++i)
    {
      C_Rho_Nominator[j] += (OverallBatchMean[j] - BatchedAverages[i][j])*(OverallBatchMean[j] - BatchedAverages[i+1][j]);
    }
    
    for (int i = 0; i < nBatches; ++i)
    {
      C_Rho_Denominator[j] += (OverallBatchMean[j] - BatchedAverages[i][j])*(OverallBatchMean[j] - BatchedAverages[i][j]);
    }
  }
  
  //KS: Final calculations of C
  #ifdef MULTITHREAD
  #pragma omp for
  #endif
  for (int j = 0; j < nDraw; ++j)
  {
    (*C_Test_Statistics)(j) = std::sqrt((nBatches*nBatches - 1)/(nBatches-2)) * ( C_Rho_Nominator[j]/C_Rho_Denominator[j] + C_Nominator[j]/ C_Denominator[j]);
  }
#ifdef MULTITHREAD
} //End parallel region
#endif

  //Save to file
  OutputFile->cd();
  BatchedVariance->Write("BatchedMeansVariance");
  C_Test_Statistics->Write("C_Test_Statistics");

  //Delete all variables
  delete BatchedVariance;
  delete C_Test_Statistics;
}

// **************************
// RC: Perform spectral analysis of MCMC based on @cite Dunkley:2004sv
void MCMCProcessor::PowerSpectrumAnalysis() {
// **************************
  TStopwatch clock;
  clock.Start();

  //KS: Store it as we go back to them at the end
  const double TopMargin  = Posterior->GetTopMargin();
  const int OptTitle = gStyle->GetOptTitle();

  Posterior->SetTopMargin(0.1);
  gStyle->SetOptTitle(1);

  MACH3LOG_INFO("Making Power Spectrum plots...");

  // This is only to reduce number of computations...
  const int N_Coeffs = std::min(10000, nEntries);
  const int start = -(N_Coeffs/2-1);
  const int end = N_Coeffs/2-1;
  const int v_size = end - start;

  int nPrams = nDraw;
  /// @todo KS: Code is awfully slow... I know how to make it faster (GPU scream in a distant) but for now just make it for two params, bit hacky sry...
  nPrams = 1;

  std::vector<std::vector<float>> k_j(nPrams, std::vector<float>(v_size, 0.0));
  std::vector<std::vector<float>> P_j(nPrams, std::vector<float>(v_size, 0.0));

  int _N = nEntries;
  if (_N % 2 != 0) _N -= 1; // N must be even

  //This is being used a lot so calculate it once to increase performance
  const double two_pi_over_N = 2 * TMath::Pi() / static_cast<double>(_N);

  // KS: This could be moved to GPU I guess
  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2)
  #endif
  // RC: equation 11: for each value of j coef, from range -N/2 -> N/2
  for (int j = 0; j < nPrams; ++j)
  {
    for (int jj = start; jj < end; ++jj)
    {
      std::complex<double> a_j = 0.0;
      for (int n = 0; n < _N; ++n)
      {
        //if(StepNumber[n] < BurnInCut) continue;
        std::complex<double> exp_temp(0, two_pi_over_N * jj * n);
        a_j += ParStep[j][n] * std::exp(exp_temp);
      }
      a_j /= std::sqrt(float(_N));
      const int _c = jj - start;

      k_j[j][_c] = two_pi_over_N * jj;
      // Equation 13
      P_j[j][_c] = std::norm(a_j);
    }
  }

  TDirectory *PowerDir = OutputFile->mkdir("PowerSpectrum");
  PowerDir->cd();

  TVectorD* PowerSpectrumStepSize = new TVectorD(nPrams);
  for (int j = 0; j < nPrams; ++j)
  {
    auto plot = std::make_unique<TGraph>(v_size, k_j[j].data(), P_j[j].data());

    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(j, Prior, PriorError, Title);

    std::string name = Form("Power Spectrum of %s;k;P(k)", Title.Data());

    plot->SetTitle(name.c_str());
    name = Form("%s_power_spectrum", Title.Data());
    plot->SetName(name.c_str());
    plot->SetMarkerStyle(7);

    // Equation 18
    TF1 *func = new TF1("power_template", "[0]*( ([1] / x)^[2] / (([1] / x)^[2] +1) )", 0.0, 1.0);
    // P0 gives the amplitude of the white noise spectrum in the k → 0 limit
    func->SetParameter(0, 10.0);
    // k* indicates the position of the turnover to a different power law behaviour
    func->SetParameter(1, 0.1);
    // alpha free parameter
    func->SetParameter(2, 2.0);

    // Set parameter limits for stability
    func->SetParLimits(0, 0.0, 100.0); // Amplitude should be non-negative
    func->SetParLimits(1, 0.001, 1.0); // k* should be within a reasonable range
    func->SetParLimits(2, 0.0, 5.0);   // alpha should be positive

    plot->Fit("power_template","Rq");

    Posterior->SetLogx();
    Posterior->SetLogy();
    Posterior->SetGrid();
    plot->Write(plot->GetName());
    plot->Draw("AL");
    func->Draw("SAME");
    if(printToPDF) Posterior->Print(CanvasName);

    //KS: I have no clue what is the reason behind this. Found this in Rick Calland code...
    (*PowerSpectrumStepSize)(j) = std::sqrt(func->GetParameter(0)/float(v_size*0.5));
    delete func;
  }

  PowerSpectrumStepSize->Write("PowerSpectrumStepSize");
  delete PowerSpectrumStepSize;
  PowerDir->Close();
  delete PowerDir;

  clock.Stop();
  MACH3LOG_INFO("Making Power Spectrum took {:.2f}s", clock.RealTime());

  Posterior->SetTopMargin(TopMargin);
  gStyle->SetOptTitle(OptTitle);
}

// **************************
// Geweke Diagnostic based on
// @cite Fang2014GewekeDiagnostics
// @cite karlsbakk2011 Chapter 3.1
void MCMCProcessor::GewekeDiagnostic() {
// **************************
  MACH3LOG_INFO("Making Geweke Diagnostic");
  //KS: Up refers to upper limit we check, it stays constant, in literature it is mostly 50% thus using 0.5 for threshold
  std::vector<double> MeanUp(nDraw, 0.0);
  std::vector<double> SpectralVarianceUp(nDraw, 0.0);
  std::vector<int> DenomCounterUp(nDraw, 0);
  const double Threshold = 0.5 * nSteps;

  //KS: Select values between which you want to scan, for example 0 means 0% burn in and 1 100% burn in.
  constexpr double LowerThreshold = 0;
  constexpr double UpperThreshold = 1.0;
  // Tells how many intervals between thresholds we want to check
  constexpr int NChecks = 100;
  constexpr double Division = (UpperThreshold - LowerThreshold)/NChecks;

  std::vector<std::unique_ptr<TH1D>> GewekePlots(nDraw);
  for (int j = 0; j < nDraw; ++j)
  {
    TString Title = "";
    double Prior = 1.0, PriorError = 1.0;
    GetNthParameter(j, Prior, PriorError, Title);
    std::string HistName = Form("%s_%s_Geweke", Title.Data(), BranchNames[j].Data());
    GewekePlots[j] = std::make_unique<TH1D>(HistName.c_str(), HistName.c_str(), NChecks, 0.0, 100 * UpperThreshold);
    GewekePlots[j]->GetXaxis()->SetTitle("Burn-In (%)");
    GewekePlots[j]->GetYaxis()->SetTitle("Geweke T score");
  }

//KS: Start parallel region
#ifdef MULTITHREAD
#pragma omp parallel
{
#endif
  //KS: First we calculate mean and spectral variance for the upper limit, this doesn't change and in literature is most often 50%
  #ifdef MULTITHREAD
  #pragma omp for
  #endif
  for (int j = 0; j < nDraw; ++j)
  {
    for(int i = 0; i < nEntries; ++i)
    {
      if(StepNumber[i] > Threshold)
      {
        MeanUp[j] += ParStep[j][i];
        DenomCounterUp[j]++;
      }
    }
    MeanUp[j] = MeanUp[j]/DenomCounterUp[j];
  }

  //KS: now Spectral variance which in this case is sample variance
  #ifdef MULTITHREAD
  #pragma omp for collapse(2)
  #endif
  for (int j = 0; j < nDraw; ++j)
  {
    for(int i = 0; i < nEntries; ++i)
    {
      if(StepNumber[i] > Threshold)
      {
        SpectralVarianceUp[j] += (ParStep[j][i] - MeanUp[j])*(ParStep[j][i] - MeanUp[j]);
      }
    }
  }

  //Loop over how many intervals we calculate
  #ifdef MULTITHREAD
  #pragma omp for
  #endif
  for (int k = 1; k < NChecks+1; ++k)
  {
    //KS each thread has it's own
    std::vector<double> MeanDown(nDraw, 0.0);
    std::vector<double> SpectralVarianceDown(nDraw, 0.0);
    std::vector<int> DenomCounterDown(nDraw, 0);

    const int ThresholsCheck = Division*k*nSteps;
    //KS: First mean
    for (int j = 0; j < nDraw; ++j)
    {
      for(int i = 0; i < nEntries; ++i)
      {
        if(StepNumber[i] < ThresholsCheck)
        {
          MeanDown[j] += ParStep[j][i];
          DenomCounterDown[j]++;
        }
      }
      MeanDown[j] = MeanDown[j]/DenomCounterDown[j];
    }
    //Now spectral variance
    for (int j = 0; j < nDraw; ++j)
    {
      for(int i = 0; i < nEntries; ++i)
      {
        if(StepNumber[i] < ThresholsCheck)
        {
          SpectralVarianceDown[j] += (ParStep[j][i] - MeanDown[j])*(ParStep[j][i] - MeanDown[j]);
        }
      }
    }
    //Lastly calc T score and fill histogram entry
    for (int j = 0; j < nDraw; ++j)
    {
      double T_score = std::fabs((MeanDown[j] - MeanUp[j])/std::sqrt(SpectralVarianceDown[j]/DenomCounterDown[j] + SpectralVarianceUp[j]/DenomCounterUp[j]));
      GewekePlots[j]->SetBinContent(k, T_score);
    }
  } //end loop over intervals
#ifdef MULTITHREAD
} //End parallel region
#endif

  //Finally save it to TFile
  OutputFile->cd();
  TDirectory *GewekeDir = OutputFile->mkdir("Geweke");
  for (int j = 0; j < nDraw; ++j)
  {
    GewekeDir->cd();
    GewekePlots[j]->Write();
  }
  for (int i = 0; i < nDraw; ++i) {
    delete[] ParStep[i];
  }
  delete[] ParStep;

  GewekeDir->Close();
  delete GewekeDir;
  OutputFile->cd();
}

// **************************
// Acceptance Probability
void MCMCProcessor::AcceptanceProbabilities() {
// **************************
  if (AccProbBatchedAverages == nullptr) PrepareDiagMCMC();

  MACH3LOG_INFO("Making AccProb plots...");

  // Set the titles and limits for TH1Ds
  auto AcceptanceProbPlot = std::make_unique<TH1D>("AcceptanceProbability", "Acceptance Probability", nEntries, 0, nEntries);
  AcceptanceProbPlot->GetXaxis()->SetTitle("Step");
  AcceptanceProbPlot->GetYaxis()->SetTitle("Acceptance Probability");

  auto BatchedAcceptanceProblot = std::make_unique<TH1D>("AcceptanceProbability_Batch", "AcceptanceProbability_Batch", nBatches, 0, nBatches);
  BatchedAcceptanceProblot->GetYaxis()->SetTitle("Acceptance Probability");
  
  for (int i = 0; i < nBatches; ++i) {
    BatchedAcceptanceProblot->SetBinContent(i+1, AccProbBatchedAverages[i]);
    const int BatchRangeLow = double(i)*double(nEntries)/double(nBatches);
    const int BatchRangeHigh = double(i+1)*double(nEntries)/double(nBatches);
    std::stringstream ss;
    ss << BatchRangeLow << " - " << BatchRangeHigh;
    BatchedAcceptanceProblot->GetXaxis()->SetBinLabel(i+1, ss.str().c_str());
  }
  
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nEntries; ++i) {
    // Set bin content for the i-th bin to the parameter values
    AcceptanceProbPlot->SetBinContent(i, AccProbValues[i]);
  }
    
  TDirectory *probDir = OutputFile->mkdir("AccProb");
  probDir->cd();
  
  AcceptanceProbPlot->Write();
  BatchedAcceptanceProblot->Write();
  delete[] AccProbValues;
  delete[] AccProbBatchedAverages;

  probDir->Close();
  delete probDir;

  OutputFile->cd();
}

// **************************
void MCMCProcessor::CheckCredibleIntervalsOrder(const std::vector<double>& CredibleIntervals, const std::vector<Color_t>& CredibleIntervalsColours) {
// **************************
  if (CredibleIntervals.size() != CredibleIntervalsColours.size()) {
    MACH3LOG_ERROR("size of CredibleIntervals is not equal to size of CredibleIntervalsColours");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (CredibleIntervals.size() > 1) {
    for (unsigned int i = 1; i < CredibleIntervals.size(); i++) {
      if (CredibleIntervals[i] > CredibleIntervals[i - 1]) {
        MACH3LOG_ERROR("Interval {} is smaller than {}", i, i - 1);
        MACH3LOG_ERROR("{:.2f} {:.2f}", CredibleIntervals[i], CredibleIntervals[i - 1]);
        MACH3LOG_ERROR("They should be grouped in decreasing order");
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }
}

// **************************
void MCMCProcessor::CheckCredibleRegionsOrder(const std::vector<double>& CredibleRegions,
                               const std::vector<Style_t>& CredibleRegionStyle,
                               const std::vector<Color_t>& CredibleRegionColor) {
// **************************
  if ((CredibleRegions.size() != CredibleRegionStyle.size()) || (CredibleRegionStyle.size() != CredibleRegionColor.size())) {
    MACH3LOG_ERROR("size of CredibleRegions is not equal to size of CredibleRegionStyle or CredibleRegionColor");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for (unsigned int i = 1; i < CredibleRegions.size(); i++) {
    if (CredibleRegions[i] > CredibleRegions[i - 1]) {
      MACH3LOG_ERROR("Interval {} is smaller than {}", i, i - 1);
      MACH3LOG_ERROR("{:.2f} {:.2f}", CredibleRegions[i], CredibleRegions[i - 1]);
      MACH3LOG_ERROR("They should be grouped in decreasing order");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
}

// **************************
int MCMCProcessor::GetGroup(const std::string& name) const {
// **************************
  // Lambda to compare strings case-insensitively
  auto caseInsensitiveCompare = [](const std::string& a, const std::string& b) {
    return std::equal(a.begin(), a.end(), b.begin(), b.end(),
                      [](char c1, char c2) { return std::tolower(c1) == std::tolower(c2); });
  };
  int numerator = 0;
  for (const auto& groupName : ParameterGroup) {
    if (caseInsensitiveCompare(groupName, name)) {
      numerator++;
    }
  }
  return numerator;
}

// **************************
void MCMCProcessor::PrintInfo() const {
// **************************
  // KS: Create a map to store the counts of unique strings
  std::unordered_map<std::string, int> paramCounts;
  std::vector<std::string> orderedKeys;

  for (const std::string& param : ParameterGroup) {
    if (paramCounts[param] == 0) {
      orderedKeys.push_back(param);  // preserve order of first appearance
    }
    paramCounts[param]++;
  }

  MACH3LOG_INFO("************************************************");
  MACH3LOG_INFO("Scanning output branches...");
  MACH3LOG_INFO("# Useful entries in tree: \033[1;32m {} \033[0m ", nDraw);
  MACH3LOG_INFO("# Model params:  \033[1;32m {} starting at {} \033[0m ", nParam[kXSecPar], ParamTypeStartPos[kXSecPar]);
  MACH3LOG_INFO("# With following groups: ");
  for (const std::string& key : orderedKeys) {
    MACH3LOG_INFO(" # {} params: {}", key, paramCounts[key]);
  }
  MACH3LOG_INFO("# ND params (legacy):    \033[1;32m {} starting at {} \033[0m ", nParam[kNDPar], ParamTypeStartPos[kNDPar]);
  MACH3LOG_INFO("# FD params (legacy):    \033[1;32m {} starting at {} \033[0m ", nParam[kFDDetPar], ParamTypeStartPos[kFDDetPar]);
  MACH3LOG_INFO("************************************************");
}

// **************************
std::vector<double> MCMCProcessor::GetMargins(const std::unique_ptr<TCanvas>& Canv) const {
// **************************
  return std::vector<double>{Canv->GetTopMargin(), Canv->GetBottomMargin(),
                              Canv->GetLeftMargin(), Canv->GetRightMargin()};
}

// **************************
void MCMCProcessor::SetMargins(std::unique_ptr<TCanvas>& Canv, const std::vector<double>& margins) {
// **************************
  if (!Canv) {
    MACH3LOG_ERROR("Canv is nullptr");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (margins.size() != 4) {
    MACH3LOG_ERROR("Margin vector must have exactly 4 elements");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  Canv->SetTopMargin(margins[0]);
  Canv->SetBottomMargin(margins[1]);
  Canv->SetLeftMargin(margins[2]);
  Canv->SetRightMargin(margins[3]);
}

// **************************
void MCMCProcessor::SetTLineStyle(TLine* Line, const Color_t Colour, const Width_t Width, const ELineStyle Style) const {
// **************************
  Line->SetLineColor(Colour);
  Line->SetLineWidth(Width);
  Line->SetLineStyle(Style);
}

// **************************
void MCMCProcessor::SetLegendStyle(TLegend* Legend, const double size) const {
// **************************
  Legend->SetTextSize(size);
  Legend->SetLineColor(0);
  Legend->SetLineStyle(0);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetBorderSize(0);
}
