#include "MCMCProcessor.h"

//Only if GPU is enabled
#ifdef CUDA
extern void InitGPU_AutoCorr(
    float **ParStep_gpu,
    float **NumeratorSum_gpu,
    float **ParamSums_gpu,
    float **DenomSum_gpu,
    int n_Entries,
    int n_Pars,
    const int n_Lags);

extern void CopyToGPU_AutoCorr(
    float *ParStep_cpu,
    float *NumeratorSum_cpu,
    float *ParamSums_cpu,
    float *DenomSum_cpu,

    float *ParStep_gpu,
    float *NumeratorSum_gpu,
    float *ParamSums_gpu,
    float *DenomSum_gpu);

extern void RunGPU_AutoCorr(
    float *ParStep_gpu,
    float *ParamSums_gpu,
    float *NumeratorSum_gpu,
    float *DenomSum_gpu,
    float *NumeratorSum_cpu,
    float *DenomSum_cpu);

extern void CleanupGPU_AutoCorr(
    float *ParStep_gpu,
    float *NumeratorSum_gpu,
    float *ParamSums_gpu,
    float *DenomSum_gpu);
#endif

// ****************************
MCMCProcessor::MCMCProcessor(const std::string &InputFile, bool MakePostfitCorr) : 
  Chain(nullptr), StepCut(""), MakeCorr(MakePostfitCorr), MadePostfit(false) {
// ****************************
  MCMCFile = InputFile;

  SetMaCh3LoggerFormat();
  MaCh3Utils::MaCh3Welcome();
  MACH3LOG_INFO("Making post-fit processor for: {}", MCMCFile);

  ParStep = nullptr;
  StepNumber = nullptr;
    
  Posterior = nullptr;
  hpost = nullptr;
  hpost2D = nullptr;
  hviolin = nullptr;
  hviolin_prior = nullptr;

  OutputFile = nullptr;
  
  ParamSums = nullptr;
  BatchedAverages = nullptr;
  LagL = nullptr;
  SampleValues = nullptr;
  SystValues = nullptr;
  AccProbValues = nullptr;
  AccProbBatchedAverages = nullptr;
    
  //KS: Warning this only work when you project from Chain, will nor work when you try SetBranchAddress etc. Turn it on only if you know how to use it
  PlotJarlskog = false;
  
  //KS:Hardcoded should be a way to get it via config or something
  plotRelativeToPrior = false;
  printToPDF = false;
  plotBinValue = false;
  PlotFlatPrior = true;
  CacheMCMC = false;
  ApplySmoothing = true;
  FancyPlotNames = true;
  doDiagMCMC = false;
  OutputSuffix = "_Process";
  Post2DPlotThreshold = 1.e-5;

  nDraw = 0;
  nFlux = 0;
  nEntries = 0;
  UpperCut = _UNDEF_;
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
  
  for(int i = 0; i < kNParameterEnum; i++)
  {
    ParamTypeStartPos[i] = 0;
    nParam[i] = 0;
  }
  //Only if GPU is enabled
  #ifdef CUDA
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
  MACH3LOG_INFO("Closing pdf in MCMCProcessor: {}", CanvasName);
  CanvasName += "]";
  if(printToPDF) Posterior->Print(CanvasName);
  if (Posterior != nullptr) delete Posterior;

  delete Gauss;
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
  
  if(hpost != nullptr)
  {
    for (int i = 0; i < nDraw; ++i) 
    {
      delete hpost[i];
    }
    delete[] hpost;
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
      delete[] hpost2D[i];
    }
    delete[] ParStep;
    delete[] hpost2D;
  }
  if(StepNumber != nullptr) delete[] StepNumber;

  if(hviolin != nullptr) delete hviolin;
  if(hviolin_prior != nullptr) delete hviolin_prior;
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
  const int ParamTypeSize = ParamType.size();
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
  Cov = (TMatrixDSym*)Covariance->Clone();
  Corr = (TMatrixDSym*)Correlation->Clone();
}

// ***************
void MCMCProcessor::MakeOutputFile() {
// ***************

  //KS: ROOT hates me... but we can create several instances of MCMC Processor, each with own TCanvas ROOT is mad and will delete if there is more than one canvas with the same name, so we add random number to avoid issue
  TRandom3* rand = new TRandom3(0);
  const int uniform = int(rand->Uniform(0, 10000));
  // Open a TCanvas to write the posterior onto
  Posterior = new TCanvas(("Posterior" + std::to_string(uniform)).c_str(), ("Posterior" + std::to_string(uniform)).c_str(), 0, 0, 1024, 1024);

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

  // Directory for posteriors
  TDirectory *PostDir = OutputFile->mkdir("Post");

  // We fit with this Gaussian
  Gauss = new TF1("Gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",-5,5);
  Gauss->SetLineWidth(2);
  Gauss->SetLineColor(kOrange-5);
  
  // nDraw is number of draws we want to do
  for (int i = 0; i < nDraw; ++i)
  {
    if (i % (nDraw/5) == 0)
    {
      MaCh3Utils::PrintProgressBar(i, nDraw);
    }
    OutputFile->cd();
    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;
    
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

    GetArithmetic(hpost[i], i);
    GetGaussian(hpost[i], i);
    GetHPD(hpost[i], i);
    
    // Write the results from the projection into the TVectors and TMatrices
    (*Covariance)(i,i) = (*Errors)(i)*(*Errors)(i);
    (*Correlation)(i,i) = 1.0;

    //KS: This need to be before SetMaximum(), this way plot is nicer as line end at the maximum
    TLine *hpd = new TLine((*Means_HPD)(i), hpost[i]->GetMinimum(), (*Means_HPD)(i), hpost[i]->GetMaximum());
    hpd->SetLineColor(kBlack);
    hpd->SetLineWidth(2);
    hpd->SetLineStyle(kSolid);
    
    hpost[i]->SetLineWidth(2);
    hpost[i]->SetLineColor(kBlue-1);
    hpost[i]->SetMaximum(hpost[i]->GetMaximum()*DrawRange);
    hpost[i]->SetTitle(Title);
    hpost[i]->GetXaxis()->SetTitle(hpost[i]->GetTitle());
    
    // Now make the TLine for the Asimov
    TLine *Asimov = new TLine(Prior, hpost[i]->GetMinimum(), Prior, hpost[i]->GetMaximum());
    Asimov->SetLineColor(kRed-3);
    Asimov->SetLineWidth(2);
    Asimov->SetLineStyle(kDashed);

    TLegend *leg = new TLegend(0.12, 0.6, 0.6, 0.97);
    leg->SetTextSize(0.04);
    leg->AddEntry(hpost[i], Form("#splitline{PDF}{#mu = %.2f, #sigma = %.2f}", hpost[i]->GetMean(), hpost[i]->GetRMS()), "l");
    leg->AddEntry(Gauss, Form("#splitline{Gauss}{#mu = %.2f, #sigma = %.2f}", Gauss->GetParameter(1), Gauss->GetParameter(2)), "l");
    leg->AddEntry(hpd, Form("#splitline{HPD}{#mu = %.2f, #sigma = %.2f (+%.2f-%.2f)}", (*Means_HPD)(i), (*Errors_HPD)(i), (*Errors_HPD_Positive)(i), (*Errors_HPD_Negative)(i)), "l");
    leg->AddEntry(Asimov, Form("#splitline{Prior}{x = %.2f , #sigma = %.2f}", Prior, PriorError), "l");
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);

    //CW: Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
    if (hpost[i]->GetMaximum() == hpost[i]->Integral()*DrawRange) 
    {
      MACH3LOG_WARN("Found fixed parameter, moving on");
      IamVaried[i] = false;
      //KS:Set mean and error to prior for fixed parameters, it looks much better when fixed parameter has mean on prior rather than on 0 with 0 error.
      (*Means_HPD)(i)  = Prior;
      (*Errors_HPD)(i) = PriorError;

      delete Asimov;
      delete hpd;
      delete leg;
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

    delete Asimov;
    delete hpd;
    delete leg;
  } // end the for loop over nDraw

  OutputFile->cd();
  TTree *SettingsBranch = new TTree("Settings", "Settings");
  int CrossSectionParameters = nParam[kXSecPar];
  SettingsBranch->Branch("CrossSectionParameters", &CrossSectionParameters);
  int CrossSectionParametersStartingPos = ParamTypeStartPos[kXSecPar];
  SettingsBranch->Branch("CrossSectionParametersStartingPos", &CrossSectionParametersStartingPos);
  int FluxParameters = nFlux;
  SettingsBranch->Branch("FluxParameters", &FluxParameters);
  
  int NDParameters = nParam[kNDPar];
  SettingsBranch->Branch("NDParameters", &NDParameters);
  int NDParametersStartingPos = ParamTypeStartPos[kNDPar];
  SettingsBranch->Branch("NDParametersStartingPos", &NDParametersStartingPos);

  int FDParameters = nParam[kFDDetPar];
  SettingsBranch->Branch("FDParameters", &FDParameters);
  int FDParametersStartingPos = ParamTypeStartPos[kFDDetPar];
  SettingsBranch->Branch("FDParametersStartingPos", &FDParametersStartingPos);
  
  int OscParameters = nParam[kOSCPar];
  SettingsBranch->Branch("OscParameters", &OscParameters);
  int OscParametersStartingPos = ParamTypeStartPos[kOSCPar];
  SettingsBranch->Branch("OscParametersStartingPos", &OscParametersStartingPos);
  
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
} // Have now written the postfit projections

// *******************
//CW: Draw the postfit
void MCMCProcessor::DrawPostfit() {
// *******************

  if (OutputFile == nullptr) MakeOutputFile();

  // Make the prefit plot
  TH1D* prefit = MakePrefit();

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

  // Same but with Gaussian output
  TH1D *paramPlot_Gauss = (TH1D*)(paramPlot->Clone());
  paramPlot_Gauss->SetMarkerColor(kOrange-5);
  paramPlot_Gauss->SetMarkerStyle(23);
  paramPlot_Gauss->SetLineWidth(2);
  paramPlot_Gauss->SetMarkerSize((prefit->GetMarkerSize())*0.75);
  paramPlot_Gauss->SetFillColor(paramPlot_Gauss->GetMarkerColor());
  paramPlot_Gauss->SetFillStyle(3244);
  paramPlot_Gauss->SetLineColor(paramPlot_Gauss->GetMarkerColor());

  // Same but with Gaussian output
  TH1D *paramPlot_HPD = (TH1D*)(paramPlot->Clone());
  paramPlot_HPD->SetMarkerColor(kBlack);
  paramPlot_HPD->SetMarkerStyle(25);
  paramPlot_HPD->SetLineWidth(2);
  paramPlot_HPD->SetMarkerSize((prefit->GetMarkerSize())*0.5);
  paramPlot_HPD->SetFillColor(0);
  paramPlot_HPD->SetFillStyle(0);
  paramPlot_HPD->SetLineColor(paramPlot_HPD->GetMarkerColor());
  
  // Set labels and data
  for (int i = 0; i < nDraw; ++i)
  {
    //Those keep which parameter type we run currently and realtive number  
    int ParamEnu = ParamType[i];
    int ParamNo = i - ParamTypeStartPos[ParameterEnum(ParamEnu)];

    //KS: Sliglthy hacky way to get realtive to prior or nominal as this is convention we use
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

  // Make a TLegend
  TLegend *CompLeg = new TLegend(0.33, 0.73, 0.76, 0.95);
  CompLeg->AddEntry(prefit, "Prefit", "fp");
  CompLeg->AddEntry(paramPlot, "Postfit PDF", "fp");
  CompLeg->AddEntry(paramPlot_Gauss, "Postfit Gauss", "fp");
  CompLeg->AddEntry(paramPlot_HPD, "Postfit HPD", "lfep");
  CompLeg->SetFillColor(0);
  CompLeg->SetFillStyle(0);
  CompLeg->SetLineWidth(0);
  CompLeg->SetLineStyle(0);
  CompLeg->SetBorderSize(0);

  const double BottomMargin = Posterior->GetBottomMargin();
  Posterior->SetBottomMargin(0.2);

  OutputFile->cd();
  //KS: Plot Xsec and Flux
  if (PlotXSec == true)
  {
    const int Start = ParamTypeStartPos[kXSecPar];
    // Plot the xsec parameters (0 to ~nXsec-nFlux) nXsec == xsec + flux, quite confusing I know
    // Have already looked through the branches earlier
    if(plotRelativeToPrior)  prefit->GetYaxis()->SetTitle("Variation rel. prior"); 
    else prefit->GetYaxis()->SetTitle("Parameter Value");
    prefit->GetYaxis()->SetRangeUser(-2.5, 2.5);
    prefit->GetXaxis()->SetTitle("");
    prefit->GetXaxis()->LabelsOption("v");

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
  if(PlotDet)
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
      prefit->GetXaxis()->SetTitle();
      prefit->GetXaxis()->LabelsOption("v");

      paramPlot->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      paramPlot->GetYaxis()->SetRangeUser(0.6, 1.4);
      paramPlot->GetXaxis()->SetRangeUser(Start, NDbinCounter);
      paramPlot->GetXaxis()->SetTitle("");
      paramPlot->SetTitle(CutPosterior1D.c_str());
      paramPlot->GetXaxis()->LabelsOption("v");

      paramPlot_Gauss->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      paramPlot_Gauss->GetYaxis()->SetRangeUser(0.6, 1.4);
      paramPlot_Gauss->GetXaxis()->SetRangeUser(Start, NDbinCounter);
      paramPlot_Gauss->GetXaxis()->SetTitle("");
      paramPlot_Gauss->SetTitle(CutPosterior1D.c_str());
      paramPlot_Gauss->GetXaxis()->LabelsOption("v");

      paramPlot_HPD->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      paramPlot_HPD->GetYaxis()->SetRangeUser(0.6, 1.4);
      paramPlot_HPD->GetXaxis()->SetRangeUser(Start, NDbinCounter);
      paramPlot_HPD->GetXaxis()->SetTitle("");
      paramPlot_HPD->SetTitle(CutPosterior1D.c_str());
      paramPlot_HPD->GetXaxis()->LabelsOption("v");

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

  delete prefit;
  delete paramPlot;
  delete paramPlot_Gauss;
  delete paramPlot_HPD;
  delete CompLeg;

  //KS: Return Margin to default one
  Posterior->SetBottomMargin(BottomMargin);
}

// *********************
// Make fancy Credible Intervals plots
void MCMCProcessor::MakeCredibleIntervals(std::vector<double> CredibleIntervals,
                                          std::vector<Color_t> CredibleIntervalsColours,
                                          bool CredibleInSigmas) {
// *********************

  if(hpost[0] == nullptr) MakePostfit();

  MACH3LOG_INFO("Making Credible Intervals ");

  const double LeftMargin = Posterior->GetLeftMargin();
  Posterior->SetLeftMargin(0.15);

  if(CredibleIntervals.size() != CredibleIntervalsColours.size())
  {
    MACH3LOG_ERROR("Size of  CredibleIntervals is not equal to size of CredibleIntervalsColours");
    throw;
  }

  if(CredibleIntervals.size() > 1)
  {
    for(unsigned int i = 1; i < CredibleIntervals.size(); i++ )
    {
      if(CredibleIntervals[i] > CredibleIntervals[i-1])
      {
        std::cerr<<" Interval "<<i<<" is smaller than "<<i-1<<std::endl;
        std::cerr<<CredibleIntervals[i] <<" "<<CredibleIntervals[i-1]<<std::endl;
        std::cerr<<" They should be grouped in decreasing order"<<std::endl;
        std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    }
  }

  const int nCredible = CredibleIntervals.size();
  TH1D** hpost_copy = new TH1D*[nDraw];
  TH1D*** hpost_cl = new TH1D**[nDraw];

  //KS: Copy all histograms to be thread safe
  for (int i = 0; i < nDraw; ++i)
  {
    hpost_copy[i] = (TH1D*) hpost[i]->Clone(Form("hpost_copy_%i", i));
    hpost_cl[i] = new TH1D*[nCredible];

    for (int j = 0; j < nCredible; ++j)
    {
      hpost_cl[i][j] = (TH1D*) hpost[i]->Clone( Form("hpost_copy_%i_CL_%f", i, CredibleIntervals[j]));

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

      //KS: We have slightly different approach depending if you passed percentage or sigmas
      if(CredibleInSigmas)
      {
        //KS: Convert sigmas into percentage
        const double CredInter = GetSigmaValue((int)std::round(CredibleIntervals[j]));
        GetCredibleInterval(hpost_copy[i], hpost_cl[i][j], CredInter);
      }
      else
      {
        GetCredibleInterval(hpost_copy[i], hpost_cl[i][j], CredibleIntervals[j]);
      }

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
    double Prior = 1.0;
    double PriorError = 1.0;

    GetNthParameter(i, Prior, PriorError, Title);

    TLine *Asimov = new TLine(Prior, hpost_copy[i]->GetMinimum(), Prior, hpost_copy[i]->GetMaximum());
    Asimov->SetLineColor(kRed-3);
    Asimov->SetLineWidth(2);
    Asimov->SetLineStyle(kDashed);

    TLegend* legend = new TLegend(0.20, 0.7, 0.4, 0.92);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetLineColor(0);
    legend->SetLineStyle(0);
    legend->SetBorderSize(0);
    hpost_copy[i]->Draw("HIST");

    for (int j = 0; j < nCredible; ++j)
        hpost_cl[i][j]->Draw("HIST SAME");
    for (int j = nCredible-1; j >= 0; --j)
    {
      if(CredibleInSigmas)
        legend->AddEntry(hpost_cl[i][j], Form("%.0f#sigma Credible Interval", CredibleIntervals[j]), "f");
      else
        legend->AddEntry(hpost_cl[i][j], Form("%.0f%% Credible Interval", CredibleIntervals[j]*100), "f");
    }
    legend->AddEntry(Asimov, Form("#splitline{Prior}{x = %.2f , #sigma = %.2f}", Prior, PriorError), "l");
    legend->Draw("SAME");
    Asimov->Draw("SAME");

    // Write to file
    Posterior->SetName(hpost[i]->GetName());
    Posterior->SetTitle(hpost[i]->GetTitle());

    if(printToPDF) Posterior->Print(CanvasName);
    // cd into directory in root file
    CredibleDir->cd();
    Posterior->Write();

    delete legend;
    delete Asimov;
  }

  //KS: Remove histograms
  for (int i = 0; i < nDraw; ++i)
  {
    delete hpost_copy[i];
    for (int j = 0; j < nCredible; ++j)
    {
      delete hpost_cl[i][j];
    }
    delete[] hpost_cl[i];
  }
  delete[] hpost_copy;
  delete[] hpost_cl;

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
    if(Chain->GetMaximum(BranchNames[i]) > maxi_y) maxi_y = Chain->GetMaximum(BranchNames[i]);
    if(Chain->GetMinimum(BranchNames[i]) < mini_y) mini_y = Chain->GetMinimum(BranchNames[i]);
  }

  const int vBins = (maxi_y-mini_y)*25;

  hviolin = new TH2D("hviolin", "hviolin", nDraw, 0, nDraw, vBins, mini_y, maxi_y);

  //KS: Prior has larger errors so we increase range and number of bins
  const int PriorFactor = 4;
  hviolin_prior = new TH2D("hviolin_prior", "hviolin_prior", nDraw, 0, nDraw, PriorFactor*vBins, PriorFactor*mini_y, PriorFactor*maxi_y);

  TRandom3* rand = new TRandom3(0);

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
  const int IntervalsSize = 10;
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
    hviolin_prior->Draw("VIOLIN");
    hviolin->Draw("VIOLIN SAME");
    if(printToPDF) Posterior->Print(CanvasName);
  }
  delete rand;
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
    if ((*Covariance)(i,i) == _UNDEF_) {
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

  int covBinning = nDraw;
  // Now we are sure we have the diagonal elements, let's make the off-diagonals
  for (int i = 0; i < covBinning; ++i) 
  {
    if (i % (covBinning/5) == 0) 
      MaCh3Utils::PrintProgressBar(i, covBinning);

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
      TH2D *hpost_2D = new TH2D(DrawMe, DrawMe, nBins, min_i, max_i, nBins, min_j, max_j);

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
          if(IsXsec[j] && IsXsec[i] && std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold)
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

      delete hpost_2D;
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
    throw;
  }

  MACH3LOG_INFO("Caching input tree...");
  MACH3LOG_INFO("Allocating {:.2f} MB", (sizeof(double)*nDraw*nEntries)/1.E6);
  TStopwatch clock;
  clock.Start();
  
  ParStep = new double*[nDraw];
  StepNumber = new int[nEntries];
  
  hpost2D = new TH2D**[nDraw]();

  for (int i = 0; i < nDraw; ++i) 
  {
    ParStep[i] = new double[nEntries];
    hpost2D[i] = new TH2D*[nDraw]();

    for (int j = 0; j < nEntries; ++j)
    {
      ParStep[i][j] = -999.99;
      //KS: Set this only once
      if(i == 0) StepNumber[j] = -999.99;
    }
  }

  // Set all the branches to off
  Chain->SetBranchStatus("*", false);
  
  // Turn on the branches which we want for parameters
  for (int i = 0; i < nDraw; ++i) 
  {
    Chain->SetBranchStatus(BranchNames[i].Data(), true);
  }
  Chain->SetBranchStatus("step", true);

  const int countwidth = nEntries/10;
  // Loop over the entries
  //KS: This is really a bottleneck right now, thus revisit with ROOT6 https://pep-root6.github.io/docs/analysis/parallell/root.html
  for (int j = 0; j < nEntries; ++j) 
  {
    if (j % countwidth == 0)
        MaCh3Utils::PrintProgressBar(j, nEntries);

    Chain->SetBranchAddress("step", &StepNumber[j]);
    // Set the branch addresses for params
    for (int i = 0; i < nDraw; ++i) 
    {
      Chain->SetBranchAddress(BranchNames[i].Data(), &ParStep[i][j]);
    }
    
    if (j % countwidth == 0) {
      MaCh3Utils::EstimateDataTransferRate(Chain, j);
    } else {
      // Fill up the ParStep array
      Chain->GetEntry(j);
    }

  }
  
  // Set all the branches to on
  Chain->SetBranchStatus("*", true);
  
  // Cache max and min in chain for covariance matrix
  for (int i = 0; i < nDraw; ++i) 
  {
    const double Min_Chain_i = Chain->GetMinimum(BranchNames[i]);
    const double Max_Chain_i = Chain->GetMaximum(BranchNames[i]);
    
    TString Title_i = "";
    double Prior_i, PriorError_i;
    GetNthParameter(i, Prior_i, PriorError_i, Title_i);
    
    for (int j = 0; j <= i; ++j)
    {
      const double Min_Chain_j = Chain->GetMinimum(BranchNames[j]);
      const double Max_Chain_j = Chain->GetMaximum(BranchNames[j]);
  
      // TH2D to hold the Correlation 
      hpost2D[i][j] = new TH2D(Form("hpost2D_%i_%i",i,j), Form("hpost2D_%i_%i",i,j), nBins, Min_Chain_i, Max_Chain_i, nBins, Min_Chain_j, Max_Chain_j);
      
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
  
  int covBinning = nDraw;

  bool HaveMadeDiagonal = false;    
  // Check that the diagonal entries have been filled
  // i.e. MakePostfit() has been called
  for (int i = 0; i < covBinning; ++i) {
    if ((*Covariance)(i,i) == _UNDEF_) {
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
  for (int i = 0; i < covBinning; ++i) 
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

  if(!Mute) clock.Stop();
  if(!Mute) MACH3LOG_INFO("Making Covariance took {:.2f}s to finish for {} steps", clock.RealTime(), nEntries);

  OutputFile->cd();
  if(printToPDF)
  {
    Posterior->cd();
    for (int i = 0; i < covBinning; ++i) 
    {    
      for (int j = 0; j <= i; ++j)
      {
        // Skip the diagonal elements which we've already done above
        if (j == i) continue;
        if (IamVaried[j] == false) continue;

        if(ParamType[i] == kXSecPar && ParamType[j] == kXSecPar)
        {
          //KS: Skip Flux Params
          if(IsXsec[j] && IsXsec[i] && std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold)
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
  if(!Mute) Covariance->Write("Covariance");
  if(!Mute) Correlation->Write("Correlation");
}


// *********************
// Based on https://www.jstor.org/stable/25651249?seq=3,
// all credits for finding and studying it goes to Henry
void MCMCProcessor::MakeSubOptimality(int NIntervals) {
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

  TH1D* SubOptimality = new TH1D("Suboptimality", "Suboptimality", NIntervals, MinStep, MaxStep);
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

  delete SubOptimality;
}

// *********************
// Make the covariance plots
void MCMCProcessor::DrawCovariance() {
// *********************
    
  int covBinning = nDraw;
  // The Covariance matrix from the fit
  TH2D* hCov = new TH2D("hCov", "hCov", covBinning, 0, covBinning, covBinning, 0, covBinning);
  hCov->GetZaxis()->SetTitle("Covariance");
  // The Covariance matrix square root, with correct sign
  TH2D* hCovSq = new TH2D("hCovSq", "hCovSq", covBinning, 0, covBinning, covBinning, 0, covBinning);
  hCovSq->GetZaxis()->SetTitle("Covariance");
  // The Correlation
  TH2D* hCorr = new TH2D("hCorr", "hCorr", covBinning, 0, covBinning, covBinning, 0, covBinning);
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
  for (int i = 0; i < covBinning; ++i)
  {
    TString titlex = "";
    double nom, err;
    GetNthParameter(i, nom, err, titlex);
    
    hCov->GetXaxis()->SetBinLabel(i+1, titlex);
    hCovSq->GetXaxis()->SetBinLabel(i+1, titlex);
    hCorr->GetXaxis()->SetBinLabel(i+1, titlex);

    for (int j = 0; j < covBinning; ++j) 
    {
      // The value of the Covariance
      double cov = (*Covariance)(i,j);
      double corr = (*Correlation)(i,j);

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
  Posterior->SetRightMargin(0.15);
  if(printToPDF) Posterior->Print(CanvasName);

  Posterior->cd();
  Posterior->Clear();
  if(plotBinValue) hCorr->Draw("colz text");
  else hCorr->Draw("colz");
  Posterior->SetRightMargin(0.15);
  if(printToPDF) Posterior->Print(CanvasName);

  hCov->Write("Covariance_plot");
  hCovSq->Write("Covariance_sq_plot");
  hCorr->Write("Correlation_plot");
  
  //Back to normal
  Posterior->SetRightMargin(0.03);
  delete hCov;
  delete hCovSq;
  delete hCorr;

  DrawCorrelations1D();
}


// *********************
//KS: Make the 1D projections of Correlations inspired by Henry's slides (page 28) https://www.t2k.org/asg/oagroup/meeting/2023/2023-07-10-oa-pre-meeting/MaCh3FDUpdate
void MCMCProcessor::DrawCorrelations1D() {
// *********************

  //KS: Store it as we go back to them at the end
  const double TopMargin  = Posterior->GetTopMargin();
  const double BottomMargin  = Posterior->GetBottomMargin();
  const int OptTitle = gStyle->GetOptTitle();

  Posterior->SetTopMargin(0.1);
  Posterior->SetBottomMargin(0.2);
  gStyle->SetOptTitle(1);

  const int Nhists = 3;
  //KS: Highest value is just meant bo be sliglhy higher than 1 to catch >,
  const double Thresholds[Nhists+1] = {0, 0.25, 0.5, 1.0001};
  const Color_t CorrColours[Nhists] = {kRed-10, kRed-6,  kRed};

  //KS: This strore neccesary entires for stripped covariance which strore only "menaingfull correlations
  std::vector<std::vector<double>> CorrOfInterest;
  CorrOfInterest.resize(nDraw);
  std::vector<std::vector<std::string>> NameCorrOfInterest;
  NameCorrOfInterest.resize(nDraw);

  TH1D ***Corr1DHist = new TH1D**[nDraw]();
  //KS: Initialising ROOT objects is never safe in MP loop
  for(int i = 0; i < nDraw; ++i)
  {
    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;
    GetNthParameter(i, Prior, PriorError, Title);

    Corr1DHist[i] = new TH1D*[Nhists]();
    for(int j = 0; j < Nhists; ++j)
    {
      Corr1DHist[i][j] = new TH1D(Form("Corr1DHist_%i_%i", i, j), Form("Corr1DHist_%i_%i", i, j), nDraw, 0, nDraw);
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
    for(int k = 1; k < Nhists; k++)
    {
      Corr1DHist[i][k]->Draw("SAME");
    }

    TLegend *leg = new TLegend(0.3, 0.75, 0.6, 0.90);
    leg->SetTextSize(0.02);
    for(int k = 0; k < Nhists; k++)
    {
      leg->AddEntry(Corr1DHist[i][k], Form("%.2f > |Corr| >= %.2f", Thresholds[k+1], Thresholds[k]), "f");
    }
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->Draw("SAME");

    Posterior->Write(Corr1DHist[i][0]->GetTitle());
    if(printToPDF) Posterior->Print(CanvasName);

    delete leg;
  }

  //KS: Plot only meaninfull correlations
  for(int i = 0; i < nDraw; i++)
  {
    const int size = CorrOfInterest[i].size();

    if(size == 0) continue;
    TH1D* Corr1DHist_Reduced = new TH1D("Corr1DHist_Reduced", "Corr1DHist_Reduced", size, 0, size);

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

    delete Corr1DHist_Reduced;
  }

  for(int i = 0; i < nDraw; i++)
  {
    for(int k = 1; k < Nhists; k++)
    {
      delete Corr1DHist[i][k];
    }
    delete[] Corr1DHist[i];
  }
  delete[] Corr1DHist;

  CorrDir->Close();
  delete CorrDir;

  OutputFile->cd();

  Posterior->SetTopMargin(TopMargin);
  Posterior->SetBottomMargin(BottomMargin);
  gStyle->SetOptTitle(OptTitle);

}

// *********************
// Make fancy Credible Intervals plots
void MCMCProcessor::MakeCredibleRegions(std::vector<double> CredibleRegions,
                                        std::vector<Style_t> CredibleRegionStyle,
                                        std::vector<Color_t> CredibleRegionColor,
                                        bool CredibleInSigmas) {
// *********************

  if(hpost2D == nullptr) MakeCovariance_MP();
  MACH3LOG_INFO("Making Credible Regions");

  if( (CredibleRegions.size() != CredibleRegionStyle.size()) || (CredibleRegionStyle.size() != CredibleRegionColor.size()) )
  {
    std::cerr<<" size of  CredibleRegions is not equat to size of CredibleRegionStyle"<<std::endl;
    std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  const int nCredible = CredibleRegions.size();
  TH2D*** hpost_2D_copy = new TH2D**[nDraw];
  TH2D**** hpost_2D_cl = new TH2D***[nDraw];

  //KS: Copy all histograms to be thread safe
  for (int i = 0; i < nDraw; ++i)
  {
    hpost_2D_copy[i] = new TH2D*[nDraw];
    hpost_2D_cl[i] = new TH2D**[nDraw];
    for (int j = 0; j <= i; ++j)
    {
      hpost_2D_copy[i][j] = (TH2D*) hpost2D[i][j]->Clone( Form("hpost_copy_%i_%i", i, j));

      hpost_2D_cl[i][j] = new TH2D*[nCredible];
      for (int k = 0; k < nCredible; ++k)
      {
        hpost_2D_cl[i][j][k] = (TH2D*)hpost2D[i][j]->Clone( Form("hpost_copy_%i_%i_CL_%f", i, j, CredibleRegions[k]));;
      }
    }
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  //Calcualte creadible histogram
  for (int i = 0; i < nDraw; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      for (int k = 0; k < nCredible; ++k)
      {
        if(CredibleInSigmas)
        {
          //KS: Convert sigmas into percentage
          double CredReg = GetSigmaValue((int)std::round(CredibleRegions[k]));
          GetCredibleRegion(hpost_2D_cl[i][j][k], CredReg);
        }
        else
        {
          GetCredibleRegion(hpost_2D_cl[i][j][k], CredibleRegions[k]);
        }
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

      TLegend* legend = new TLegend(0.20, 0.7, 0.4, 0.92);
      legend->SetTextColor(kRed);
      legend->SetTextSize(0.03);
      legend->SetFillColor(0);
      legend->SetFillStyle(0);
      legend->SetLineColor(0);
      legend->SetLineStyle(0);
      legend->SetBorderSize(0);

      //Get Best point
      TGraph *bestfitM = new TGraph(1);
      const int MaxBin = hpost_2D_copy[i][j]->GetMaximumBin();
      int Mbx, Mby, Mbz;
      hpost_2D_copy[i][j]->GetBinXYZ(MaxBin, Mbx, Mby, Mbz);
      const double Mx = hpost_2D_copy[i][j]->GetXaxis()->GetBinCenter(Mbx);
      const double My = hpost_2D_copy[i][j]->GetYaxis()->GetBinCenter(Mby);

      bestfitM->SetPoint(0, Mx, My);
      bestfitM->SetMarkerStyle(22);
      bestfitM->SetMarkerSize(1);
      bestfitM->SetMarkerColor(kMagenta);
      legend->AddEntry(bestfitM,"Best Fit","p");

      //Plot default 2D posterior
      hpost_2D_copy[i][j]->Draw("COLZ");

      //Now credible regions
      for (int k = 0; k < nCredible; ++k)
        hpost_2D_cl[i][j][k]->Draw("CONT3 SAME");
      for (int k = nCredible-1; k >= 0; --k)
      {
        if(CredibleInSigmas)
          legend->AddEntry(hpost_2D_cl[i][j][k], Form("%.0f#sigma Credible Interval", CredibleRegions[k]), "l");
        else
          legend->AddEntry(hpost_2D_cl[i][j][k], Form("%.0f%% Credible Region", CredibleRegions[k]*100), "l");
      }

      legend->Draw("SAME");
      bestfitM->Draw("SAME.P");

      // Write to file
      Posterior->SetName(hpost2D[i][j]->GetName());
      Posterior->SetTitle(hpost2D[i][j]->GetTitle());

      //KS: Print only regions with correlation greater than specified value, by defualt 0.2. This is done to avoid dumping thousands of plots
      if(printToPDF && std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold) Posterior->Print(CanvasName);
      // Write it to root file
      //OutputFile->cd();
      //if( std::fabs((*Correlation)(i,j)) > Post2DPlotThreshold ) Posterior->Write();

      delete legend;
      delete bestfitM;
    }
  }

  OutputFile->cd();
  //KS: Remove histograms
  for (int i = 0; i < nDraw; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      delete hpost_2D_copy[i][j];
      for (int k = 0; k < nCredible; ++k)
      {
        delete hpost_2D_cl[i][j][k];
      }
      delete[] hpost_2D_cl[i][j];
    }
    delete[] hpost_2D_copy[i];
    delete[] hpost_2D_cl[i];
  }
  delete[] hpost_2D_copy;
  delete[] hpost_2D_cl;
}


// *********************
// Make fancy triangle plot for selected parameters
void MCMCProcessor::MakeTrianglePlot(std::vector<std::string> ParNames,
                                     // 1D
                                     std::vector<double> CredibleIntervals,
                                     std::vector<Color_t> CredibleIntervalsColours,
                                     //2D
                                     std::vector<double> CredibleRegions,
                                     std::vector<Style_t> CredibleRegionStyle,
                                     std::vector<Color_t> CredibleRegionColor,
                                     // Other
                                     bool CredibleInSigmas) {
// *********************

  if(hpost2D == nullptr) MakeCovariance_MP();
  MACH3LOG_INFO("Making Triangle Plot");

  const int nParamPlot = ParNames.size();
  std::vector<int> ParamNumber;
  for(int j = 0; j < nParamPlot; ++j)
  {
    //KS: First we need to find parameter number based on name
    int ParamNo = _UNDEF_;
    for (int i = 0; i < nDraw; ++i)
    {
      TString Title = "";
      double Prior = 1.0;
      double PriorError = 1.0;

      GetNthParameter(i, Prior, PriorError, Title);

      if(ParNames[j] == Title) ParamNo = i;
    }
    if(ParamNo == _UNDEF_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Will not plot Triangle plot", ParNames[j]);
      return;
    }
    ParamNumber.push_back(ParamNo);
  }

  //KS: Store it as we go back to them at the end
  const double TopMargin    = Posterior->GetTopMargin();
  const double BottomMargin = Posterior->GetBottomMargin();
  const double LeftMargin   = Posterior->GetLeftMargin();
  const double RighMargin   = Posterior->GetRightMargin();
  Posterior->SetTopMargin(0.001);
  Posterior->SetBottomMargin(0.001);
  Posterior->SetLeftMargin(0.001);
  Posterior->SetRightMargin(0.001);

  Posterior->cd();
  Posterior->Clear();
  Posterior->Update();

  //KS: We sort to have parameters from highest to lowest, this is related to how we make 2D projections in MakeCovariance_MP
  std::sort(ParamNumber.begin(), ParamNumber.end(),  std::greater<int>());

  //KS: Calculate how many pads/plots we need
  int Npad = 0;
  for(int j = 1; j < nParamPlot+1; j++) Npad += j;
  Posterior->cd();

  if(CredibleIntervals.size() != CredibleIntervalsColours.size())
  {
    std::cerr<<" size of  CredibleIntervals is not equat to size of CredibleIntervalsColours"<<std::endl;
    std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  if(CredibleIntervals.size() > 1)
  {
    for(unsigned int i = 1; i < CredibleIntervals.size(); i++ )
    {
      if(CredibleIntervals[i] > CredibleIntervals[i-1])
      {
        std::cerr<<" Interval "<<i<<" is smaller than "<<i-1<<std::endl;
        std::cerr<<CredibleIntervals[i] <<" "<<CredibleIntervals[i-1]<<std::endl;
        std::cerr<<" They should be grouped in decreasing order"<<std::endl;
        std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    }
  }
  if( (CredibleRegions.size() != CredibleRegionStyle.size()) || (CredibleRegionStyle.size() != CredibleRegionColor.size()) )
  {
    std::cerr<<" size of  CredibleRegions is not equat to size of CredibleRegionStyle"<<std::endl;
    std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  for(unsigned int i = 1; i < CredibleRegions.size(); i++ )
  {
    if(CredibleRegions[i] > CredibleRegions[i-1])
    {
      std::cerr<<" Interval "<<i<<" is smaller than "<<i-1<<std::endl;
      std::cerr<<CredibleRegions[i] <<" "<<CredibleRegions[i-1]<<std::endl;
      std::cerr<<" They should be grouped in decreasing order"<<std::endl;
      std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
  }
  const int nCredibleIntervals = CredibleIntervals.size();
  const int nCredibleRegions = CredibleRegions.size();

  //KS: Initialise Tpad histograms etc we will need
  TPad** TrianglePad = new TPad*[Npad];
  //KS: 1D copy of posterior, we need it as we modify them
  TH1D** hpost_copy = new TH1D*[nParamPlot];
  TH1D*** hpost_cl = new TH1D**[nParamPlot];
  TText **TriangleText = new TText *[nParamPlot*2];
  TH2D** hpost_2D_copy = new TH2D*[Npad-nParamPlot];
  TH2D*** hpost_2D_cl = new TH2D**[Npad-nParamPlot];
  gStyle->SetPalette(51);

  //KS: Super convoluted way of calculating ranges for our pads, trust me it works...
  double* X_Min = new double[nParamPlot];
  double* X_Max = new double[nParamPlot];

  X_Min[0] = 0.10;
  double xScale = (0.95 - (X_Min[0]+0.05))/nParamPlot;
  //KS: 0.05 is becasue we need additional offset for labels
  X_Max[0] = X_Min[0]+xScale+0.05;
  for(int i = 1; i < nParamPlot; i++)
  {
    X_Min[i] = X_Max[i-1];
    X_Max[i] = X_Min[i]+xScale;
  }
  double* Y_Min = new double[nParamPlot];
  double* Y_Max = new double[nParamPlot];
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
  int counterPad = 0;
  int counterText = 0;
  int counterPost = 0;
  int counter2DPost = 0;
  //KS: We start from top of the plot, might be confusing but works very well
  for(int y = 0; y < nParamPlot; y++)
  {
    //KS: start from left and go right, depending on y
    for(int x = 0; x <= y; x++)
    {
      //KS: Need to go to canvas every time to have our pads in the same canvas, not pads in the pads
      Posterior->cd();
      TrianglePad[counterPad] = new TPad(Form("TPad_%i", counterPad), Form("TPad_%i", counterPad), X_Min[x], Y_Min[y], X_Max[x], Y_Max[y]);

      TrianglePad[counterPad]->SetTopMargin(0);
      TrianglePad[counterPad]->SetRightMargin(0);

      TrianglePad[counterPad]->SetGrid();
      TrianglePad[counterPad]->SetFrameBorderMode(0);
      TrianglePad[counterPad]->SetBorderMode(0);
      TrianglePad[counterPad]->SetBorderSize(0);

      //KS: Corresponds to bottom part of the plot, need marings for lables
      if(y == (nParamPlot-1)) TrianglePad[counterPad]->SetBottomMargin(0.1);
      else TrianglePad[counterPad]->SetBottomMargin(0);

      //KS: Corresponds to left part, need marings for lables
      if(x == 0) TrianglePad[counterPad]->SetLeftMargin(0.15);
      else TrianglePad[counterPad]->SetLeftMargin(0);

      TrianglePad[counterPad]->Draw();
      TrianglePad[counterPad]->cd();

      //KS:if diagonal plot main posterior
      if(x == y)
      {
        hpost_copy[counterPost] = (TH1D*) hpost[ParamNumber[x]]->Clone(Form("hpost_copy_%i", ParamNumber[x]));
        hpost_cl[counterPost] = new TH1D*[nCredibleIntervals];
        /// Scale the histograms so it shows the posterior probability
        hpost_copy[counterPost]->Scale(1. / hpost_copy[counterPost]->Integral());
        for (int j = 0; j < nCredibleIntervals; ++j)
        {
          hpost_cl[counterPost][j] = (TH1D*) hpost[ParamNumber[x]]->Clone( Form("hpost_copy_%i_CL_%f", ParamNumber[x], CredibleIntervals[j]));
          //KS: Reset to get rid to TF1 otherwise we run into segfault :(
          hpost_cl[counterPost][j]->Reset("");
          hpost_cl[counterPost][j]->Fill(0.0, 0.0);

          // Scale the histograms before gettindg credible intervals
          hpost_cl[counterPost][j]->Scale(1. / hpost_cl[counterPost][j]->Integral());
          //KS: Slightly different approach depending if interavls are in pecentage or sigmas
          if(CredibleInSigmas)
          {
            //KS: Convert sigmas into percentage
            double CredReg = GetSigmaValue((int)std::round(CredibleIntervals[j]));
            GetCredibleInterval(hpost_copy[counterPost], hpost_cl[counterPost][j], CredReg);
          }
          else
          {
            GetCredibleInterval(hpost_copy[counterPost], hpost_cl[counterPost][j], CredibleIntervals[j]);
          }

          hpost_cl[counterPost][j]->SetFillColor(CredibleIntervalsColours[j]);
          hpost_cl[counterPost][j]->SetLineWidth(1);
        }

        hpost_copy[counterPost]->SetMaximum(hpost_copy[counterPost]->GetMaximum()*1.2);
        hpost_copy[counterPost]->SetLineWidth(2);
        hpost_copy[counterPost]->SetLineColor(kBlack);

        //KS: Don't want any titles
        hpost_copy[counterPost]->GetXaxis()->SetTitle("");
        hpost_copy[counterPost]->GetYaxis()->SetTitle("");
        hpost_copy[counterPost]->SetTitle("");

        hpost_copy[counterPost]->GetXaxis()->SetLabelSize(0.1);
        hpost_copy[counterPost]->GetYaxis()->SetLabelSize(0.1);

        hpost_copy[counterPost]->GetXaxis()->SetNdivisions(4);
        hpost_copy[counterPost]->GetYaxis()->SetNdivisions(4);

        hpost_copy[counterPost]->Draw("HIST");
        for (int j = 0; j < nCredibleIntervals; ++j)
          hpost_cl[counterPost][j]->Draw("HIST SAME");
        counterPost++;
      }
      //KS: Here we plot 2D credible regions
      else
      {
        hpost_2D_copy[counter2DPost] = (TH2D*) hpost2D[ParamNumber[x]][ParamNumber[y]]->Clone( Form("hpost_copy_%i_%i", ParamNumber[x], ParamNumber[y]));
        hpost_2D_cl[counter2DPost] = new TH2D*[nCredibleRegions];
        //KS: Now copy for every credible region
        for (int k = 0; k < nCredibleRegions; ++k)
        {
          hpost_2D_cl[counter2DPost][k] = (TH2D*)hpost2D[ParamNumber[x]][ParamNumber[y]]->Clone( Form("hpost_copy_%i_%i_CL_%f", ParamNumber[x], ParamNumber[y], CredibleRegions[k]));

          if(CredibleInSigmas)
          {
            //KS: Convert sigmas into percentage
            double CredReg = GetSigmaValue((int)std::round(CredibleRegions[k]));
            GetCredibleRegion(hpost_2D_cl[counter2DPost][k], CredReg);
          }
          else
          {
            GetCredibleRegion(hpost_2D_cl[counter2DPost][k], CredibleRegions[k]);
          }

          hpost_2D_cl[counter2DPost][k]->SetLineColor(CredibleRegionColor[k]);
          hpost_2D_cl[counter2DPost][k]->SetLineWidth(2);
          hpost_2D_cl[counter2DPost][k]->SetLineStyle(CredibleRegionStyle[k]);
        }

        //KS: Don't want any titles
        hpost_2D_copy[counter2DPost]->GetXaxis()->SetTitle("");
        hpost_2D_copy[counter2DPost]->GetYaxis()->SetTitle("");
        hpost_2D_copy[counter2DPost]->SetTitle("");

        hpost_2D_copy[counter2DPost]->GetXaxis()->SetLabelSize(0.1);
        hpost_2D_copy[counter2DPost]->GetYaxis()->SetLabelSize(0.1);

        hpost_2D_copy[counter2DPost]-> GetXaxis()->SetNdivisions(4);
        hpost_2D_copy[counter2DPost]-> GetYaxis()->SetNdivisions(4);
        hpost_2D_copy[counter2DPost]->Draw("COL");
        //Now credible regions
        for (int k = 0; k < nCredibleRegions; ++k)
          hpost_2D_cl[counter2DPost][k]->Draw("CONT3 SAME");
        counter2DPost++;
      }
      //KS: Corresponds to bottom part of the plot
      if(y == (nParamPlot-1))
      {
        Posterior->cd();
        TriangleText[counterText] = new TText(X_Min[x]+ (X_Max[x]-X_Min[x])/4, 0.04, hpost[ParamNumber[x]]->GetTitle());
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
        TriangleText[counterText] = new TText(0.04, Y_Min[y] + (Y_Max[y]-Y_Min[y])/4, hpost[ParamNumber[y]]->GetTitle());
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
  TLegend* legend = new TLegend(0.60, 0.7, 0.9, 0.9);
  legend->SetTextSize(0.03);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetLineStyle(0);
  legend->SetBorderSize(0);
  //KS: Legend is shared so just take first histograms
  for (int j = nCredibleIntervals-1; j >= 0; --j)
  {
    if(CredibleInSigmas)
      legend->AddEntry(hpost_cl[0][j], Form("%.0f#sigma Credible Interval", CredibleIntervals[j]), "f");
    else
      legend->AddEntry(hpost_cl[0][j], Form("%.0f%% Credible Interval", CredibleRegions[j]*100), "f");
  }
  for (int k = nCredibleRegions-1; k >= 0; --k)
  {
    if(CredibleInSigmas)
      legend->AddEntry(hpost_2D_cl[0][k], Form("%.0f#sigma Credible Region", CredibleRegions[k]), "l");
    else
      legend->AddEntry(hpost_2D_cl[0][k], Form("%.0f%% Credible Region", CredibleRegions[k]*100), "l");
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
  for(int i = 0; i < nParamPlot*2; i++) delete TriangleText[i];
  for(int i = 0; i < nParamPlot; i++)
  {
    delete hpost_copy[i];
    for (int j = 0; j < nCredibleIntervals; ++j)
    {
      delete hpost_cl[i][j];
    }
    delete[] hpost_cl[i];
  }
  for(int i = 0; i < Npad - nParamPlot; i++)
  {
    delete hpost_2D_copy[i];
    for (int j = 0; j < nCredibleRegions; ++j)
    {
      delete hpost_2D_cl[i][j];
    }
    delete[] hpost_2D_cl[i];
  }

  delete[] hpost_copy;
  delete[] hpost_cl;
  delete[] hpost_2D_copy;
  delete[] hpost_2D_cl;
  delete[] TrianglePad;
  delete[] TriangleText;
  delete[] X_Min;
  delete[] X_Max;
  delete[] Y_Min;
  delete[] Y_Max;
  delete legend;

  //KS: Restore margin
  Posterior->SetTopMargin(TopMargin);
  Posterior->SetLeftMargin(BottomMargin);
  Posterior->SetLeftMargin(LeftMargin);
  Posterior->SetRightMargin(RighMargin);
}

// **************************
// Scan the input trees
void MCMCProcessor::ScanInput() {
// **************************
  // Open the Chain
  Chain = new TChain("posteriors","posteriors");
  Chain->Add(MCMCFile.c_str());

  nEntries = Chain->GetEntries();
  
  //Only is suboptimality we might want to change it, therefore set it high enough so it doesn't affect other functionality
  UpperCut = nEntries+1;

  // Get the list of branches
  TObjArray* brlis = (TObjArray*)(Chain->GetListOfBranches());

  // Get the number of branches
  nBranches = brlis->GetEntries();

  BranchNames.reserve(nBranches);
  IamVaried.reserve(nBranches);
  ParamType.reserve(nBranches);

  // Set all the branches to off
  Chain->SetBranchStatus("*", false);

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  for (int i = 0; i < nBranches; i++)
  {
    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    //KS: Exclude paramer types
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

    // If we're on beam systematics
    if(bname.BeginsWith("xsec_")) 
    {
      BranchNames.push_back(bname);
      ParamType.push_back(kXSecPar);
      PlotXSec = true;
      nParam[kXSecPar]++;
    }
    else if (bname.BeginsWith("ndd_"))
    {
      BranchNames.push_back(bname);
      ParamType.push_back(kNDPar);
      PlotDet = true;
      nParam[kNDPar]++;
    }
    else if (bname.BeginsWith("skd_joint_"))
    {
      BranchNames.push_back(bname);
      ParamType.push_back(kFDDetPar);
      nParam[kFDDetPar]++;
    }
    else if (bname.BeginsWith("sin2th_") ||
          bname.BeginsWith("delm2_")  ||
          bname.BeginsWith("delta_")    )
    {
      BranchNames.push_back(bname);
      ParamType.push_back(kOSCPar);
      nParam[kOSCPar]++;
    }

    //KS: as a bonus get LogL systeamtic
    if (bname.BeginsWith("LogL_sample_")) {
      SampleName_v.push_back(bname);
      nSamples++;
    }
    else if (bname.BeginsWith("LogL_systematic_")) {
      SystName_v.push_back(bname);
      nSysts++;
    }
  }
  nDraw = BranchNames.size();
  // Read the input Covariances
  ReadInputCov();
  
  // Check order of parameter types
  ScanParameterOrder();
  MACH3LOG_INFO("************************************************");
  MACH3LOG_INFO("Scanning output branches...");
  MACH3LOG_INFO("# useful entries in tree: \033[1;32m {} \033[0m ", nDraw);
  MACH3LOG_INFO("# XSec params:  \033[1;32m {} starting at {} \033[0m ", nParam[kXSecPar] - nFlux, ParamTypeStartPos[kXSecPar]);
  MACH3LOG_INFO("# Flux params:   {}", nFlux);
  MACH3LOG_INFO("# Flux params:   {}", nFlux);
  MACH3LOG_INFO("# ND params:    \033[1;32m {} starting at {} \033[0m ", nParam[kNDPar] - nFlux, ParamTypeStartPos[kNDPar]);
  MACH3LOG_INFO("# FD params:    \033[1;32m {} starting at {} \033[0m ", nParam[kFDDetPar] - nFlux, ParamTypeStartPos[kFDDetPar]);
  MACH3LOG_INFO("# Osc params:   \033[1;32m {} starting at {} \033[0m ", nParam[kOSCPar] - nFlux, ParamTypeStartPos[kOSCPar]);
  MACH3LOG_INFO("************************************************");

  nSteps = Chain->GetMaximum("step");
  // Set the step cut to be 20%
  int cut = nSteps/5;
  SetStepCut(cut);
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
  Gauss = new TF1("gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",   -5, 5);

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
    (*Central_Value)(i) = _UNDEF_;
    (*Means)(i) = _UNDEF_;
    (*Errors)(i) = _UNDEF_;
    (*Means_Gauss)(i) = _UNDEF_;
    (*Errors_Gauss)(i) = _UNDEF_;
    (*Means_HPD)(i) = _UNDEF_;
    (*Errors_HPD)(i) = _UNDEF_;
    (*Errors_HPD_Positive)(i) = _UNDEF_;
    (*Errors_HPD_Negative)(i) = _UNDEF_;
    for (int j = 0; j < nDraw; ++j) {
      (*Covariance)(i, j) = _UNDEF_;
      (*Correlation)(i, j) = _UNDEF_;
    }
  } 
  
  hpost = new TH1D*[nDraw]();
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
TH1D* MCMCProcessor::MakePrefit() {
// *****************************
  
  if (OutputFile == nullptr) MakeOutputFile();

  TH1D *PreFitPlot = new TH1D("Prefit", "Prefit", nDraw, 0, nDraw);
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
      if ( CentralValueTemp != 0)
      {
        Central = ParamCentral[ParamEnum][ParamNo] / CentralValueTemp;
        Error = ParamErrors[ParamEnum][ParamNo]/CentralValueTemp;
      } else 
      {
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
    if(!PlotFlatPrior && ParamFlat[ParamEnum][ParamNo])
    {
      Error = 0.;
    }

    PreFitPlot->SetBinContent(i+1, Central);
    PreFitPlot->SetBinError(i+1, Error);
    PreFitPlot->GetXaxis()->SetBinLabel(i+1, ParamNames[ParamEnum][ParamNo]);
  }
  PreFitPlot->SetDirectory(0);

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
  if(nParam[kXSecPar] > 0)  ReadXSecFile();
  if(nParam[kNDPar] > 0) ReadNDFile();
  if(nParam[kFDDetPar] > 0) ReadFDFile();
  if(nParam[kOSCPar] > 0)   ReadOSCFile();
  //KS: Remove parameters which were removed
  RemoveParameters();
}

// **************************
// Read the output MCMC file and find what inputs were used
void MCMCProcessor::FindInputFiles() {
// **************************

  // Now read the MCMC file
  TFile *TempFile = new TFile(MCMCFile.c_str(), "open");

  // Get the settings for the MCMC
  TMacro *Config = (TMacro*)(TempFile->Get("MaCh3_Config"));
  if (Config == nullptr) {
    MACH3LOG_ERROR("Didn't find MaCh3_Config tree in MCMC file! {}", MCMCFile);
    TempFile->ls();
    throw;
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

  //CW: And the ND Covariance matrix
  CovPos[kNDPar].push_back(GetFromManager<std::string>(Settings["General"]["Systematics"]["NDCovFile"], "none"));
  if(CovPos[kNDPar].back() == "none")
  {
    MACH3LOG_WARN("Couldn't find NDCov branch in output");
    InputNotFound = true;
  }

  //CW: And the FD Covariance matrix
  CovPos[kFDDetPar].push_back(GetFromManager<std::string>(Settings["General"]["Systematics"]["FDCovFile"], "none"));
  if(CovPos[kFDDetPar].back() == "none")
  {
    MACH3LOG_WARN("Couldn't find FDCov branch in output");
    InputNotFound = true;
  }

  //CW: And the Osc Covariance matrix
  CovPos[kOSCPar].push_back(GetFromManager<std::string>(Settings["General"]["Systematics"]["OscCovFile"], "none"));
  if(CovPos[kOSCPar].back() == "none")
  {
    MACH3LOG_WARN("Couldn't find OscCov branch in output");
    InputNotFound = true;
  }
  if(InputNotFound) MaCh3Utils::PrintConfig(Settings);

  if (std::getenv("MACH3") != nullptr)
  {
    for(unsigned int i = 0; i < CovPos[kXSecPar].size(); i++)
      CovPos[kXSecPar][i].insert(0, std::string(std::getenv("MACH3"))+"/");

    for(unsigned int i = 0; i < CovPos[kNDPar].size(); i++)
      CovPos[kNDPar][i].insert(0, std::string(std::getenv("MACH3"))+"/");

    for(unsigned int i = 0; i < CovPos[kFDDetPar].size(); i++)
      CovPos[kFDDetPar][i].insert(0, std::string(std::getenv("MACH3"))+"/");

    for(unsigned int i = 0; i < CovPos[kOSCPar].size(); i++)
      CovPos[kOSCPar][i].insert(0, std::string(std::getenv("MACH3"))+"/");
  }

  // Delete the TTrees and the input file handle since we've now got the settings we need
  delete Config;

  // Delete the MCMCFile pointer we're reading
  TempFile->Close();
  delete TempFile;
}


// ***************
// Read the xsec file and get the input central values and errors
void MCMCProcessor::ReadXSecFile() {
// ***************
  YAML::Node XSecFile;
  XSecFile["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for(unsigned int i = 0; i < CovPos[kXSecPar].size(); i++)
  {
    YAML::Node YAMLDocTemp = YAML::LoadFile(CovPos[kXSecPar][i]);
    for (const auto& item : YAMLDocTemp["Systematics"]) {
      XSecFile["Systematics"].push_back(item);
    }
  }

  auto systematics = XSecFile["Systematics"];
  int i = 0;
  for (auto it = systematics.begin(); it != systematics.end(); ++it, ++i)
  {
    auto const &param = *it;

    // Push back the name
    std::string TempString = (param["Systematic"]["Names"]["FancyName"].as<std::string>());

    //KS:Reject particular parameter names, noticed that sometimes string comparison doesn't work becasue of some weird casting of TObjString into std::string. This is rare and sooner or later we move away from TObjString so this is fine
    bool rejected = false;
    for (unsigned int ik = 0; ik < ExcludedNames.size(); ++ik)
    {
      if (TempString.rfind(ExcludedNames.at(ik), 0) == 0)
      {
        const int Tracker = ParamTypeStartPos[kXSecPar] + i;
        BranchNames[Tracker] = "delete";
        nParam[kXSecPar]--;
        rejected = true;
        break;
      }
    }
    if(rejected) continue;
    ParamNames[kXSecPar].push_back(TempString);
    if(ParamNames[kXSecPar][i].BeginsWith("b_"))
    {
      IsXsec.push_back(false);
      nFlux++;
    } 
    else IsXsec.push_back(true);  
    
    ParamCentral[kXSecPar].push_back( param["Systematic"]["ParameterValues"]["PreFitValue"].as<double>() );
    ParamNom[kXSecPar].push_back( param["Systematic"]["ParameterValues"]["Generated"].as<double>() );
    ParamErrors[kXSecPar].push_back( param["Systematic"]["Error"].as<double>() );

    bool flat = false;
    if (param["Systematic"]["FlatPrior"]) { flat = param["Systematic"]["FlatPrior"].as<bool>(); }
    ParamFlat[kXSecPar].push_back( flat );
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
    throw;
  }
  NDdetFile->cd();

  TMatrixDSym *NDdetMatrix = (TMatrixDSym*)(NDdetFile->Get("nddet_cov"));
  TVectorD *NDdetNominal = (TVectorD*)(NDdetFile->Get("det_weights"));
  TDirectory *BinningDirectory = (TDirectory*)NDdetFile->Get("Binning")->Clone();

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
  while ((key = (TKey*)next()))
  {
    std::string name = std::string(key->GetName());
    TH2Poly* RefPoly = (TH2Poly*)BinningDirectory->Get((name).c_str())->Clone();
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
    MACH3LOG_ERROR("Couldn't find NDdetFile {}", CovPos[kFDDetPar].back());
    throw;
  }
  FDdetFile->cd();

  TMatrixDSym *FDdetMatrix = (TMatrixDSym*)(FDdetFile->Get("SKJointError_Erec_Total"));

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
  if(FancyPlotNames) ParamNames[kFDDetPar].back() = "Momentum Scale";

  FDdetFile->Close();
  delete FDdetFile;
  delete FDdetMatrix;
}

// ***************
// Read the Osc cov file and get the input central values and errors
void MCMCProcessor::ReadOSCFile() {
// ***************

  // Do the same for the ND280
  TFile *OscFile = new TFile(CovPos[kOSCPar].back().c_str(), "open");
  if (OscFile->IsZombie()) {
    std::cerr << "Couldn't find OSCFile " << CovPos[kOSCPar].back() << std::endl;
    throw;
  }
  OscFile->cd();

  TMatrixDSym *OscMatrix = (TMatrixDSym*)(OscFile->Get("osc_cov"));
  //KS: Osc nominal we can also set via config so there is danger that this will nor corrspond to what was used in the fit
  TVectorD *OscNominal = (TVectorD*)(OscFile->Get("osc_nom"));
  TObjArray* osc_param_names = (TObjArray*)(OscFile->Get("osc_param_names"));
  TVectorD* osc_flat_prior = (TVectorD*)OscFile->Get("osc_flat_prior");

  for (int i = 0; i < osc_flat_prior->GetNrows(); ++i)
  {
    ParamNom[kOSCPar].push_back( (*OscNominal)(i) );
    ParamCentral[kOSCPar].push_back( (*OscNominal)(i) );

    ParamErrors[kOSCPar].push_back( std::sqrt((*OscMatrix)(i,i)) );
    // Push back the name
    std::string TempString = std::string(((TObjString*)osc_param_names->At(i))->GetString());
    ParamNames[kOSCPar].push_back(TempString);

    ParamFlat[kOSCPar].push_back( (bool)((*osc_flat_prior)(i)) );
  }
  if(PlotJarlskog)
  {
    Chain->SetAlias("J_cp", "TMath::Sqrt(sin2th_13)*TMath::Sqrt(1.-sin2th_13)*TMath::Sqrt(1.-sin2th_13)*TMath::Sqrt(sin2th_12)*TMath::Sqrt(1.-sin2th_12)*TMath::Sqrt(sin2th_23)*TMath::Sqrt(1.-sin2th_23)*TMath::Sin(delta_cp)");
    BranchNames.push_back("J_cp");
    ParamType.push_back(kOSCPar);
    nParam[kOSCPar]++;
    nDraw++;

    //TODO we should actually calculate central value and prior error but leave it for now...
    ParamNom[kOSCPar].push_back( 0. );
    ParamCentral[kOSCPar].push_back( 0. );
    ParamErrors[kOSCPar].push_back( 1. );
    // Push back the name
    ParamNames[kOSCPar].push_back("J_cp");
    ParamFlat[kOSCPar].push_back( false );
  }
  OscFile->Close();
  delete OscFile;
}

// ***************
//This is bit messy as currently BranchNames is taken from actual Chain.root file while proper parameter name from matrix.root and right now we don't access them at the same time.
void MCMCProcessor::RemoveParameters() {
// ***************

  for(int i = 0; i < nDraw; i++)
  {
    if(BranchNames[i] == "delete")
    {
      BranchNames.erase(BranchNames.begin() + i);
      ParamType.erase(ParamType.begin() + i);
      nDraw--;
      i = 0;
    }
  }
}


// ***************
// Make the step cut from a string
void MCMCProcessor::SetStepCut(std::string Cuts) {
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

// **************************
//CW: Get the mean and RMS of a 1D posterior
void MCMCProcessor::GetArithmetic(TH1D * const hist, const int i) {
// **************************
  (*Means)(i) = hist->GetMean();
  (*Errors)(i) = hist->GetRMS();
}

// **************************
//CW: Get Gaussian characteristics
void MCMCProcessor::GetGaussian(TH1D *& hist , const int i) {
// **************************

  const double mean = hist->GetMean();
  const double err = hist->GetRMS();
  const double peakval = hist->GetBinCenter(hist->GetMaximumBin());

  // Set the range for the Gaussian fit
  Gauss->SetRange(mean - 1.5*err , mean + 1.5*err);
  // Set the starting parameters close to RMS and peaks of the histograms
  Gauss->SetParameters(hist->GetMaximum()*err*std::sqrt(2*3.14), peakval, err);

  // Perform the fit
  hist->Fit(Gauss->GetName(),"Rq");
  hist->SetStats(0);

  (*Means_Gauss)(i) = Gauss->GetParameter(1);
  (*Errors_Gauss)(i) = Gauss->GetParameter(2);
}


// ***************
//CW: Get the highest posterior density from a TH1D
void MCMCProcessor::GetHPD(TH1D* const hist, const int i, const double coverage) {
// ***************
  // Get the bin which has the largest posterior density
  const int MaxBin = hist->GetMaximumBin();
  // And it's value
  const double peakval = hist->GetBinCenter(MaxBin);

  // The total integral of the posterior
  const long double Integral = hist->Integral();
  //KS: and integral of left handed and right handed parts
  const long double LowIntegral = hist->Integral(1, MaxBin-1) + hist->GetBinContent(MaxBin)/2.0;
  const long double HighIntegral = hist->Integral(MaxBin+1, hist->GetNbinsX()) + hist->GetBinContent(MaxBin)/2.0;

  // Keep count of how much area we're covering
  //KS: Take only half content of HPD bin as one half goes for right handed error and the other for left handed error
  long double sum = hist->GetBinContent(MaxBin)/2.0;

  // Counter for current bin
  int CurrBin = MaxBin;
  while (sum/HighIntegral < coverage && CurrBin < hist->GetNbinsX()) {
    CurrBin++;
    sum += hist->GetBinContent(CurrBin);
  }
  const double sigma_p = std::fabs(hist->GetBinCenter(MaxBin)-hist->GetXaxis()->GetBinUpEdge(CurrBin));
  // Reset the sum
  //KS: Take only half content of HPD bin as one half goes for right handed error and the other for left handed error
  sum = hist->GetBinContent(MaxBin)/2.0;

  // Reset the bin counter
  CurrBin = MaxBin;
  // Counter for current bin
  while (sum/LowIntegral < coverage && CurrBin > 1) {
    CurrBin--;
    sum += hist->GetBinContent(CurrBin);
  }
  const double sigma_m = std::fabs(hist->GetBinCenter(CurrBin)-hist->GetBinLowEdge(MaxBin));

  // Now do the double sided HPD
  //KS: Start sum from the HPD
  sum = hist->GetBinContent(MaxBin);
  int LowBin = MaxBin;
  int HighBin = MaxBin;
  long double LowCon = 0.0;
  long double HighCon = 0.0;

  while (sum/Integral < coverage && (LowBin > 0 || HighBin < hist->GetNbinsX()+1))
  {
    LowCon = 0.0;
    HighCon = 0.0;
    //KS:: Move further only if you haven't reached histogram end
    if(LowBin > 1)
    {
      LowBin--;
      LowCon = hist->GetBinContent(LowBin);
    }
    if(HighBin < hist->GetNbinsX())
    {
      HighBin++;
      HighCon = hist->GetBinContent(HighBin);
    }

    // If we're on the last slice and the lower contour is larger than the upper
    if ((sum+LowCon+HighCon)/Integral > coverage && LowCon > HighCon) {
      sum += LowCon;
      break;
      // If we're on the last slice and the upper contour is larger than the lower
    } else if ((sum+LowCon+HighCon)/Integral > coverage && HighCon >= LowCon) {
      sum += HighCon;
      break;
    } else {
      sum += LowCon + HighCon;
    }
  }

  double sigma_hpd = 0.0;
  if (LowCon > HighCon) {
    sigma_hpd = std::fabs(hist->GetBinLowEdge(LowBin)-hist->GetBinCenter(MaxBin));
  } else {
    sigma_hpd = std::fabs(hist->GetXaxis()->GetBinUpEdge(HighBin)-hist->GetBinCenter(MaxBin));
  }

  (*Means_HPD)(i) = peakval;
  (*Errors_HPD)(i) = sigma_hpd;
  (*Errors_HPD_Positive)(i) = sigma_p;
  (*Errors_HPD_Negative)(i) = sigma_m;
}

// ***************
//KS: Get 1D histogram within credible interval, hpost_copy has to have the same binning, I don't do Copy() as this will lead to problems if this is used under multithreading
void MCMCProcessor::GetCredibleInterval(TH1D* const hist, TH1D* hpost_copy, const double coverage) {
// ***************

  if(coverage > 1)
  {
    MACH3LOG_ERROR("Specified Credible Interval is greater that 1 and equal to {} Should be between 0 and 1", coverage);
    throw;
  }
  //KS: Reset first copy of histogram
  hpost_copy->Reset("");
  hpost_copy->Fill(0.0, 0.0);

  //KS: Temporary structure to be thread save
  double *hist_copy = new double[hist->GetXaxis()->GetNbins()+1];
  bool *hist_copy_fill = new bool[hist->GetXaxis()->GetNbins()+1];
  for (int i = 0; i <= hist->GetXaxis()->GetNbins(); ++i)
  {
    hist_copy[i] = hist->GetBinContent(i);
    hist_copy_fill[i] = false;
  }

  /// Loop over histogram bins with highest number of entries until covered 90 or 68.3%
  const long double Integral = hist->Integral();
  long double sum = 0;

  while ((sum / Integral) < coverage)
  {
    /// Get bin of highest content and save the number of entries reached so far
    int max_entry_bin = 0;
    double max_entries = 0.;
    for (int i = 0; i <= hist->GetXaxis()->GetNbins(); ++i)
    {
      if (hist_copy[i] > max_entries)
      {
        max_entries = hist_copy[i];
        max_entry_bin = i;
      }
    }
    /// Replace bin value by -1 so it is not looped over as being maximum bin again
    hist_copy[max_entry_bin] = -1.;
    hist_copy_fill[max_entry_bin] = true;

    sum += max_entries;
  }
  //KS: Now fill our copy only for bins which got included in coverage region
  for(int i = 0; i <= hist->GetXaxis()->GetNbins(); ++i)
  {
    if(hist_copy_fill[i]) hpost_copy->SetBinContent(i, hist->GetBinContent(i));
  }

  delete[] hist_copy;
  delete[] hist_copy_fill;

  return;
}

// ***************
//KS: Set 2D contour within some coverage
void MCMCProcessor::GetCredibleRegion(TH2D* const hist2D, const double coverage) {
// ***************

  if(coverage > 1)
  {
    std::cerr<<"Specified Credible Region is greater than 1 and equal to "<< coverage <<" Should be between 0 and 1"<<std::endl;
    throw;
  }

  //KS: Temporary structure to be thread save
  double **hist_copy = new double*[hist2D->GetXaxis()->GetNbins()+1];
  for (int i = 0; i <= hist2D->GetXaxis()->GetNbins(); ++i)
  {
    hist_copy[i] = new double[hist2D->GetYaxis()->GetNbins()+1];
    for (int j = 0; j <= hist2D->GetYaxis()->GetNbins(); ++j)
    {
      hist_copy[i][j] = hist2D->GetBinContent(i,j);
    }
  }

  /// Loop over histogram bins with highest number of entries until covered 90 or 68.3%
  const long double Integral = hist2D->Integral();
  long double sum = 0;

  //We need to as ROOT requires array to set to contour
  double Contour[1];
  while ((sum / Integral) < coverage)
  {
    /// Get bin of highest content and save the number of entries reached so far
    int max_entry_bin_x = 0;
    int max_entry_bin_y = 0;
    double max_entries = 0.;
    for (int i = 0; i <= hist2D->GetXaxis()->GetNbins(); ++i)
    {
      for (int j = 0; j <= hist2D->GetYaxis()->GetNbins(); ++j)
      {
        if (hist_copy[i][j] > max_entries)
        {
          max_entries = hist_copy[i][j];
          max_entry_bin_x = i;
          max_entry_bin_y = j;
        }
      }
    }
    /// Replace bin value by -1 so it is not looped over as being maximum bin again
    hist_copy[max_entry_bin_x][max_entry_bin_y] = -1.;

    sum += max_entries;
    Contour[0] = max_entries;
  }
  hist2D->SetContour(1, Contour);

  //Delete temporary arrays
  for (int i = 0; i <= hist2D->GetXaxis()->GetNbins(); ++i)
  {
    delete[] hist_copy[i];
  }
  delete[] hist_copy;

  return;
}

// ***************
// Pass central value
void MCMCProcessor::GetNthParameter(const int param, double &Prior, double &PriorError, TString &Title){
// **************************

  ParameterEnum ParType = ParamType[param];
  int ParamNo = _UNDEF_;
  ParamNo = param - ParamTypeStartPos[ParType];

  Prior = ParamCentral[ParType][ParamNo];
  PriorError = ParamErrors[ParType][ParamNo];
  Title = ParamNames[ParType][ParamNo];
  return;
}

// ***************
// Find Param Index based on name
int MCMCProcessor::GetParamIndexFromName(const std::string Name){
// **************************

  int ParamNo = _UNDEF_;
  for (int i = 0; i < nDraw; ++i)
  {
    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;

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
void MCMCProcessor::GetPolarPlot(std::vector<std::string> ParNames){
// **************************

  if(hpost[0] == nullptr) MakePostfit();

  const double TopMargin = Posterior->GetTopMargin();
  const double BottomMargin = Posterior->GetBottomMargin();
  const double LeftMargin = Posterior->GetLeftMargin();
  const double RightMargin = Posterior->GetRightMargin();

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
    bool skip = false;
    if(ParamNo == _UNDEF_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Will not calculate Polar Plot", ParNames[k]);
      skip = true;
    }
    if(skip) continue;

    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;
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

    TGraphPolar * PolarGraph = new TGraphPolar(nBins, x_val.data(), y_val.data());
    PolarGraph->SetLineWidth(2);
    PolarGraph->SetFillStyle(3001);
    PolarGraph->SetLineColor(kRed);
    PolarGraph->SetFillColor(kRed);
    PolarGraph->Draw("AFL");

    TText* Text = new TText(0.6, 0.1, Title);
    Text->SetTextSize(0.04);
    Text->SetNDC(true);
    Text->Draw("");

    Posterior->Print(CanvasName);
    Posterior->Write(Title);

    delete PolarGraph;
    delete Text;
  } //End loop over parameters

  PolarDir->Close();
  delete PolarDir;

  OutputFile->cd();

  Posterior->SetTopMargin(TopMargin);
  Posterior->SetBottomMargin(BottomMargin);
  Posterior->SetLeftMargin(LeftMargin);
  Posterior->SetRightMargin(RightMargin);
}

// **************************
// Get Bayes Factor for particualar parameter
void MCMCProcessor::GetBayesFactor(std::vector<std::string> ParNames, std::vector<std::vector<double>> Model1Bounds, std::vector<std::vector<double>> Model2Bounds, std::vector<std::vector<std::string>> ModelNames){
// **************************

  if(hpost[0] == nullptr) MakePostfit();

  MACH3LOG_INFO("Calculating Bayes Factor");
  if((ParNames.size() != Model1Bounds.size()) || (Model2Bounds.size() != Model1Bounds.size())  || (Model2Bounds.size() != ModelNames.size()))
  {
    MACH3LOG_ERROR("Size doesn't match");
    throw;
  }
  for(unsigned int k = 0; k < ParNames.size(); ++k)
  {
    //KS: First we need to find parameter number based on name
    int ParamNo = GetParamIndexFromName(ParNames[k]);
    bool skip = false;
    if(ParamNo == _UNDEF_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Will not calculate Bayes Factor", ParNames[k]);
      skip = true;
    }
    if(skip) continue;

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
    MACH3LOG_INFO("Following Jeffreys Scale = ", JeffreysScale);
    MACH3LOG_INFO("Following Dunne-Kaboth Scale = ", DunneKabothScale);
    std::cout<<std::endl;
  }
  return;
}


// **************************
// KS: Get Savage Dockey point hypothesis test
void MCMCProcessor::GetSavageDickey(std::vector<std::string> ParNames, std::vector<double> EvaluationPoint, std::vector<std::vector<double>> Bounds){
// **************************

  if((ParNames.size() != EvaluationPoint.size()) || (Bounds.size() != EvaluationPoint.size()))
  {
    MACH3LOG_ERROR("Size doesn't match");
    throw;
  }
  
  if(hpost[0] == nullptr) MakePostfit();

  MACH3LOG_INFO("Calculating Savage Dickey");
  TDirectory *SavageDickeyDir = OutputFile->mkdir("SavageDickey");
  SavageDickeyDir->cd();
  
  for(unsigned int k = 0; k < ParNames.size(); ++k)
  {
    //KS: First we need to find parameter number based on name
    int ParamNo = GetParamIndexFromName(ParNames[k]);
    bool skip = false;
    if(ParamNo == _UNDEF_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Will not calculate SavageDickey", ParNames[k]);
      skip = true;
    }
    if(skip) continue;
    
    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;
    bool FlatPrior = false;
    GetNthParameter(ParamNo, Prior, PriorError, Title);
    
    ParameterEnum ParType = ParamType[ParamNo];
    int ParamTemp = ParamNo - ParamTypeStartPos[ParType];
    FlatPrior = ParamFlat[ParType][ParamTemp];
    
    TH1D* PosteriorHist = (TH1D*) hpost[ParamNo]->Clone(Title);
    RemoveFitter(PosteriorHist, "Gauss");
            
    TH1D* PriorHist = nullptr;
    //KS: If flat prior we need to have well defined bounds otherwise Prior distriution will not make sense
    if(FlatPrior)
    {
      int NBins = PosteriorHist->GetNbinsX();
      if(Bounds[k][0] > Bounds[k][1])
      {
        MACH3LOG_ERROR("Lower bound is higher than upper bound");
        throw;
      }
      PriorHist = new TH1D("PriorHist", Title, NBins, Bounds[k][0], Bounds[k][1]);
      
      double FlatProb = ( Bounds[k][1] - Bounds[k][0]) / NBins;
      for (int g = 0; g < NBins + 1; ++g) 
      {
        PriorHist->SetBinContent(g+1, FlatProb);
      }
    }
    else //KS: Otherwise throw from Gaussian
    {
      PriorHist = (TH1D*) PosteriorHist->Clone("Prior");
      PriorHist->Reset("");
      PriorHist->Fill(0.0, 0.0);
      
      TRandom3* rand = new TRandom3(0);
      //KS: Throw nice gaussian, just need big number to have smooth distribution
      for(int g = 0; g < 1000000; ++g)
      {
        PriorHist->Fill(rand->Gaus(Prior, PriorError));
      }
      delete rand;
    }
    // Area normalise the distributions
    PriorHist->Scale(1./PriorHist->Integral(), "width");
    PosteriorHist->Scale(1./PosteriorHist->Integral(), "width");
      
    PriorHist->SetLineColor(kRed);
    PriorHist->SetMarkerColor(kRed);
    PriorHist->SetFillColorAlpha(kRed, 0.35);
    PriorHist->SetFillStyle(1001);
    PriorHist->GetXaxis()->SetTitle(Title);
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

    double ProbPrior = PriorHist->GetBinContent(PriorHist->FindBin(EvaluationPoint[k]));
    //KS: In case we go so far away that prior is 0, set this to small value to avoid dividing by 0
    if(ProbPrior < 0) ProbPrior = 0.00001;
    double ProbPosterior = PosteriorHist->GetBinContent(PosteriorHist->FindBin(EvaluationPoint[k]));
    double SavageDickey = ProbPosterior/ProbPrior;
    
    std::string DunneKabothScale = GetDunneKaboth(SavageDickey);
    //Get Best point
    TGraph *PostPoint = new TGraph(1);
    PostPoint->SetPoint(0, EvaluationPoint[k], ProbPosterior);
    PostPoint->SetMarkerStyle(20);
    PostPoint->SetMarkerSize(1);
    PostPoint->Draw("P same");
    
    TGraph *PriorPoint = new TGraph(1);
    PriorPoint->SetPoint(0, EvaluationPoint[k], ProbPrior);
    PriorPoint->SetMarkerStyle(20);
    PriorPoint->SetMarkerSize(1);
    PriorPoint->Draw("P same");
    
    TLegend *legend = new TLegend(0.12, 0.6, 0.6, 0.97);
    legend->SetTextSize(0.04);
    legend->AddEntry(PriorHist, "Prior", "l");
    legend->AddEntry(PosteriorHist, "Posterior", "l");
    legend->AddEntry(PostPoint, Form("SavageDickey = %.2f, (%s)", SavageDickey, DunneKabothScale.c_str()),"");
    legend->SetLineColor(0);
    legend->SetLineStyle(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw("same");
  
    Posterior->Print(CanvasName);
    Posterior->Write(Title);
    
    delete PosteriorHist;
    delete PriorHist;
    delete PostPoint;
    delete PriorPoint;
    delete legend;
  } //End loop over parameters

  SavageDickeyDir->Close();
  delete SavageDickeyDir;

  OutputFile->cd();
}


// **************************
// KS: Reweight prior of MCMC chain to another
void MCMCProcessor::ReweightPrior(std::vector<std::string> Names, std::vector<double> NewCentral, std::vector<double> NewError){
// **************************

  MACH3LOG_INFO("Reweighting Prior");

  if( (Names.size() != NewCentral.size()) || (NewCentral.size() != NewError.size()))
  {
    MACH3LOG_ERROR("Size of passed vectors doesn't match in ReweightPrior");
    throw;
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
    if(ParamNo == _UNDEF_)
    {
      MACH3LOG_WARN("Couldn't find param {}. Can't reweight Prior", Names[k]);
      continue;
    }

    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;
    GetNthParameter(ParamNo, Prior, PriorError, Title);

    Param.push_back(ParamNo);
    OldCentral.push_back(Prior);
    OldError.push_back(PriorError);

    ParameterEnum ParType = ParamType[ParamNo];
    int ParamTemp = ParamNo - ParamTypeStartPos[ParType];

    FlatPrior.push_back(ParamFlat[ParType][ParamTemp]);
  }

  double* ParameterPos = new double[Names.size()];

  std::string InputFile = MCMCFile+".root";
  std::string OutputFilename = MCMCFile + "_reweighted.root";

  //KS: Simply create copy of file and add there new branch
  system(("cp "+InputFile+" "+OutputFilename).c_str());

  TFile *OutputChain = new TFile(OutputFilename.c_str(), "UPDATE");
  OutputChain->cd();
  TTree *post = (TTree *)OutputChain->Get("posteriors");

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
      if(FlatPrior[j])
      {
        old_prior = 1.0;
      }
      else
      {
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
  OutputChain->Close();
  delete OutputChain;
  delete[] ParameterPos;

  OutputFile->cd();
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
  AutoCorrelation();

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
    
  if(ParStep != nullptr)
  {
    MACH3LOG_ERROR("It look like ParStep was already filled ");
    MACH3LOG_ERROR("Eventhough it is used for MakeCovariance_MP and for DiagMCMC");
    MACH3LOG_ERROR("it has differnt structure in both for cache hits, sorry ");
    throw;
  }
  if(nBatches == 0)
  {
    MACH3LOG_ERROR("nBatches is equal to 0");
    MACH3LOG_ERROR("please use SetnBatches to set other value fore example 20");
    throw;      
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

  // Turn on the branches which we want for parameters
  for (int i = 0; i < nDraw; ++i) {
    Chain->SetBranchStatus(BranchNames[i].Data(), true);
  }
  
  // Turn on the branches which we want for LogL sample
  for (int i = 0; i < nSamples; ++i) {
    Chain->SetBranchStatus(SampleName_v[i].Data(), true);
  }

  // Turn on the branches which we want for LogL systs
  for (int i = 0; i < nSysts; ++i) {
    Chain->SetBranchStatus(SystName_v[i].Data(), true);
  }

  // Turn on the branches which we want for acc prob
  Chain->SetBranchStatus("accProb", true);
  
  // Only needed for Geweke right now
  Chain->SetBranchStatus("step", true);

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

  // Loop over the entries
  //KS: This is really a bottleneck right now, thus revisit with ROOT6 https://pep-root6.github.io/docs/analysis/parallell/root.html
  for (int i = 0; i < nEntries; ++i) {

    if (i % countwidth == 0)
      MaCh3Utils::PrintProgressBar(i, nEntries);

    // Set the branch addresses for params
    for (int j = 0; j < nDraw; ++j) {
      Chain->SetBranchAddress(BranchNames[j].Data(), &ParStep[j][i]);
    }

    // Set the branch addresses for samples
    for (int j = 0; j < nSamples; ++j) {
      Chain->SetBranchAddress(SampleName_v[j].Data(), &SampleValues[i][j]);
    }

    // Set the branch addresses for systematics
    for (int j = 0; j < nSysts; ++j) {
      Chain->SetBranchAddress(SystName_v[j].Data(), &SystValues[i][j]);
    }
      
    // Set the branch addresses for Acceptance Probability
    Chain->SetBranchAddress("accProb", &AccProbValues[i]);

    Chain->SetBranchAddress("step", &StepNumber[i]);

    // Fill up the arrays
    Chain->GetEntry(i);

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
    
    //KS: Could easyli add this to above loop but I accProb is different beast so better keep it like this
    AccProbBatchedAverages[BatchNumber] += AccProbValues[i];
  }

  clock.Stop();
  MACH3LOG_INFO("Took {:.2f}s to finish caching statistic for Diag MCMC with {} steps", clock.RealTime(), nEntries);

  // Make the sums into average
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nDraw; ++i) {
    ParamSums[i] /= nEntries;
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
  TH1D** TraceParamPlots = new TH1D*[nDraw];
  TH1D** TraceSamplePlots = new TH1D*[nSamples];
  TH1D** TraceSystsPlots = new TH1D*[nSysts];

  // Set the titles and limits for TH2Ds
  for (int j = 0; j < nDraw; ++j) {

    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;
    
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
    TF1 *Fitter = new TF1("Fitter","[0]", int(nEntries/2), nEntries);
    Fitter->SetLineColor(kRed);
    TraceParamPlots[j]->Fit("Fitter","Rq");
    TraceParamPlots[j]->Write();
    delete Fitter;
    delete TraceParamPlots[j];
  }
  delete[] TraceParamPlots;

  TDirectory *LLDir = OutputFile->mkdir("LogL");
  LLDir->cd();
  for (int j = 0; j < nSamples; ++j) {
    TraceSamplePlots[j]->Write();
    delete TraceSamplePlots[j];
    delete[] SampleValues[j];
  }
  delete[] TraceSamplePlots;
  delete[] SampleValues;

  for (int j = 0; j < nSysts; ++j) {
    TraceSystsPlots[j]->Write();
    delete TraceSystsPlots[j];
    delete SystValues[j];
  }
  delete[] TraceSystsPlots;
  delete[] SystValues;

  TraceDir->Close();
  delete TraceDir;

  OutputFile->cd();
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
  double **DenomSum = new double*[nDraw]();
  double **NumeratorSum = new double*[nDraw]();
  LagL = new double*[nDraw];
  for (int j = 0; j < nDraw; ++j) {
    DenomSum[j] = new double[nLags];
    NumeratorSum[j] = new double[nLags];
    LagL[j] = new double[nLags];
  }
  TH1D** LagKPlots = new TH1D*[nDraw];
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
    double Prior = 1.0;
    double PriorError = 1.0;
    
    GetNthParameter(j, Prior, PriorError, Title);
    std::string HistName = Form("%s_%s_Lag", Title.Data(), BranchNames[j].Data());
    LagKPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nLags, 0.0, nLags);
    LagKPlots[j]->GetXaxis()->SetTitle("Lag");
    LagKPlots[j]->GetYaxis()->SetTitle("Auto-correlation function");
  }
//KS: If CUDA is not enabled do calculations on CPU
#ifndef CUDA
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
  delete[] LagKPlots;

  //KS: This is different diagnostic however it relies on calculated Lag, thus we call it before we delete LagKPlots
  CalculateESS(nLags);

  for (int j = 0; j < nDraw; ++j) {
    delete[] NumeratorSum[j];
    delete[] DenomSum[j];
    delete[] LagL[j];
  }
  delete[] NumeratorSum;
  delete[] DenomSum;
  delete[] LagL;
  delete[] ParamSums;

  AutoCorrDir->Close();
  delete AutoCorrDir;

  OutputFile->cd();

  clock.Stop();
  MACH3LOG_INFO("Making auto-correlations took {:.2f}s", clock.RealTime());
}

#ifdef CUDA
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
// KS: calc Effective Sample Size Following https://mc-stan.org/docs/2_18/reference-manual/effective-sample-size-section.html
// Furthermore we calculate Sampling efficiency following https://kmh-lanl.hansonhub.com/talks/maxent00b.pdf
// Rule of thumb is to have efficiency above 25%
void MCMCProcessor::CalculateESS(const int nLags) {
// **************************

  if(LagL == nullptr)
  {
    MACH3LOG_ERROR("Trying to call CalculateESS before LagL was calcauted, this will not work");
    std::cerr <<__FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  MACH3LOG_INFO("Making ESS plots...");
    
  TVectorD* EffectiveSampleSize = new TVectorD(nDraw);
  TVectorD* SamplingEfficiency = new TVectorD(nDraw);
  double *TempDenominator = new double[nDraw]();

  const int Nhists = 5;
  const double Thresholds[Nhists+1] = {1, 0.02, 0.005, 0.001, 0.0001, 0.0};
  const Color_t ESSColours[Nhists] = {kGreen, kGreen+2, kYellow, kOrange, kRed};

  //KS: This histogram is inspired by the following: https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
  TH1D **EffectiveSampleSizeHist = new TH1D*[Nhists]();
  for(int i = 0; i < Nhists; ++i)
  {
    EffectiveSampleSizeHist[i] = new TH1D(Form("EffectiveSampleSizeHist_%i",i), Form("EffectiveSampleSizeHist_%i",i), nDraw, 0, nDraw);
    EffectiveSampleSizeHist[i]->GetYaxis()->SetTitle("N_{eff}/N");
    EffectiveSampleSizeHist[i]->SetFillColor(ESSColours[i]);
    EffectiveSampleSizeHist[i]->SetLineColor(ESSColours[i]);
    EffectiveSampleSizeHist[i]->Sumw2();
    for (int j = 0; j < nDraw; ++j)
    {
      TString Title = "";
      double Prior = 1.0;
      double PriorError = 1.0;
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
    (*EffectiveSampleSize)(j) = _UNDEF_;
    (*SamplingEfficiency)(j) = _UNDEF_;
    TempDenominator[j] = 0.;
    //KS: Firs sum over all Calculated autocorrelations
    for (int k = 0; k < nLags; ++k)
    {
      TempDenominator[j] += LagL[j][k];
    }
    TempDenominator[j] = 1+2*TempDenominator[j];
    (*EffectiveSampleSize)(j) = nEntries/TempDenominator[j];
    // 100 because we convert to percentage
    (*SamplingEfficiency)(j) = 100 * 1/TempDenominator[j];

    for(int i = 0; i < Nhists; ++i)
    {
      EffectiveSampleSizeHist[i]->SetBinContent(j+1, 0);
      EffectiveSampleSizeHist[i]->SetBinError(j+1, 0);

      const double TempEntry = std::fabs((*EffectiveSampleSize)(j)) / nEntries;
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

  TLegend *leg = new TLegend(0.2, 0.7, 0.6, 0.95);
  leg->SetTextSize(0.03);
  for(int i = 0; i < Nhists; ++i)
  {
    leg->AddEntry(EffectiveSampleSizeHist[i], Form("%.4f >= N_{eff}/N > %.4f", Thresholds[i], Thresholds[i+1]), "f");
  }
  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");

  Posterior->Write("EffectiveSampleSizeCanvas");

  //Delete all variables
  delete EffectiveSampleSize;
  delete SamplingEfficiency;
  for(int i = 0; i < Nhists; ++i)
  {
    delete EffectiveSampleSizeHist[i];
  }
  delete leg;
  delete[] EffectiveSampleSizeHist;
  //KS Remove auxiliary arrays
  delete[] TempDenominator;
}

// **************************
//CW: Batched means, literally read from an array and chuck into TH1D
void MCMCProcessor::BatchedMeans() {
// **************************

  if (BatchedAverages == nullptr) PrepareDiagMCMC();
  MACH3LOG_INFO("Making BatchedMeans plots...");
  
  TH1D ** BatchedParamPlots = new TH1D*[nDraw];
  for (int j = 0; j < nDraw; ++j) {
    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;
    
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
  delete[] BatchedParamPlots;

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
    throw;
  }

  // Calculate variance estimator using batched means following https://arxiv.org/pdf/1911.00915.pdf see Eq. 1.2
  TVectorD* BatchedVariance = new TVectorD(nDraw);
  //KS: The hypothesis is rejected if C > z α for a given confidence level α. If the batch means do not pass the test, Correlated is reported for the half-width on the statistical reports following https://rossetti.github.io/RossettiArenaBook/ch5-BatchMeansMethod.html alternatively for more old-school see Alexopoulos and Seila 1998 section 3.4.3
  TVectorD* C_Test_Statistics = new TVectorD(nDraw);
 
  double* OverallBatchMean = new double[nDraw]();
  double* C_Rho_Nominator = new double[nDraw]();
  double* C_Rho_Denominator = new double[nDraw]();
  double* C_Nominator = new double[nDraw]();
  double* C_Denominator = new double[nDraw]();
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
  delete[] OverallBatchMean;
  delete[] C_Rho_Nominator;
  delete[] C_Rho_Denominator;
  delete[] C_Nominator;
  delete[] C_Denominator;
}

// **************************
// RC: Perform spectral analysis of MCMC based on http://arxiv.org/abs/astro-ph/0405462
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
  //KS: WARNING Code is awfully slow... I know how to make it faster (GPU scream in a distant) but for now just make it for two params, bit hacky sry...
  nPrams = 2;

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
      a_j /= float(std::sqrt(float(_N)));
      const int _c = jj - start;

      k_j[j][_c] = two_pi_over_N * jj;
      // Equation 13
      P_j[j][_c] = std::norm(a_j);
    }
  }

  TDirectory *PowerDir = OutputFile->mkdir("PowerSpectrum");
  PowerDir->cd();

  TGraph **plot = new TGraph*[nPrams];
  TVectorD* PowerSpectrumStepSize = new TVectorD(nPrams);
  for (int j = 0; j < nPrams; ++j)
  {
    plot[j] = new TGraph(v_size, k_j[j].data(), P_j[j].data());

    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;

    GetNthParameter(j, Prior, PriorError, Title);

    std::string name = Form("Power Spectrum of %s;k;P(k)", Title.Data());

    plot[j]->SetTitle(name.c_str());
    name = Form("%s_power_spectrum", Title.Data());
    plot[j]->SetName(name.c_str());
    plot[j]->SetMarkerStyle(7);

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

    plot[j]->Fit("power_template","Rq");

    Posterior->SetLogx();
    Posterior->SetLogy();
    Posterior->SetGrid();
    plot[j]->Write(plot[j]->GetName());
    plot[j]->Draw("AL");
    func->Draw("SAME");
    if(printToPDF) Posterior->Print(CanvasName);

    //KS: I have no clue what is the reason behind this. Found this in Rick Calland code...
    (*PowerSpectrumStepSize)(j) = std::sqrt(func->GetParameter(0)/float(v_size*0.5));

    delete func;
    delete plot[j];
  }
  delete [] plot;

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
// https://www.math.arizona.edu/~piegorsch/675/GewekeDiagnostics.pdf
// https://www2.math.su.se/matstat/reports/master/2011/rep2/report.pdf Chapter 3.1
void MCMCProcessor::GewekeDiagnostic() {
// **************************

  MACH3LOG_INFO("Making Geweke Diagnostic");

  //KS: Up refers to upper limit we check, it stays constnt, in literature it is moslty 50% thus using 0.5 for threshold
  double* MeanUp = new double[nDraw]();
  double* SpectralVarianceUp = new double[nDraw]();
  int* DenomCounterUp = new int[nDraw]();
  const double Threshold = 0.5 * nSteps;

  //KS: Select values betwen which you want to scan, for example 0 means 0% burn in and 1 100% burn in.
  const double LowerThreshold = 0;
  const double UpperThreshold = 1.0;
  // Tells how many intervals between thresholds we want to check
  const int NChecks = 100;
  const double Division = (UpperThreshold - LowerThreshold)/NChecks;

  TH1D** GewekePlots = new TH1D*[nDraw];
  for (int j = 0; j < nDraw; ++j)
  {
    MeanUp[j] = 0;
    SpectralVarianceUp[j] = 0;
    DenomCounterUp[j] = 0;

    TString Title = "";
    double Prior = 1.0;
    double PriorError = 1.0;
    GetNthParameter(j, Prior, PriorError, Title);
    std::string HistName = Form("%s_%s_Geweke", Title.Data(), BranchNames[j].Data());
    GewekePlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), NChecks, 0.0, 100*UpperThreshold);
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
    double* MeanDown = new double[nDraw]();
    double* SpectralVarianceDown = new double[nDraw]();
    int* DenomCounterDown = new int[nDraw]();

    //set to 0
    for (int j = 0; j < nDraw; ++j)
    {
      MeanDown[j] = 0;
      SpectralVarianceDown[j] = 0;
      DenomCounterDown[j] = 0;
    }

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
    //KS: delete for each thread
    delete[] MeanDown;
    delete[] SpectralVarianceDown;
    delete[] DenomCounterDown;
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
    delete GewekePlots[j];
  }
  delete[] GewekePlots;

  //Free memory
  delete[] MeanUp;
  delete[] DenomCounterUp;
  delete[] SpectralVarianceUp;

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
  TH1D* AcceptanceProbPlot = new TH1D("AcceptanceProbability", "Acceptance Probability", nEntries, 0, nEntries);
  AcceptanceProbPlot->GetXaxis()->SetTitle("Step");
  AcceptanceProbPlot->GetYaxis()->SetTitle("Acceptance Probability");

  TH1D* BatchedAcceptanceProblot = new TH1D("AcceptanceProbability_Batch", "AcceptanceProbability_Batch", nBatches, 0, nBatches);
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
    // Set bin content for the ith bin to the parameter values
    AcceptanceProbPlot->SetBinContent(i, AccProbValues[i]);
  }
    
  TDirectory *probDir = OutputFile->mkdir("AccProb");
  probDir->cd();
  
  AcceptanceProbPlot->Write();
  BatchedAcceptanceProblot->Write();
    
  delete AcceptanceProbPlot;  
  delete BatchedAcceptanceProblot; 
  delete[] AccProbValues;
  delete[] AccProbBatchedAverages;

  probDir->Close();
  delete probDir;

  OutputFile->cd();

}
