#include "MCMCProcessor.h"
MCMCProcessor::MCMCProcessor(const std::string &InputFile, bool MakePostfitCorr) : 
  Chain(NULL), StepCut(""), MakeCorr(MakePostfitCorr), MadePostfit(false) {

  MCMCFile = InputFile;

  std::cout << "Making post-fit processor for " << MCMCFile << std::endl;

  Posterior = NULL;
  
  //KS:Hardcoded should be a way to get it via config or something
  PlotDet = false;
  MakeOnlyXsecCorr = false;
  MakeOnlyXsecCorrFlux = false;
  plotRelativeToPrior = false;
  printToPDF = false;
  plotBinValue = false;
  CacheMCMCM = false;
  FancyPlotNames = true;
  doDiagMCMC = false;
  OutputSuffix = "_Process";

  nDraw = 0;
  nFlux = 0;
  nEntries = 0;
  nBatches = 0;
  nSysts = 0;
  nSamples = 0;
  
  nBins = 70;
  DrawRange = 1.5;
  
  //KS:Those keep basic information for ParameterEnum
  ParamNames.resize(kNParameterEnum);
  ParamCentral.resize(kNParameterEnum);
  ParamNom.resize(kNParameterEnum);
  ParamErrors.resize(kNParameterEnum);
  ParamTypeStartPos.resize(kNParameterEnum);
  nParam.resize(kNParameterEnum);
  CovPos.resize(kNParameterEnum);
  
  for(int i =0; i<kNParameterEnum; i++)
  {
     ParamTypeStartPos[i] = 0;
     nParam[i] = 0;
     CovPos[i] = "";
  }
}

// ****************************
// The destructor
MCMCProcessor::~MCMCProcessor() {
// ****************************

  // Close the pdf file
  std::cout << "Closing pdf in MCMCProcessor " << CanvasName << std::endl;
  CanvasName += "]";
  if(printToPDF) Posterior->Print(CanvasName);
  if (Posterior != NULL)  delete Posterior;

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
  
  if(CacheMCMCM)
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
    delete[] Min_Chain;
    delete[] Max_Chain;
    delete[] hpost;
    delete[] hpost2D;
  }
  if (OutputFile != NULL) OutputFile->Close();
  if (OutputFile != NULL) delete OutputFile;
  delete Chain;
}


void MCMCProcessor::Initialise(){
// ***************
  // Scan the ROOT file for useful branches
  ScanInput();

  // Setup the output
  SetupOutput();
}

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
// Get post-fits for the ParameterEnum type, e.g. xsec params, ND280 params or flux params
void MCMCProcessor::GetPostfit_Ind(TVectorD *&PDF_Central, TVectorD *&PDF_Errors, TVectorD *&Peak_Values, ParameterEnum kParam) {
// ***************
  // Make the post fit
  MakePostfit();

  // Loop over the loaded param types
  int ParamTypeSize = ParamType.size();
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
  if(MakeOnlyXsecCorr) std::cout<<"Will plot only xsec covariance"<<std::endl;
  if(MakeOnlyXsecCorrFlux) std::cout<<"Will plot only xsec covariances and flux"<<std::endl;

  MakeCovariance();
  Cov = Covariance;
  Corr = Correlation;
}

// ***************
void MCMCProcessor::MakeOutputFile() {
// ***************

  // Open a TCanvas to write the posterior onto
  Posterior = new TCanvas("Posterior", "Posterior", 0, 0, 1024, 1024);
  Posterior->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Posterior->SetTickx();
  Posterior->SetTicky();
  Posterior->SetBottomMargin(0.1);
  Posterior->SetTopMargin(0.05);
  Posterior->SetRightMargin(0.03);
  Posterior->SetLeftMargin(0.10);
  
  // Output file to write to
  OutputName = MCMCFile + OutputSuffix +".root";

  // Output file
  OutputFile = new TFile(OutputName.c_str(), "recreate");
  OutputFile->cd();
}


// ****************************
// Function to make the post-fit
void MCMCProcessor::MakePostfit() {
// ****************************

  // Check if we've already made post-fit
  if (MadePostfit == true) return;
  MadePostfit = true;

  // Check if the output file is ready
  if (OutputFile == NULL) MakeOutputFile();
  
  std::cout << "MCMCProcessor is making post-fit plots..." << std::endl;

  PostDir = OutputFile->mkdir("Post");

  
  // We fit with this Gaussian
  Gauss = new TF1("Gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",-5,5);
  Gauss->SetLineWidth(2);
  Gauss->SetLineColor(kOrange-5);
  
  // nDraw is number of draws we want to do
  for (int i = 0; i < nDraw; ++i)
  {
    if (i % (nDraw/5) == 0) {
      std::cout << "  " << i << "/" << nDraw << " (" << int((double(i)/double(nDraw)*100.0))+1 << "%)" << std::endl;
    }

    OutputFile->cd();
    TString Title = BranchNames[i];
    double Nominal = 1.0;
    double NominalError = 1.0;
    
    GetNthParameter(i, Nominal, NominalError, Title);

    OutputFile->cd();
    // This holds the posterior density
    double maxi = Chain->GetMaximum(BranchNames[i]);
    double mini = Chain->GetMinimum(BranchNames[i]);
    // This holds the posterior density
    hpost[i] = new TH1D(BranchNames[i], BranchNames[i], nBins, mini, maxi);
    hpost[i]->SetMinimum(0);
    hpost[i]->GetYaxis()->SetTitle("Steps");
    hpost[i]->GetYaxis()->SetNoExponent(false);
    // Project BranchNames[i] onto hpost, applying stepcut
    Chain->Project(BranchNames[i], BranchNames[i], StepCut.c_str());

    hpost[i]->Smooth();

    for(int ik = 0; ik < kNParameterEnum; ik++)
    {
        if (ParamType[i] == ParameterEnum(ik))
        {
            int ParamNo = __UNDEF__;
            ParamNo = i - ParamTypeStartPos[ParameterEnum(ik)];

            (*Central_Value)(i)  = ParamCentral[ParameterEnum(ik)][ParamNo];
        }
    }

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
    TLine *Asimov = new TLine(Nominal, hpost[i]->GetMinimum(), Nominal, hpost[i]->GetMaximum());
    Asimov->SetLineColor(kRed-3);
    Asimov->SetLineWidth(2);
    Asimov->SetLineStyle(kDashed);

    TLegend *leg = new TLegend(0.12, 0.6, 0.6, 0.97);
    leg->SetTextSize(0.04);
    leg->AddEntry(hpost[i], Form("#splitline{PDF}{#mu = %.2f, #sigma = %.2f}", hpost[i]->GetMean(), hpost[i]->GetRMS()), "l");
    leg->AddEntry(Gauss, Form("#splitline{Gauss}{#mu = %.2f, #sigma = %.2f}", Gauss->GetParameter(1), Gauss->GetParameter(2)), "l");
    leg->AddEntry(hpd, Form("#splitline{HPD}{#mu = %.2f, #sigma = %.2f (+%.2f-%.2f)}", (*Means_HPD)(i), (*Errors_HPD)(i), (*Errors_HPD_Positive)(i), (*Errors_HPD_Negative)(i)), "l");
    leg->AddEntry(Asimov, Form("#splitline{Asimov}{x = %.2f , #sigma = %.2f}", Nominal, NominalError), "l");
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);

    // Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
    if (hpost[i]->GetMaximum() == hpost[i]->Integral()*DrawRange) 
    {
        std::cout << "Found fixed parameter, moving on" << std::endl;
        IamVaried[i] = false;
        //KS:Set mean and error to prior for fixed parameters, this is mostly for comaprison with BANFF
        //but it looks much better when fixed parameter has mean on prior rather than on 0 with 0 error.
        for(int ik = 0; ik < kNParameterEnum; ik++)
        {
            if (ParamType[i] == ParameterEnum(ik)) 
            {
                int ParamNo = __UNDEF__;
                ParamNo = i - ParamTypeStartPos[ParameterEnum(ik)];
                
                (*Means_HPD)(i)  = ParamCentral[ParameterEnum(ik)][ParamNo];
                (*Errors_HPD)(i) = ParamErrors[ParameterEnum(ik)][ParamNo];
            }
        }
        delete Asimov;
        delete hpd;
        delete leg;
        continue;
    }

    // Store that this parameter is indeed being varied
    IamVaried[i] = true;

    // Write to file
    Posterior->SetName(hpost[i]->GetName());
    Posterior->SetTitle(hpost[i]->GetTitle());

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
  int FluxParameters = nFlux;
  SettingsBranch->Branch("FluxParameters", &FluxParameters);
  int NDParameters = nParam[kND280Par];
  SettingsBranch->Branch("NDParameters", &NDParameters);
  int FDParameters = nParam[kFDDetPar];
  SettingsBranch->Branch("FDParameters", &FDParameters);
  int OscParameters = nParam[kOSCPar];
  SettingsBranch->Branch("OscParameters", &OscParameters);

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

} // Have now written the postfit projections

// *******************
// Draw the postfit
void MCMCProcessor::DrawPostfit() {
// *******************
//
  if (OutputFile == NULL) MakeOutputFile();

  // Make the prefit plot
  TH1D* prefit = MakePrefit();

  // cd into the output file
  OutputFile->cd();
 
  // Make a TH1D of the central values and the errors
  TH1D *paramPlot = new TH1D("paramPlot", "paramPlot", nDraw, 0, nDraw);
  paramPlot->SetName("mach3params");
  paramPlot->SetTitle(StepCut.c_str());
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
    int ParamNo = __UNDEF__;
    int ParamEnu = __UNDEF__;
    for(int ik = 0; ik < kNParameterEnum; ik++)
    {
        if (ParamType[i] == ParameterEnum(ik)) 
        {
            ParamNo = i - ParamTypeStartPos[ParameterEnum(ik)];
            ParamEnu = ParameterEnum(ik);
        }
    }
      
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
    //KS: Just get value of each parmeter without dividing by prior
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

  Posterior->SetBottomMargin(0.2);

  OutputFile->cd();
  //KS: Plot Xsec and Flux
  if (PlotXSec == true) {
      
    int Start = ParamTypeStartPos[kXSecPar];
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
    int Start = ParamTypeStartPos[kND280Par];
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
      paramPlot->SetTitle(StepCut.c_str());
      paramPlot->GetXaxis()->LabelsOption("v");

      paramPlot_Gauss->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      paramPlot_Gauss->GetYaxis()->SetRangeUser(0.6, 1.4);
      paramPlot_Gauss->GetXaxis()->SetRangeUser(Start, NDbinCounter);
      paramPlot_Gauss->GetXaxis()->SetTitle("");
      paramPlot_Gauss->SetTitle(StepCut.c_str());
      paramPlot_Gauss->GetXaxis()->LabelsOption("v");

      paramPlot_HPD->GetYaxis()->SetTitle(("Variation for "+NDname).c_str());
      paramPlot_HPD->GetYaxis()->SetRangeUser(0.6, 1.4);
      paramPlot_HPD->GetXaxis()->SetRangeUser(Start, NDbinCounter);
      paramPlot_HPD->GetXaxis()->SetTitle("");
      paramPlot_HPD->SetTitle(StepCut.c_str());
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

  delete paramPlot;
  delete CompLeg;

  //KS: Return Margin to default one
  Posterior->SetBottomMargin(0.1);
}

// *********************
// Make the post-fit covariance matrix in all dimensions
void MCMCProcessor::MakeCovariance() {
// *********************

  if (OutputFile == NULL) MakeOutputFile();

  bool HaveMadeDiagonal = false;
  std::cout << "Making post-fit covariances..." << std::endl;

  // Check that the diagonal entries have been filled
  // i.e. MakePostfit() has been called
  for (int i = 0; i < nDraw; ++i) {
    if ((*Covariance)(i,i) == __UNDEF__) {
      HaveMadeDiagonal = false;
      std::cout << "Have not run diagonal elements in covariance, will do so now by calling MakePostfit()" << std::endl;
      break;
    } else {
      HaveMadeDiagonal = true;
    }
  }

  if (HaveMadeDiagonal == false) {
    MakePostfit();
  }

  int covBinning = nDraw;
  //If we only plot correlation/covariance between xsec (without flux)
  if(MakeOnlyXsecCorr)
  {
     covBinning = nParam[kXSecPar] - nFlux; 
  }
    
    
  // Now we are sure we have the diagonal elements, let's make the off-diagonals
  for (int i = 0; i < covBinning; ++i) {

    if (i % (covBinning/5) == 0) {
      std::cout << "  " << i << "/" << covBinning << " (" << int((double(i)/double(covBinning)*100.0))+1 << "%)" << std::endl;
    }
      
    TString Title_i = BranchNames[i];
    double Nominal_i, NominalError;

    GetNthParameter(i, Nominal_i, NominalError, Title_i);
    
    double min_i = Chain->GetMinimum(BranchNames[i]);
    double max_i = Chain->GetMaximum(BranchNames[i]);

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

      TString Title_j = BranchNames[j];
      double Nominal_j, NominalError_j;
      GetNthParameter(j, Nominal_j, NominalError_j, Title_j);

      OutputFile->cd();

      // The draw which we want to perform
      TString DrawMe = BranchNames[j]+":"+BranchNames[i];

      double max_j = Chain->GetMaximum(BranchNames[j]);
      double min_j = Chain->GetMinimum(BranchNames[j]);

      // TH2F to hold the Correlation 
      TH2D *hpost_2D = new TH2D(DrawMe, DrawMe, nBins, min_i, max_i, nBins, min_j, max_j);

      hpost_2D->SetMinimum(0);
      hpost_2D->GetXaxis()->SetTitle(Title_i);
      hpost_2D->GetYaxis()->SetTitle(Title_j);
      hpost_2D->GetZaxis()->SetTitle("Steps");

      // The draw command we want, i.e. draw param j vs param i
      Chain->Project(DrawMe, DrawMe, StepCut.c_str());
      
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
          if(IsXsec[j] && IsXsec[i])
          {
              hpost_2D->Draw("colz");
              Posterior->SetName(hpost_2D->GetName());
              Posterior->SetTitle(hpost_2D->GetTitle());
              Posterior->Print(CanvasName);
          }
          }
      }
      // Write it to root file
      //OutputFile->cd();
      //hpost_2D->Write();

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
    CacheMCMCM = true;
    std::cout << "Caching input tree..." << std::endl;
    TStopwatch clock;
    clock.Start();
    
    ParStep = new double*[nDraw];
    StepNumber = new int[nEntries];
    
    Min_Chain = new double[nDraw];
    Max_Chain= new double[nDraw];
    
    hpost2D = new TH2D**[nDraw]();

    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < nDraw; ++i) 
    {
        ParStep[i] = new double[nEntries];
        hpost2D[i] = new TH2D*[nDraw]();
        Min_Chain[i] = -999.99;
        Max_Chain[i] = -999.99;

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

    int countwidth = nEntries/10;
    // Loop over the entries
    for (int j = 0; j < nEntries; ++j) {
        if (j % countwidth == 0) {
        std::cout << j << "/" << nEntries << " (" << double(j)/double(nEntries)*100. << "%)" << std::endl;
        }

        Chain->SetBranchAddress("step", &StepNumber[j]);
        // Set the branch addresses for params
        for (int i = 0; i < nDraw; ++i) 
        {
            Chain->SetBranchAddress(BranchNames[i].Data(), &ParStep[i][j]);
        }
        
        // Fill up the ParStep array
        Chain->GetEntry(j);
    }
    
    // Set all the branches to on
    Chain->SetBranchStatus("*", true);
    
    // Cache max and min in chain for covariance matrix
    for (int i = 0; i < nDraw; ++i) 
    {
        Min_Chain[i] = Chain->GetMinimum(BranchNames[i]);
        Max_Chain[i] = Chain->GetMaximum(BranchNames[i]);
        for (int j = 0; j <= i; ++j)
        {
            // TH2D to hold the Correlation 
            hpost2D[i][j] = new TH2D(Form("hpost2D_%i_%i",i,j), Form("hpost2D_%i_%i",i,j), nBins, Min_Chain[i], Max_Chain[i], nBins, Min_Chain[j], Max_Chain[j]);
        }
    }
        
    clock.Stop();
    std::cout << "Caching steps took " << clock.RealTime() << "s to finish for " << nEntries << " events" << std::endl;
}


// *********************
// Make the post-fit covariance matrix in all dimensions
void MCMCProcessor::MakeCovariance_MP() {
// *********************
    
  if (OutputFile == NULL) MakeOutputFile();
    
  int covBinning = nDraw;
  //If we only plot correlation/covariance between xsec (without flux)
  if(MakeOnlyXsecCorr)
  {
     covBinning = nParam[kXSecPar] - nFlux; 
  }

  bool HaveMadeDiagonal = false;
  std::cout << "Making post-fit covariances..." << std::endl;
    
  // Check that the diagonal entries have been filled
  // i.e. MakePostfit() has been called
  for (int i = 0; i < covBinning; ++i) {
    if ((*Covariance)(i,i) == __UNDEF__) {
      HaveMadeDiagonal = false;
      std::cout << "Have not run diagonal elements in covariance, will do so now by calling MakePostfit()" << std::endl;
      break;
    } else {
      HaveMadeDiagonal = true;
    }
  }
    
 if (HaveMadeDiagonal == false) MakePostfit();
  
std::cout << "Calculating covaraince matrix" << std::endl;
TStopwatch clock;
clock.Start();
  
// Now we are sure we have the diagonal elements, let's make the off-diagonals
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
for (int i = 0; i < covBinning; ++i) 
{
    TString Title_i = BranchNames[i];
    double Nominal_i, NominalError_i;
    GetNthParameter(i, Nominal_i, NominalError_i, Title_i);
    
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
      
        TString Title_j = BranchNames[j];
        double Nominal_j, NominalError_j;
        GetNthParameter(j, Nominal_j, NominalError_j, Title_j);
      
        //OutputFile->cd();

        hpost2D[i][j]->SetMinimum(0);
        hpost2D[i][j]->GetXaxis()->SetTitle(Title_i);
        hpost2D[i][j]->GetYaxis()->SetTitle(Title_j);
        hpost2D[i][j]->GetZaxis()->SetTitle("Steps");
            
        for (int k = 0; k < nEntries; ++k) 
        {
            //KS: Burn in cut
            if(StepNumber[k] < BurnInCut) continue;

            //KS: Fill histogram with cached steps
            hpost2D[i][j]->Fill(ParStep[i][k], ParStep[j][k]);
        }
        // Get the Covariance for these two parameters
        (*Covariance)(i,j) = hpost2D[i][j]->GetCovariance();
        (*Covariance)(j,i) = (*Covariance)(i,j);

        (*Correlation)(i,j) = hpost2D[i][j]->GetCorrelationFactor();
        (*Correlation)(j,i) = (*Correlation)(i,j);

        /*
        if(printToPDF)
        {
            //KS: Skip Flux Params
            if(IsXsec[j] && IsXsec[i])
            {
                hpost2D->Draw("colz");
                Posterior->SetName(hpost2D->GetName());
                Posterior->SetTitle(hpost2D->GetTitle());
                Posterior->Print(CanvasName);
            }
        }
        */
        // Write it to root file
        //OutputFile->cd();
        //hpost2D->Write();

        //delete hpost2D;
    }// End j loop
}// End i loop

clock.Stop();
std::cout << "Making Covariance took " << clock.RealTime() << "s to finish for " << nEntries << " events" << std::endl;
    
  OutputFile->cd();
  Covariance->Write("Covariance");
  Correlation->Write("Correlation");
}

// *********************
// Make the covariance plots
void MCMCProcessor::DrawCovariance() {
// *********************
    
  int covBinning = nDraw;
  //If we only plot correlation/covariance between xsec (without flux)
  if(MakeOnlyXsecCorr)
  {
     covBinning = nParam[kXSecPar] - nFlux; 
  }
  
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
  for (int i = 0; i < covBinning; i++)
  {
    TString titlex = "";
    double nom, err;
    GetNthParameter(i, nom, err, titlex);
    
    hCov->GetXaxis()->SetBinLabel(i+1, titlex);
    hCovSq->GetXaxis()->SetBinLabel(i+1, titlex);
    hCorr->GetXaxis()->SetBinLabel(i+1, titlex);

    for (int j = 0; j < covBinning; j++) {

      // The value of the Covariance
      double cov = (*Covariance)(i,j);
      double corr = (*Correlation)(i,j);

      hCov->SetBinContent(i+1, j+1, cov);
      hCovSq->SetBinContent(i+1, j+1, ((cov > 0) - (cov < 0))*sqrt(fabs(cov)));
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
  
  delete hCov;
  delete hCovSq;
  delete hCorr;
}

// **************************
// Scan the input trees
void MCMCProcessor::ScanInput() {
  // **************************
  // Open the Chain
  Chain = new TChain("posteriors","posteriors");
  Chain->Add(MCMCFile.c_str());

  nEntries = Chain->GetEntries();

  // Get the list of branches
  TObjArray* brlis = (TObjArray*)(Chain->GetListOfBranches());

  // Get the number of branches
  nBranches = brlis->GetEntries();

  BranchNames.reserve(nBranches);
  IamVaried.reserve(nBranches);
  ParamType.reserve(nBranches);


  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  for (int i = 0; i < nBranches; i++) {

    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    // If we're on beam systematics
    if(bname.BeginsWith("xsec_")) 
      {
	BranchNames.push_back(bname);
	ParamType.push_back(kXSecPar);
	PlotXSec = true;
	nParam[kXSecPar]++;
      }
    if(!MakeOnlyXsecCorrFlux) //HI a bit of a dodgy way to do this, for now just want to get the work flow sorted 
      {
      if (bname.BeginsWith("ndd_") && PlotDet) 
	{
	  BranchNames.push_back(bname);
	  ParamType.push_back(kND280Par);
	  nParam[kND280Par]++;
	}
      else if (bname.BeginsWith("skd_joint_") && PlotDet)
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
    }
    //KS: as a bonus get LogL systeamtic
    else if (bname.BeginsWith("LogL_sample_")) {
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
  
  std::cout << "************************************************" << std::endl;
  std::cout << "Scanning output branches..." << std::endl;
  std::cout << "# useful entries in tree: " << nDraw  << std::endl;
  std::cout << "# XSec params:  " << nParam[kXSecPar] - nFlux <<" starting at "<<ParamTypeStartPos[kXSecPar] << std::endl;
  std::cout << "# Flux params:  " << nFlux << std::endl;
  std::cout << "# ND280 params: " << nParam[kND280Par] <<" starting at  "<<ParamTypeStartPos[kND280Par] << std::endl;
  std::cout << "# FD params:    " << nParam[kFDDetPar] <<" starting at  "<<ParamTypeStartPos[kFDDetPar] << std::endl;
  std::cout << "# Osc params:   " << nParam[kOSCPar]   <<" starting at  "<<ParamTypeStartPos[kOSCPar]   << std::endl;
  std::cout << "************************************************" << std::endl;


  // Set the step cut to be 20%
  int cut = Chain->GetMaximum("step")/5;
  SetStepCut(cut);
}

// ****************************
// Set up the output files and canvases
void MCMCProcessor::SetupOutput() {
  // ****************************

  // Make sure we can read files located anywhere and strip the .root ending
  MCMCFile = MCMCFile.substr(0, MCMCFile.find(".root"));

  // Check if the output file is ready
  if (OutputFile == NULL) MakeOutputFile();
  
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
    (*Central_Value)(i) = __UNDEF__;
    (*Means)(i) = __UNDEF__;
    (*Errors)(i) = __UNDEF__;
    (*Means_Gauss)(i) = __UNDEF__;
    (*Errors_Gauss)(i) = __UNDEF__;
    (*Means_HPD)(i) = __UNDEF__;
    (*Errors_HPD)(i) = __UNDEF__;
    (*Errors_HPD_Positive)(i) = __UNDEF__;
    (*Errors_HPD_Negative)(i) = __UNDEF__;
    for (int j = 0; j < nDraw; ++j) {
      (*Covariance)(i, j) = __UNDEF__;
      (*Correlation)(i, j) = __UNDEF__;
    }
  } 
  
  hpost = new TH1D*[nDraw]();

  OutputFile = NULL;
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
  
  if (OutputFile == NULL) MakeOutputFile();

  TH1D *PreFitPlot = new TH1D("Prefit", "Prefit", nDraw, 0, nDraw);
  for (int i = 0; i < PreFitPlot->GetNbinsX() + 1; ++i) {
    PreFitPlot->SetBinContent(i+1, 0);
    PreFitPlot->SetBinError(i+1, 0);
  }

  //KS: Sliglthy hacky way to get realtive to prior or nominal as this is convention we use,
  //Only applies for xsec, for other systematic it make no difference
  double CentralValueTemp, Central, Error;

  // Set labels and data
  for (int i = 0; i < nDraw; ++i)
  {
    //Those keep which parameter type we run currently and realtive number  
    int ParamNo = __UNDEF__;
    int ParamEnu = __UNDEF__;
    for(int ik = 0; ik < kNParameterEnum; ik++)
    {
        if (ParamType[i] == ParameterEnum(ik)) 
        {
            ParamNo = i - ParamTypeStartPos[ParameterEnum(ik)];
            ParamEnu = ParameterEnum(ik);
        }
    }
    CentralValueTemp = ParamCentral[ParamEnu][ParamNo];
    if(plotRelativeToPrior) 
    {
        // Normalise the prior relative the nominal/prior, just the way we get our fit results in MaCh3
        if ( CentralValueTemp != 0)
        {
            Central = ParamCentral[ParamEnu][ParamNo] / CentralValueTemp;
            Error = ParamErrors[ParamEnu][ParamNo]/CentralValueTemp;
        } else 
        {
            Central = CentralValueTemp + 1.0;
            Error = ParamErrors[ParamEnu][ParamNo];
        }
    }
    else
    {
        Central = CentralValueTemp;
        Error = ParamErrors[ParamEnu][ParamNo]; 
    }
    PreFitPlot->SetBinContent(i+1, Central);
    PreFitPlot->SetBinError(i+1, Error);
    PreFitPlot->GetXaxis()->SetBinLabel(i+1, ParamNames[ParamEnu][ParamNo]);
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
// Read the input Covariance matrix entries
// Get stuff like parameter input errors, names, and so on
void MCMCProcessor::ReadInputCov() {
  // **************************
  FindInputFiles();
  if(nParam[kXSecPar] > 0)  ReadXSecFile();
  if(nParam[kND280Par] > 0) ReadND280File();
  if(nParam[kFDDetPar] > 0) ReadFDFile();
  if(nParam[kOSCPar] > 0)   ReadOSCFile();
}

// **************************
// Read the output MCMC file and find what inputs were used
void MCMCProcessor::FindInputFiles() {
  // **************************

  // Now read the MCMC file
  TFile *TempFile = new TFile(MCMCFile.c_str(), "open");

  // Get the settings for the MCMC
  TTree *Settings = (TTree*)(TempFile->Get("Settings"));
  if (Settings == NULL) {
    std::cerr << "Didn't find Settings tree in MCMC file " << MCMCFile << std::endl;
    std::cerr << "Will try lowercase" << std::endl;
    TempFile->ls();
    Settings = (TTree*)(TempFile->Get("settings"));
    if (Settings == NULL) throw;
  }

  // Get the xsec Covariance matrix
  std::string *XSecInput = 0;
  if (Settings->SetBranchAddress("XsecCov", &XSecInput) < 0) {
    Settings->Print();
    std::cerr << "Couldn't find XsecCov branch in output" << std::endl;
    Settings->Print();
    throw;
  }

  // And the ND Covariance matrix
  std::string *ND280Input = 0;
  if (Settings->SetBranchAddress("NDCov", &ND280Input) < 0) {
    std::cerr << "Couldn't find NDCov branch in output" << std::endl;
    Settings->Print();
    throw;
  }

  // And the FD Covariance matrix
  std::string *FDInput = 0;
  if (Settings->SetBranchAddress("SKCov", &FDInput) < 0) {
    std::cerr << "Couldn't find SKCov branch in output" << std::endl;
    Settings->Print();
    throw;
  }
  
  // And the Osc Covariance matrix
  std::string *OscInput = 0;
  if (Settings->SetBranchAddress("oscCov", &OscInput) < 0) {
    std::cerr << "Couldn't find oscCov branch in output" << std::endl;
    Settings->Print();
    throw;
  }
  
  // Get the ND runs (needed for post fit distributions)
  std::string *NDrunsInput = 0;
  if (Settings->SetBranchAddress("NDruns", &NDrunsInput) < 0) {
    std::cerr << "Couldn't find NDruns branch in output" << std::endl;
    Settings->Print();
    throw;
  }

  // Get the vector of ND280 selections
  std::vector<std::string> *NDselInput = 0;
  if (Settings->SetBranchAddress("ND_Sel", &NDselInput) < 0) {
    std::cerr << "Couldn't find ND_Sel branch in output" << std::endl;
    Settings->Print();
    throw;
  }

  // Write the XSecCov and TempFile
  Settings->GetEntry(0);

  // Delete the TTrees and the input file handle since we've now got the settings we need
  delete Settings;

  // Delete the MCMCFile pointer we're reading
  TempFile->Close();
  delete TempFile;

  // Save the variables
  CovPos[kXSecPar]  = *XSecInput;
  CovPos[kND280Par] = *ND280Input;
  CovPos[kFDDetPar] = *FDInput;
  CovPos[kOSCPar]   = *OscInput;
  NDruns = *NDrunsInput;
  NDsel = *NDselInput;
}


// ***************
// Read the xsec file and get the input central values and errors
void MCMCProcessor::ReadXSecFile() {
  // ***************

  //KS:Most inputs are in ${MACH3}/inputs/blarb.root
  if (std::getenv("MACH3") != NULL) {
     std::cout << "Found MACH3 environment variable: " << std::getenv("MACH3") << std::endl;
      CovPos[kXSecPar].insert(0, std::string(std::getenv("MACH3"))+"/");
   }
   
  // Do the same for the cross-section
  TFile *XSecFile = new TFile(CovPos[kXSecPar].c_str(), "open");
  if (XSecFile->IsZombie()) {
    std::cerr << "Couldn't find XSecFile " << CovPos[kXSecPar] << std::endl;
    throw;
  }
  XSecFile->cd();

  // Get the matrix
  TMatrixDSym *XSecMatrix = (TMatrixDSym*)(XSecFile->Get("xsec_cov"));
  // Central priors
  TVectorD *XSecPrior = (TVectorD*)(XSecFile->Get("xsec_param_prior"));
  // Nominal values
  TVectorD *XSecNominal = (TVectorD*)(XSecFile->Get("xsec_param_nom"));
  // IDs
  TMatrixT<double> *XSecID = (TMatrixD*)(XSecFile->Get("xsec_param_id"));
  // Names
  TObjArray* xsec_param_names = (TObjArray*)(XSecFile->Get("xsec_param_names"));

  // Now make a TH1D of it
  ParamNames[kXSecPar].reserve(nParam[kXSecPar]);
  ParamCentral[kXSecPar].reserve(nParam[kXSecPar]);
  ParamNom[kXSecPar].reserve(nParam[kXSecPar]);
  ParamErrors[kXSecPar].reserve(nParam[kXSecPar]);

  for (int i = 0; i < nParam[kXSecPar]; ++i) {

    // Push back the name
    std::string TempString = std::string(((TObjString*)xsec_param_names->At(i))->GetString());
    ParamNames[kXSecPar].push_back(TempString);
    if(ParamNames[kXSecPar][i].BeginsWith("b_"))
    {
      IsXsec.push_back(false);
      nFlux++;
    } 
    else IsXsec.push_back(true);  
    
    ParamCentral[kXSecPar].push_back( ((*XSecPrior)(i)) );
    ParamNom[kXSecPar].push_back( ((*XSecNominal)(i)) );
    ParamErrors[kXSecPar].push_back( sqrt((*XSecMatrix)(i,i)) );
  }

  XSecFile->Close();
  delete XSecFile;
  delete XSecMatrix;
  delete XSecPrior;
  delete XSecNominal;
  delete XSecID;
}

// ***************
// Read the ND cov file and get the input central values and errors
void MCMCProcessor::ReadND280File() {
// ***************
  //KS:Most inputs are in ${MACH3}/inputs/blarb.root
  if (std::getenv("MACH3") != NULL) {
      std::cout << "Found MACH3 environment variable: " << std::getenv("MACH3") << std::endl;
      CovPos[kND280Par].insert(0, std::string(std::getenv("MACH3"))+"/");
   }
   
    // Do the same for the ND280
    TFile *NDdetFile = new TFile(CovPos[kND280Par].c_str(), "open");
    if (NDdetFile->IsZombie()) {
        std::cerr << "Couldn't find NDdetFile " << CovPos[kND280Par] << std::endl;
        throw;
    }
    NDdetFile->cd();
    
    TMatrixDSym *NDdetMatrix = (TMatrixDSym*)(NDdetFile->Get("nddet_cov"));
    TVectorD *NDdetNominal = (TVectorD*)(NDdetFile->Get("det_weights"));
    TObjArray* det_poly = (TObjArray*)(NDdetFile->Get("det_polys")->Clone());

    ParamNames[kND280Par].reserve(nParam[kND280Par]);
    ParamCentral[kND280Par].reserve(nParam[kND280Par]);
    ParamNom[kND280Par].reserve(nParam[kND280Par]);
    ParamErrors[kND280Par].reserve(nParam[kND280Par]);

    for (int i = 0; i < nParam[kND280Par]; ++i) 
    {
        ParamNom[kND280Par].push_back( (*NDdetNominal)(i) );
        ParamCentral[kND280Par].push_back( (*NDdetNominal)(i) );
        
        ParamErrors[kND280Par].push_back( sqrt((*NDdetMatrix)(i,i)) );
        ParamNames[kND280Par].push_back( Form("ND Det %i", i) );
    }  

    for (int j = 0; j <det_poly->GetLast()+1; ++j)
    {
        if( det_poly->At(j) != NULL)
        {
            TH2Poly *RefPoly = (TH2Poly*)det_poly->At(j);
            int size = RefPoly->GetNumberOfBins();
            NDSamplesBins.push_back(size);
            NDSamplesNames.push_back(RefPoly->GetTitle());
        }
    }

    NDdetFile->Close();
    delete NDdetFile;
    delete NDdetMatrix;
    delete NDdetNominal;
}


// ***************
// Read the FD cov file and get the input central values and errors
void MCMCProcessor::ReadFDFile() {
// ***************
  //KS:Most inputs are in ${MACH3}/inputs/blarb.root
  if (std::getenv("MACH3") != NULL) {
      std::cout << "Found MACH3 environment variable: " << std::getenv("MACH3") << std::endl;
      CovPos[kFDDetPar].insert(0, std::string(std::getenv("MACH3"))+"/");
   }
   
    // Do the same for the FD
    TFile *FDdetFile = new TFile(CovPos[kFDDetPar].c_str(), "open");
    if (FDdetFile->IsZombie()) {
        std::cerr << "Couldn't find FDdetFile " << CovPos[kFDDetPar] << std::endl;
        throw;
    }
    FDdetFile->cd();
    
    TMatrixDSym *FDdetMatrix = (TMatrixDSym*)(FDdetFile->Get("SKJointError_Erec_Total"));
    
    ParamNames[kFDDetPar].reserve(nParam[kFDDetPar]);
    ParamCentral[kFDDetPar].reserve(nParam[kFDDetPar]);
    ParamNom[kFDDetPar].reserve(nParam[kFDDetPar]);
    ParamErrors[kFDDetPar].reserve(nParam[kFDDetPar]);

    for (int i = 0; i < nParam[kFDDetPar]; ++i) 
    {
        //KS: FD parameters start at 1. in contrary to ND280
        ParamNom[kFDDetPar].push_back(1.);
        ParamCentral[kFDDetPar].push_back(1.);
        
        ParamErrors[kFDDetPar].push_back( sqrt((*FDdetMatrix)(i,i)) );
        ParamNames[kFDDetPar].push_back( Form("FD Det %i", i) );
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
  //KS:Most inputs are in ${MACH3}/inputs/blarb.root
  if (std::getenv("MACH3") != NULL) {
      std::cout << "Found MACH3 environment variable: " << std::getenv("MACH3") << std::endl;
      CovPos[kOSCPar].insert(0, std::string(std::getenv("MACH3"))+"/");
   }
   
    // Do the same for the ND280
    TFile *OscFile = new TFile(CovPos[kOSCPar].c_str(), "open");
    if (OscFile->IsZombie()) {
        std::cerr << "Couldn't find OSCFile " << CovPos[kOSCPar] << std::endl;
        throw;
    }
    OscFile->cd();
    
    TMatrixDSym *OscMatrix = (TMatrixDSym*)(OscFile->Get("osc_cov"));
    //KS: Osc nominal we can also set via config so there is danger that this will nor corrspond to what was used in the fit
    TVectorD *OscNominal = (TVectorD*)(OscFile->Get("osc_nom"));
    TObjArray* osc_param_names = (TObjArray*)(OscFile->Get("osc_param_names"));

    for (int i = 0; i < nParam[kOSCPar]; ++i) 
    {
        ParamNom[kOSCPar].push_back( (*OscNominal)(i) );
        ParamCentral[kOSCPar].push_back( (*OscNominal)(i) );
        
        ParamErrors[kOSCPar].push_back( sqrt((*OscMatrix)(i,i)) );
        // Push back the name
        std::string TempString = std::string(((TObjString*)osc_param_names->At(i))->GetString());
        ParamNames[kOSCPar].push_back(TempString);
    }  

    OscFile->Close();
    delete OscFile;
    delete OscMatrix;
    delete OscNominal;
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
void MCMCProcessor::SetStepCut(int Cuts) {
// ***************
  std::stringstream TempStream;
  TempStream << "step > " << Cuts;
  StepCut = TempStream.str();
  BurnInCut = Cuts;
}


// **************************
// Get the mean and RMS of a 1D posterior
void MCMCProcessor::GetArithmetic(TH1D * const hpost, int i) {
  // **************************
  (*Means)(i) = hpost->GetMean();
  (*Errors)(i) = hpost->GetRMS();
}

// **************************
// Get Gaussian characteristics
void MCMCProcessor::GetGaussian(TH1D *& hpost , int i) {
// **************************

  double mean = hpost->GetMean();
  double err = hpost->GetRMS();
  double peakval = hpost->GetBinCenter(hpost->GetMaximumBin());

  // Set the range for the Gaussian fit
  Gauss->SetRange(mean - 1.5*err , mean + 1.5*err);
  // Set the starting parameters close to RMS and peaks of the histograms
  Gauss->SetParameters(hpost->GetMaximum()*err*sqrt(2*3.14), peakval, err);

  // Perform the fit
  hpost->Fit(Gauss->GetName(),"Rq");
  hpost->SetStats(0);

  (*Means_Gauss)(i) = Gauss->GetParameter(1);
  (*Errors_Gauss)(i) = Gauss->GetParameter(2);
}


// ***************
// Get the highest posterior density from a TH1D
void MCMCProcessor::GetHPD(TH1D * const hpost, int i) {
// ***************
  // Get the bin which has the largest posterior density
  int MaxBin = hpost->GetMaximumBin();
  // And it's value
  double peakval = hpost->GetBinCenter(MaxBin);

  // The total integral of the posterior
  double integral = hpost->Integral();

  // Keep count of how much area we're covering
  double sum = 0.0;

  // Counter for current bin
  int CurrBin = MaxBin;
  while (sum/integral < 0.6827/2.0 && CurrBin < hpost->GetNbinsX()+1) {
    sum += hpost->GetBinContent(CurrBin);
    CurrBin++;
  }
  double sigma_p = fabs(hpost->GetBinCenter(MaxBin)-hpost->GetBinCenter(CurrBin));
  // Reset the sum
  sum = 0.0;

  // Reset the bin counter
  CurrBin = MaxBin;
  // Counter for current bin
  while (sum/integral < 0.6827/2.0 && CurrBin >= 0) {
    sum += hpost->GetBinContent(CurrBin);
    CurrBin--;
  }
  double sigma_m = fabs(hpost->GetBinCenter(CurrBin)-hpost->GetBinCenter(MaxBin));

  // Now do the double sided HPD
  sum = 0.0;
  int LowBin = MaxBin-1;
  int HighBin = MaxBin+1;
  double LowCon = 0.0;
  double HighCon = 0.0;

  while (sum/integral < 0.6827 && (LowBin >= 0 || HighBin < hpost->GetNbinsX()+1)) 
  {
    // Get the slice
    //KS:: If each slice reached histogram end then set value to 0, then other slice will be able to move further
    if(LowBin >= 0)LowCon = hpost->GetBinContent(LowBin);
    else LowCon = 0.0;
        
    if(HighBin < hpost->GetNbinsX()+1){HighCon = hpost->GetBinContent(HighBin);}
    else HighCon = 0.0;

    // If we're on the last slice and the lower contour is larger than the upper
    if ((sum+LowCon+HighCon)/integral > 0.6827 && LowCon > HighCon) {
      sum += LowCon;
      break;
      // If we're on the last slice and the upper contour is larger than the lower
    } else if ((sum+LowCon+HighCon)/integral > 0.6827 && HighCon >= LowCon) {
      sum += HighCon;
      break;
    } else {
      sum += LowCon + HighCon;
    }

    //KS:: Move further only if you haven't reached histogram end
    if(LowBin >= 0) LowBin--;
    if(HighBin < hpost->GetNbinsX()+1) HighBin++;
  }

  double sigma_hpd = 0.0;
  if (LowCon > HighCon) {
    sigma_hpd = fabs(hpost->GetBinCenter(LowBin)-hpost->GetBinCenter(MaxBin));
  } else {
    sigma_hpd = fabs(hpost->GetBinCenter(HighBin)-hpost->GetBinCenter(MaxBin));
  }

  (*Means_HPD)(i) = peakval;
  (*Errors_HPD)(i) = sigma_hpd;
  (*Errors_HPD_Positive)(i) = sigma_p;
  (*Errors_HPD_Negative)(i) = sigma_m;
}


// ***************
// Pass central value
void MCMCProcessor::GetNthParameter(int param, double &Nominal, double &NominalError, TString &Title){
// **************************
    for(int i = 0; i < kNParameterEnum; i++)
    {
        if (ParamType[param] == ParameterEnum(i)) 
        {
            int ParamNo = __UNDEF__;
            ParamNo = param - ParamTypeStartPos[ParameterEnum(i)];
            
            Nominal = ParamCentral[ParameterEnum(i)][ParamNo];
            NominalError = ParamErrors[ParameterEnum(i)][ParamNo];
            Title = ParamNames[ParameterEnum(i)][ParamNo];
        }
    }
}


// **************************************************
// Helper function to reset  histograms
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
// Diagnose the MCMC
void MCMCProcessor::DiagMCMC() {
  // **************************

// MCMC stuff to implement:
// Trace plots            -- DONE
// LogL vs step plots     -- DONE
// Acceptance probability -- DONE
// Autocorrelation        -- DONE
// _Batched Means_        -- DONE

    
  // Prepare branches etc for DiagMCMC
  PrepareDiagMCMC();

  // Draw the simple trace matrices
  ParamTraces();

  // Get the batched means
  BatchedMeans();

  // Draw the auto-correlations
  AutoCorrelation();

  // Draw acceptance Probability
  AcceptanceProbabilities();
}


// **************************
// Prepare branches etc. for DiagMCMC
void MCMCProcessor::PrepareDiagMCMC() {
  // **************************
  
   //bool doDiagMCMC = true;
    
    if(ParStep != NULL)
    {
        std::cout<<"It look like ParStep was already filled "<<std::endl;
        std::cout<<"Eventhough it is used for MakeCovariance_MP and for DiagMCMC "<<std::endl; 
        std::cout<<"it has differnt structure in both for cache hits, sorry "<<std::endl;
        throw;
    }
  // Initialise ParStep
  ParStep = new double*[nEntries]();
  SampleValues = new double*[nEntries]();
  SystValues = new double*[nEntries]();
  AccProbValues = new double[nEntries]();
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nEntries; ++i) {
    ParStep[i] = new double[nDraw]();
    SampleValues[i] = new double[nSamples]();
    SystValues[i] = new double[nSysts]();
    for (int j = 0; j < nDraw; ++j) {
      ParStep[i][j] = -999.99;
    }
    for (int j = 0; j < nSamples; ++j) {
      SampleValues[i][j] = -999.99;
    }
    for (int j = 0; j < nSysts; ++j) {
      SystValues[i][j] = -999.99;
    }
    AccProbValues[i] = -999.99;
  }

  // Initialise the sums
  ParamSums = new double[nDraw]();
  for (int i = 0; i < nDraw; ++i) {
    ParamSums[i] = 0.0;
  }
  
  std::cout << "Reading input tree..." << std::endl;
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
  
  // 10 entries output
  int countwidth = nEntries/10;

  // Can also do the batched means here to minimize excessive loops
  // The length of each batch
  int BatchLength = nEntries/nBatches+1;
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
  for (int i = 0; i < nEntries; ++i) {

    if (i % countwidth == 0) {
      std::cout << i << "/" << nEntries << " (" << double(i)/double(nEntries)*100. << "%)" << std::endl;
    }

    // Set the branch addresses for params
    for (int j = 0; j < nDraw; ++j) {
      Chain->SetBranchAddress(BranchNames[j].Data(), &ParStep[i][j]);
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
      ParamSums[j] += ParStep[i][j];
      BatchedAverages[BatchNumber][j] += ParStep[i][j];
    }
    
      //KS: Could easyli add this to above loop but I accProb is different beast so better keep it like this
      AccProbBatchedAverages[BatchNumber] += AccProbValues[i];
  }

  clock.Stop();

  std::cout << "Took " << clock.RealTime() << "s to finish cachin statysitc for Diag MCMC with " << nEntries << " events" << std::endl;

  // Make the sums into average
  for (int i = 0; i < nDraw; ++i) {
    ParamSums[i] /= nEntries;
    for (int j = 0; j < nBatches; ++j) {
      // Divide by the total number of events in the batch
      BatchedAverages[j][i] /= BatchLength;
      if(i==0) AccProbBatchedAverages[j] /= BatchLength; //KS: we have only one accProb, keep it like this for now
    }
  }

  // And make our sweet output file
  if (OutputFile == NULL) MakeOutputFile();
  
}


// *****************
// Draw trace plots of the parameters
// i.e. parameter vs step
void MCMCProcessor::ParamTraces() {
  // *****************

  std::cout << "Making trace plots..." << std::endl;

  // Make the TH1Ds
  TraceParamPlots = new TH1D*[nDraw];
  TraceSamplePlots = new TH1D*[nSamples];
  TraceSystsPlots = new TH1D*[nSysts];

  // Set the titles and limits for TH2Ds
  for (int j = 0; j < nDraw; ++j) {

    TString Title = BranchNames[j];
    double Nominal = 1.0;
    double NominalError = 1.0;
    
    GetNthParameter(j, Nominal, NominalError, Title);
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

  // Have now made the empty TH2Ds, now for writing content to them!

  // Loop over the number of parameters to draw their traces
  // Each histogram
#ifdef MULTITHREAD
  std::cout << "Using multi-threading..." << std::endl;
#pragma omp parallel for
#endif
  for (int i = 0; i < nEntries; ++i) {
    // Set bin content for the ith bin to the parameter values
    for (int j = 0; j < nDraw; ++j) {
      TraceParamPlots[j]->SetBinContent(i, ParStep[i][j]);
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
    delete SampleValues[j];
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
}


// *********************************
void MCMCProcessor::AutoCorrelation() {
  // *********************************

  std::cout << "Making auto-correlations..." << std::endl;
  TStopwatch clock;
  clock.Start();
  const int nLags = 25000;

  // The sum of (Y-Ymean)^2 over all steps for each parameter
  double **DenomSum = new double*[nDraw]();
  double **NumeratorSum = new double*[nDraw]();
  for (int i = 0; i < nDraw; ++i) {
    DenomSum[i] = new double[nLags];
    NumeratorSum[i] = new double[nLags];
  }
  LagKPlots = new TH1D*[nDraw];

  // Loop over the parameters of interacts
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int j = 0; j < nDraw; ++j) {

    // Loop over each lag
    for (int k = 0; k < nLags; ++k) {
      NumeratorSum[j][k] = 0.0;
      DenomSum[j][k] = 0.0;
    }

    // Make TH1Ds for each parameter which hold the lag
    TString Title = BranchNames[j];
    double Nominal = 1.0;
    double NominalError = 1.0;
    
    GetNthParameter(j, Nominal, NominalError, Title);
    std::string HistName = Form("%s_%s_Lag", Title.Data(), BranchNames[j].Data());
    LagKPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nLags, 0.0, nLags);
    LagKPlots[j]->GetXaxis()->SetTitle("Lag");
    LagKPlots[j]->GetYaxis()->SetTitle("Auto-correlation function");
  }

  // Loop over the lags
  // Each lag is indepdent so might as well multi-thread them!
#ifdef MULTITHREAD
  std::cout << "Using multi-threading..." << std::endl;
//KS: Consider moving it to GPU, and have huge loop over nLag*nDraw, should work
#pragma omp parallel for
#endif
  for (int k = 0; k < nLags; ++k) {

    // Loop over the number of entries
    for (int i = 0; i < nEntries; ++i) {

      // Loop over the number of parameters
      for (int j = 0; j < nDraw; ++j) {

        double Diff = ParStep[i][j]-ParamSums[j];

        // Only sum the numerator up to i = N-k
        if (i < nEntries-k) {
          double LagTerm = ParStep[i+k][j]-ParamSums[j];
          double Product = Diff*LagTerm;
          NumeratorSum[j][k] += Product;
        }

        // Square the difference to form the denominator
        double Denom = Diff*Diff;
        DenomSum[j][k] += Denom;
      }
    }
  }

  OutputFile->cd();
  TDirectory *AutoCorrDir = OutputFile->mkdir("Auto_corr");
  // Now fill the LagK auto-correlation plots
  for (int j = 0; j < nDraw; ++j) {
    for (int k = 0; k < nLags; ++k) {
      LagKPlots[j]->SetBinContent(k, NumeratorSum[j][k]/DenomSum[j][k]);
    }
    AutoCorrDir->cd();
    LagKPlots[j]->Write();
    delete LagKPlots[j];
  }
  delete[] LagKPlots;
  for (int i = 0; i < nEntries; ++i) {
    delete ParStep[i];
  }
  delete[] ParStep;
  
  delete ParamSums;

  clock.Stop();
  std::cout << "It took " << clock.RealTime() << std::endl;
}


// **************************
// Batched means, literally read from an array and chuck into TH1D
void MCMCProcessor::BatchedMeans() {
  // **************************

  BatchedParamPlots = new TH1D*[nDraw];
  for (int j = 0; j < nDraw; ++j) {
    TString Title = BranchNames[j];
    double Nominal = 1.0;
    double NominalError = 1.0;
    
    GetNthParameter(j, Nominal, NominalError, Title);
    
    std::string HistName = Form("%s_%s_batch", Title.Data(), BranchNames[j].Data());
    BatchedParamPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nBatches, 0, nBatches);
  }

  for (int i = 0; i < nBatches; ++i) {
    for (int j = 0; j < nDraw; ++j) {
      BatchedParamPlots[j]->SetBinContent(i+1, BatchedAverages[i][j]);
      int BatchRangeLow = double(i)*double(nEntries)/double(nBatches);
      int BatchRangeHigh = double(i+1)*double(nEntries)/double(nBatches);
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

  for (int i = 0; i < nBatches; ++i) {
    delete BatchedAverages[i];
  }

  delete[] BatchedAverages;

}

// **************************
// Acceptance Probability
void MCMCProcessor::AcceptanceProbabilities() {
  // **************************
    std::cout << "Making AccProb plots..." << std::endl;

    // Set the titles and limits for TH2Ds

    AcceptanceProbPlot = new TH1D("AcceptanceProbability", "Acceptance Probability", nEntries, 0, nEntries);
    AcceptanceProbPlot->GetXaxis()->SetTitle("Step");
    AcceptanceProbPlot->GetYaxis()->SetTitle("Acceptance Probability");

    BatchedAcceptanceProblot = new TH1D("AcceptanceProbability_Batch", "AcceptanceProbability_Batch", nBatches, 0, nBatches);
    BatchedAcceptanceProblot->GetYaxis()->SetTitle("Acceptance Probability");
    
  for (int i = 0; i < nBatches; ++i) {
      BatchedAcceptanceProblot->SetBinContent(i+1, AccProbBatchedAverages[i]);
      int BatchRangeLow = double(i)*double(nEntries)/double(nBatches);
      int BatchRangeHigh = double(i+1)*double(nEntries)/double(nBatches);
      std::stringstream ss;
      ss << BatchRangeLow << " - " << BatchRangeHigh;
      BatchedAcceptanceProblot->GetXaxis()->SetBinLabel(i+1, ss.str().c_str());
  }
  
  
#ifdef MULTITHREAD
  std::cout << "Using multi-threading..." << std::endl;
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
}

