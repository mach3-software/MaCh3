#include "MCMCProcessor.h"

MCMCProcessor::MCMCProcessor(const std::string &InputFile, bool MakePostfitCorr) : 
  Chain(NULL), StepCut(""), MakeCorr(MakePostfitCorr), MadePostfit(false) {

  MCMCFile = InputFile;

  std::cout << "Making post-fit processor for " << MCMCFile << std::endl;

  // Read the input Covariances
  ReadInputCov();

  // Scan the ROOT file for useful branches
  ScanInput();

  // Setup the output
  SetupOutput();
}

// ****************************
// The destructor
MCMCProcessor::~MCMCProcessor() {
// ****************************

  // Close the pdf file
  std::cout << "Closing pdf in MCMCProcessor " << CanvasName << std::endl;
  CanvasName += "]";
  //Posterior->Print(CanvasName);

  if (OutputFile != NULL) OutputFile->Close();
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
  Peak_Values = Peaks;
}

// ***************
// Get post-fits for the ParameterEnum type, e.g. xsec params, Near Det params or flux params
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
    (*Peak_Values)(ParamNumber) = (*Peaks)(i);
    ++ParamNumber;
  }
}


// ***************
void MCMCProcessor::GetCovariance(TMatrixDSym *&Cov, TMatrixDSym *&Corr) {
// ***************
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
  OutputName = MCMCFile + "_Process.root";

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
  nBins = 50;
  DrawRange = 1.3;

  // nDraw is number of draws we want to do
  for (int i = 0; i < nDraw; ++i) {

    if (i % (nDraw/5) == 0) {
      std::cout << "  " << i << "/" << nDraw << " (" << int((double(i)/double(nDraw)*100.0))+1 << "%)" << std::endl;
    }

    OutputFile->cd();
    TString Title = BranchNames[i];
    double Nominal = 1.0;

    if (ParamType[i] == kXSecPar) {
      int ParamNo = __UNDEF__;
      if (PlotFlux) {
        ParamNo = i - nFlux;
      } 
      Nominal = XSecCentral[ParamNo];
      Title = XSecNames[ParamNo];
    }

    OutputFile->cd();
    // This holds the posterior density
    TH1D *hpost = new TH1D(BranchNames[i], BranchNames[i], nBins, 0., 2.);
    hpost->SetMinimum(0);
    hpost->GetYaxis()->SetTitle("Steps");
    hpost->GetYaxis()->SetNoExponent(false);
    double maximum = Chain->GetMaximum(BranchNames[i]);
    double minimum = Chain->GetMinimum(BranchNames[i]);
    hpost->SetBins(nBins, minimum, maximum);
    hpost->SetTitle(BranchNames[i]);

    // Project BranchNames[i] onto hpost, applying stepcut
    Chain->Project(BranchNames[i], BranchNames[i], StepCut.c_str());

    // Get the characteristics of the histogram
    double mean = hpost->GetMean();
    double rms = hpost->GetRMS();
    double peakval = hpost->GetBinCenter(hpost->GetMaximumBin());

    // Set the range for the Gaussian fit
    Gauss->SetRange(mean - 1.5*rms , mean + 1.5*rms);
    // Set the starting parameters close to RMS and peaks of the histograms
    Gauss->SetParameters(hpost->GetMaximum()*rms*sqrt(2*3.14), peakval, rms);

    // Perform the fit
    hpost->Fit("Gauss","Rq");
    hpost->SetStats(0);

    // Write the results from the projection into the TVectors and TMatrices
    (*Means)(i) = mean;
    (*Errors)(i) = rms;
    (*Means_Gauss)(i) = Gauss->GetParameter(1);
    (*Errors_Gauss)(i) = Gauss->GetParameter(2);
    (*Peaks)(i) = peakval;
    (*Covariance)(i,i) = rms*rms;
    (*Correlation)(i,i) = 1.0;

    hpost->SetLineWidth(2);
    hpost->SetMaximum(hpost->GetMaximum()*DrawRange);
    hpost->SetTitle(Title);
    hpost->GetXaxis()->SetTitle(hpost->GetTitle());
    Gauss->SetLineWidth(2);

    // Now make the TLine for the asimov
    TLine *Asimov = new TLine(Nominal, hpost->GetMinimum(), Nominal, hpost->GetMaximum());
    Asimov->SetLineColor(kBlack);
    Asimov->SetLineWidth(2);
    Asimov->SetLineStyle(kDashed);

    // Make the legend
    TLegend *leg = new TLegend(0.12, 0.6, 0.6, 0.97);
    leg->SetTextSize(0.04);
    leg->AddEntry(hpost, Form("#splitline{PDF}{#mu = %.2f, #sigma = %.2f}", hpost->GetMean(), hpost->GetRMS()), "l");
    leg->AddEntry(Gauss, Form("#splitline{Gauss}{#mu = %.2f, #sigma = %.2f}", Gauss->GetParameter(1), Gauss->GetParameter(2)), "l");
    leg->AddEntry(Asimov, Form("#splitline{Asimov}{x = %.2f}", Nominal), "l");
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);

    // Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
    if (hpost->GetMaximum() == hpost->Integral()*DrawRange) {
      std::cout << "Found fixed parameter, moving on" << std::endl;
      IamVaried[i] = false;
      delete hpost;
      delete Asimov;
      delete leg;
      continue;
    }

    // Store that this parameter is indeed being varied
    IamVaried[i] = true;

    // Write to file
    Posterior->SetName(hpost->GetName());
    Posterior->SetTitle(hpost->GetTitle());
    //Posterior->Print(CanvasName);

    // Draw onto the TCanvas
    hpost->Draw();
    Asimov->Draw("same");
    leg->Draw("same");

    // cd into params directory in root file
    PostDir->cd();
    Posterior->Write();

    delete hpost;
    delete Asimov;
    delete leg;
  } // end the for loop over nParams

  OutputFile->cd();
  TDirectory *Names = OutputFile->mkdir("Names");
  Names->cd();
  for (std::vector<TString>::iterator it = BranchNames.begin(); it != BranchNames.end(); ++it) {
    TObjString((*it)).Write();
  }
  OutputFile->cd();
  Means->Write("PDF_Means");
  Errors->Write("PDF_Error");
  Means_Gauss->Write("Gauss_Means");
  Errors_Gauss->Write("Gauss_Errors");
  Peaks->Write("Peaks");

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
  paramPlot->SetFillStyle(3144);
  paramPlot->SetFillColor(kBlue-3);
  paramPlot->SetMarkerColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerStyle(20);
  paramPlot->SetLineColor(paramPlot->GetFillColor());
  paramPlot->SetMarkerSize(prefit->GetMarkerSize());

  // Same but with Gaussian output
  TH1D *paramPlot_Gauss = (TH1D*)(paramPlot->Clone());
  paramPlot_Gauss->SetMarkerColor(kYellow-7);
  paramPlot_Gauss->SetMarkerStyle(21);
  paramPlot_Gauss->SetLineWidth(2);
  paramPlot_Gauss->SetMarkerSize(1.3);
  paramPlot_Gauss->SetFillColor(0);
  paramPlot_Gauss->SetFillStyle(0);
  paramPlot_Gauss->SetLineColor(paramPlot_Gauss->GetMarkerColor());

  // Set labels and data
  for (int i = 0; i < nDraw; ++i) {

    paramPlot->SetBinContent(i+1, (*Means)(i));
    paramPlot->SetBinError(i+1, (*Errors)(i));

    paramPlot_Gauss->SetBinContent(i+1, (*Means_Gauss)(i));
    paramPlot_Gauss->SetBinError(i+1, (*Errors_Gauss)(i));

    paramPlot->GetXaxis()->SetBinLabel(i+1, prefit->GetXaxis()->GetBinLabel(i+1));
    paramPlot_Gauss->GetXaxis()->SetBinLabel(i+1, prefit->GetXaxis()->GetBinLabel(i+1));
  }

  // Make a TLegend
  TLegend *CompLeg = new TLegend(0.33, 0.73, 0.76, 0.95);
  CompLeg->AddEntry(prefit, "Prefit", "fp");
  CompLeg->AddEntry(paramPlot, "Postfit PDF", "fp");
  CompLeg->AddEntry(paramPlot_Gauss, "Postfit Gauss", "lfep");
  CompLeg->SetFillColor(0);
  CompLeg->SetFillStyle(0);
  CompLeg->SetLineWidth(0);
  CompLeg->SetLineStyle(0);
  CompLeg->SetBorderSize(0);

  OutputFile->cd();
  // Plot the flux parameters (0 to 100) if enabled
  // Have already looked through the branches earlier
  if (PlotFlux == true) {
    prefit->GetYaxis()->SetTitle("Variation");
    prefit->GetYaxis()->SetRangeUser(0.7, 1.3);
    prefit->GetXaxis()->SetRangeUser(0, nFlux);
    prefit->GetXaxis()->SetTitle("");
    prefit->GetXaxis()->LabelsOption("v");

    paramPlot->GetYaxis()->SetTitle("Variation");
    paramPlot->GetYaxis()->SetRangeUser(0.7, 1.3);
    paramPlot->GetXaxis()->SetRangeUser(0, nFlux);
    paramPlot->GetXaxis()->SetTitle("");
    paramPlot->SetTitle(StepCut.c_str());
    paramPlot->GetXaxis()->LabelsOption("v");

    paramPlot_Gauss->GetYaxis()->SetTitle("Variation");
    paramPlot_Gauss->GetYaxis()->SetRangeUser(0.7, 1.3);
    paramPlot_Gauss->GetXaxis()->SetRangeUser(0, nFlux);
    paramPlot_Gauss->GetXaxis()->SetTitle("");
    paramPlot_Gauss->SetTitle(StepCut.c_str());
    paramPlot_Gauss->GetXaxis()->LabelsOption("v");

    prefit->Write("param_flux_prefit");
    paramPlot->Write("param_flux");
    paramPlot_Gauss->Write("param_flux_gaus");

    prefit->Draw("e2");
    paramPlot->Draw("e2, same");
    paramPlot_Gauss->Draw("e1, same");

    CompLeg->Draw("same");
    Posterior->Write("param_flux_canv");
    Posterior->Clear();
  }

  Posterior->SetBottomMargin(0.2);
  // Plot the xsec parameters (100 to ~125)
  // Have already looked through the branches earlier
  OutputFile->cd();
  if (PlotXSec == true) {
    prefit->GetYaxis()->SetTitle("Variation relative prefit");
    prefit->GetYaxis()->SetRangeUser(-0.0, 2.0);
    prefit->GetXaxis()->SetTitle("");
    prefit->GetXaxis()->LabelsOption("v");

    prefit->GetXaxis()->SetRangeUser(nFlux, nFlux+nXSec);
    paramPlot->GetXaxis()->SetRangeUser(nFlux, nFlux+nXSec);
    paramPlot_Gauss->GetXaxis()->SetRangeUser(nFlux, nFlux+nXSec);

    // Write the individual ones
    prefit->Write("param_xsec_prefit");
    paramPlot->Write("param_xsec");
    paramPlot_Gauss->Write("param_xsec_gaus");

    // And the combined
    prefit->Draw("e2");
    paramPlot->Draw("e2, same");
    paramPlot_Gauss->Draw("same");
    CompLeg->Draw("same");
    Posterior->Write("param_xsec_canv");
    Posterior->Clear();
  }

  delete CompLeg;
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

  // Now we are sure we have the diagonal elements, let's make the off-diagonals
  for (int i = 0; i < nDraw; ++i) {

    if (i % (nDraw/5) == 0) {
      std::cout << "  " << i << "/" << nDraw << " (" << int((double(i)/double(nDraw)*100.0))+1 << "%)" << std::endl;
    }

    TString Title_i = BranchNames[i];

    // For xsec parameters get the parameter name
    if (ParamType[i] == kXSecPar) {
      int ParamNo = 0;
      if (PlotFlux && i >= nFlux) {
        ParamNo = i - nFlux;
      }
      Title_i = XSecNames[ParamNo];
    }

    double min_i = Chain->GetMinimum(BranchNames[i]);
    double max_i = Chain->GetMaximum(BranchNames[i]);

    // Loop over the other parameters to get the correlations
    for (int j = 0; j <= i; ++j) {

      // Skip the diagonal elements which we've already done above
      if (j == i) continue;

      // If this parameter isn't varied
      if (IamVaried[j] == false) {
        (*Covariance)(i,j) = 0.0;
        (*Correlation)(i,j) = 0.0;
        continue;
      }

      TString Title_j = BranchNames[j];

      // For xsec parameters get the parameter name
      if (ParamType[j] == kXSecPar) {
        int ParamNo = 0;
        if (PlotFlux && j >= nFlux) {
          ParamNo = j - nFlux;
        }
        Title_j = XSecNames[ParamNo];
      }

      OutputFile->cd();

      // The draw which we want to perform
      TString DrawMe = BranchNames[j]+":"+BranchNames[i];

      double max_j = Chain->GetMaximum(BranchNames[j]);
      double min_j = Chain->GetMinimum(BranchNames[j]);

      // TH2F to hold the Correlation 
      TH2D *hpost = new TH2D(DrawMe, DrawMe, nBins, min_i, max_i, nBins, min_j, max_j);

      hpost->SetMinimum(0);
      hpost->GetXaxis()->SetTitle(Title_i);
      hpost->GetYaxis()->SetTitle(Title_j);
      hpost->GetZaxis()->SetTitle("Steps");

      // The draw command we want, i.e. draw param j vs param i
      Chain->Project(DrawMe, DrawMe, StepCut.c_str());

      // Get the Covariance for these two parameters
      (*Covariance)(i,j) = hpost->GetCovariance();
      (*Covariance)(j,i) = (*Covariance)(i,j);

      (*Correlation)(i,j) = hpost->GetCorrelationFactor();
      (*Correlation)(j,i) = (*Correlation)(i,j);

      Posterior->SetName(hpost->GetName());
      Posterior->SetTitle(hpost->GetTitle());
      //Posterior->Print(CanvasName);

      // Write it to root file
      //OutputFile->cd();
      //hpost->Write();

      delete hpost;
    } // End j loop
  } // End i loop
  OutputFile->cd();
  Covariance->Write("Covariance");
  Correlation->Write("Correlation");
}

// *********************
// Make the covariance plots
void MCMCProcessor::DrawCovariance() {
// *********************
  // The Covariance matrix from the fit
  TH2D* hCov = new TH2D("hCov", "hCov", nDraw, 0, nDraw, nDraw, 0, nDraw);
  hCov->GetZaxis()->SetTitle("Covariance");
  // The Covariance matrix square root, with correct sign
  TH2D* hCovSq = new TH2D("hCovSq", "hCovSq", nDraw, 0, nDraw, nDraw, 0, nDraw);
  hCovSq->GetZaxis()->SetTitle("Covariance");
  // The Correlation
  TH2D* hCorr = new TH2D("hCorr", "hCorr", nDraw, 0, nDraw, nDraw, 0, nDraw);
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
  for (int i = 0; i < nDraw; i++) {

    TString titlex = "";
    if (ParamType[i] == kFluxPar){
      titlex = FluxNames[i];
    } else if (ParamType[i] == kXSecPar) {
      titlex = XSecNames[i-nFlux];
    } else if (ParamType[i] == kNearDetPar) {
      titlex = "";
    }

    hCov->GetXaxis()->SetBinLabel(i+1, titlex);
    hCovSq->GetXaxis()->SetBinLabel(i+1, titlex);
    hCorr->GetXaxis()->SetBinLabel(i+1, titlex);

    for (int j = 0; j < nDraw; j++) {

      // The value of the Covariance
      double cov = (*Covariance)(i,j);
      double corr = (*Correlation)(i,j);

      hCov->SetBinContent(i+1, j+1, cov);
      hCovSq->SetBinContent(i+1, j+1, ((cov > 0) - (cov < 0))*sqrt(fabs(cov)));
      hCorr->SetBinContent(i+1, j+1, corr);

      TString titley = "";
      if (ParamType[j] == kFluxPar){
        titley = FluxNames[j];
      } else if (ParamType[j] == kXSecPar) {
        titley = XSecNames[j-nFlux];
      } else if (ParamType[j] == kNearDetPar) {
        titley = "";
      }

      hCov->GetYaxis()->SetBinLabel(j+1, titley);
      hCovSq->GetYaxis()->SetBinLabel(j+1, titley);
      hCorr->GetYaxis()->SetBinLabel(j+1, titley);
    }

  }

  // Take away the stat box
  gStyle->SetOptStat(0);
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
  hCov->Draw("colz");
  Posterior->SetRightMargin(0.15);
  //Posterior->Print(CanvasName);

  Posterior->cd();
  Posterior->Clear();
  hCorr->Draw("colz");
  Posterior->SetRightMargin(0.15);
  //Posterior->Print(CanvasName);

  hCov->Write("Covariance_plot");
  hCovSq->Write("Covariance_sq_plot");
  hCorr->Write("Correlation_plot");
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

  // Number of flux, cross-section and ND280 parameters
  nFlux = 0;
  nXSec = 0;
  nNear = 0;

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  for (int i = 0; i < nBranches; i++) {

    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    // If we're on beam systematics
    if (bname.BeginsWith("b_")) {
      BranchNames.push_back(bname);
      ParamType.push_back(kFluxPar);
      PlotFlux = true;
      nFlux++;
    } else if(bname.BeginsWith("xsec_")) {
      BranchNames.push_back(bname);
      ParamType.push_back(kXSecPar);
      PlotXSec = true;
      nXSec++;
    } else if (bname.BeginsWith("ndd_")) {
      BranchNames.push_back(bname);
      ParamType.push_back(kNearDetPar);
      PlotDet = true;
      nNear++;
    }
  }

  nDraw = BranchNames.size();

  std::cout << "************************************************" << std::endl;
  std::cout << "Scanning output branches..." << std::endl;
  std::cout << "# useful entries in tree: " << nDraw << std::endl;
  std::cout << "# Flux params: " << nFlux << std::endl;
  std::cout << "# XSec params: " << nXSec << std::endl;
  std::cout << "# ND280 params: " << nNear << std::endl;
  std::cout << "************************************************" << std::endl;

  // Also read what input Covariances were used
  ReadInputCov();

  // Set the step cut to be 25%
  SetStepCut(nEntries/4);
}

// ****************************
// Set up the output files and canvases
void MCMCProcessor::SetupOutput() {
  // ****************************

  // Make sure we can read files located anywhere and strip the .root ending
  MCMCFile = MCMCFile.substr(0, MCMCFile.find(".root"));

  CanvasName = MCMCFile + "_Process.pdf[";
  //Posterior->Print(CanvasName);

  // Once the pdf file is open no longer need to bracket
  CanvasName.ReplaceAll("[","");

  // We fit with this Gaussian
  Gauss = new TF1("Gauss", "[0]/(sqrt(2.0*3.14159)*[2])*TMath::Exp(-0.5*pow(x-[1],2)/pow([2],2))", -2.0, 3.0);

  // Declare the TVectors
  Covariance = new TMatrixDSym(nDraw);
  Correlation = new TMatrixDSym(nDraw);
  Means = new TVectorD(nDraw);
  Errors = new TVectorD(nDraw);
  Means_Gauss = new TVectorD(nDraw);
  Errors_Gauss = new TVectorD(nDraw);
  Peaks = new TVectorD(nDraw);

  // Initialise to something silly
  for (int i = 0; i < nDraw; ++i) {
    (*Means)(i) = __UNDEF__;
    (*Errors)(i) = __UNDEF__;
    (*Means_Gauss)(i) = __UNDEF__;
    (*Errors_Gauss)(i) = __UNDEF__;
    (*Peaks)(i) = __UNDEF__;
    for (int j = 0; j < nDraw; ++j) {
      (*Covariance)(i, j) = __UNDEF__;
      (*Correlation)(i, j) = __UNDEF__;
    }
  }

  OutputFile = NULL;
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

  // Set the flux values
  for (int i = 0; i < nFlux; ++i) {
    PreFitPlot->SetBinContent(i+1, FluxCentral[i]);
    PreFitPlot->SetBinError(i+1, FluxErrors[i]);

    // Set title for every fifth parameter
    if (i % 5 == 0) {
      PreFitPlot->GetXaxis()->SetBinLabel(i+1, Form("Flux param %i", i));
    }
  }

  // Set the xsec values
  for (int i = nFlux; i < nFlux + nXSec; ++i) {
    PreFitPlot->SetBinContent(i+1, XSecCentral[i-nFlux]);
    PreFitPlot->SetBinError(i+1, XSecErrors[i-nFlux]);
    PreFitPlot->GetXaxis()->SetBinLabel(i+1, XSecNames[i-nFlux]);
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
  ReadFluxFile();
  ReadXSecFile();
  ReadNearFile();
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

  // And the flux Covariance matrix
  std::string *FluxInput = 0;
  if (Settings->SetBranchAddress("FluxCov", &FluxInput) < 0) {
    std::cerr << "Couldn't find FluxCov branch in output" << std::endl;
    Settings->Print();
    throw;
  }

  // And the flux Covariance matrix
  std::string *NearInput = 0;
  if (Settings->SetBranchAddress("NDCov", &NearInput) < 0) {
    std::cerr << "Couldn't find NDCov branch in output" << std::endl;
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

  // Get the vector of Near selections
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
  FluxCov = *FluxInput;
  XSecCov = *XSecInput;
  NearCov = *NearInput;
  NDruns = *NDrunsInput;
  NDsel = *NDselInput;
}


// ***************
// Read the flux file and get the input central values and errors
void MCMCProcessor::ReadFluxFile() {
  // ***************
  // Get the flux Covariance matrix
  TFile *FluxFile = new TFile(FluxCov.c_str(), "open");
  if (FluxFile->IsZombie()) {
    std::cerr << "Couldn't find FluxFile " << FluxCov << std::endl;
    throw;
  }
  FluxFile->cd();

  TMatrixDSym *FluxMatrix = (TMatrixDSym*)(FluxFile->Get("total_flux_cov"));
  int nFlux = FluxMatrix->GetNrows();

  FluxCentral.reserve(nFlux);
  FluxErrors.reserve(nFlux);

  for (int i = 0; i < nFlux; ++i) {
    FluxCentral.push_back(1.0);
    FluxErrors.push_back((*FluxMatrix)(i,i));
    FluxNames.push_back(Form("Flux %i", i));
  }

  FluxFile->Close();
  delete FluxFile;
  delete FluxMatrix;
}

// ***************
// Read the xsec file and get the input central values and errors
void MCMCProcessor::ReadXSecFile() {
  // ***************

  // Do the same for the cross-section
  TFile *XSecFile = new TFile(XSecCov.c_str(), "open");
  if (XSecFile->IsZombie()) {
    std::cerr << "Couldn't find XSecFile " << XSecCov << std::endl;
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
  int nXSec = XSecPrior->GetNrows();
  XSecCentral.reserve(nXSec);
  XSecErrors.reserve(nXSec);

  for (int i = 0; i < nXSec; ++i) {

    double central = __UNDEF__;
    double error = __UNDEF__;

    // Push back the name
    std::string TempString = std::string(((TObjString*)xsec_param_names->At(i))->GetString());
    XSecNames.push_back(TempString);

    // Normalise the prior relative the nominal, just the way we get our fit results in MaCh3
    if (((*XSecNominal)(i)) != 0 && ((*XSecID)(i,0)) != -2) {
      central = ((*XSecPrior)(i)) / ((*XSecNominal)(i));
      error = sqrt((*XSecMatrix)(i,i))/((*XSecNominal)(i));
      // If the nominal is zero or a function we set it to the prior
    } else {
      central = ((*XSecPrior)(i));
      error = sqrt((*XSecMatrix)(i,i));
    }

    XSecCentral.push_back(central);
    XSecErrors.push_back(error);
  }

  XSecFile->Close();
  delete XSecMatrix;
  delete XSecPrior;
  delete XSecNominal;
  delete XSecID;
}

// ***************
// TODO SET NEAR DET PARS
void MCMCProcessor::ReadNearFile() {
// ***************
}


// ***************
// Make the step cut from a string
void MCMCProcessor::SetStepCut(std::string Cuts) {
// ***************
  StepCut = Cuts;
}

// ***************
// Make the step cut from an int
void MCMCProcessor::SetStepCut(int Cuts) {
// ***************
  std::stringstream TempStream;
  TempStream << "step > " << Cuts;
  StepCut = TempStream.str();
}

