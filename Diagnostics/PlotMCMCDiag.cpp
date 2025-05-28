// MaCh3 includes
#include "Fitters/MCMCProcessor.h"
#include "Samples/HistogramUtils.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"

/// @file PlotMCMCDiag.cpp
/// @brief KS: This script is used to analyse output form DiagMCMC.
/// @warning This script support comparing up to 4 files, there is easy way to expand it up to five or six,
/// @todo this need serious refactor

TString DUMMYFILE = "KillMePlease";

/// @brief KS: function which looks for minimum in given range
double GetMinimumInRange(TH1D *hist, double minRange, double maxRange)
{
  double MinVale = 1234567890.;
  TAxis* xaxis = hist->GetXaxis();
  for(int x = 1; x <= hist->GetNbinsX(); x++)
  {
    if(xaxis->GetBinLowEdge(x)>minRange && xaxis->GetBinUpEdge(x)<maxRange){
      if(MinVale > hist->GetBinContent(x))MinVale = hist->GetBinContent(x);
    }
  }
  return MinVale;
}

void MakePlot(TString fname1, TString fname2,TString fname3, TString fname4)
{
  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 800,630);
  gStyle->SetOptStat(0); //Set 0 to disable statistic box
  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  TKey *key;
  TFile *infile = TFile::Open(fname1.Data());

  TFile *infile2 = NULL;
  if(fname2 != DUMMYFILE)infile2 = TFile::Open(fname2.Data());
  TFile *infile3 = NULL;
  if(fname3 != DUMMYFILE) infile3 = TFile::Open(fname3.Data());
  TFile *infile4 = NULL;
  if(fname4 != DUMMYFILE) infile4 = TFile::Open(fname4.Data());

  TIter next(infile->GetListOfKeys());
  while ((key = static_cast<TKey*>(next()))) {
    std::string dirname = std::string(key->GetName());
    if (std::string(key->GetClassName()) != "TDirectoryFile") continue;
    //KS: Script will work with LogL and Batched_means, you can comment it if you are interested in it
    if( (dirname == "LogL") || (dirname == "Batched_means") ) continue;
    //KS: Trace wo longer chains is super big, the way to avoid is to plot as png but I don't like png,
    //keep possibility to skip it
    //if( (dirname == "Trace") ) continue;
    infile->cd(dirname.c_str());
    TIter nextsub(gDirectory->GetListOfKeys());
    c1->Print(Form("%s.pdf[",dirname.c_str()), "pdf");
    TKey *subkey;
    while ((subkey = static_cast<TKey*>(nextsub())))
    {
      std::string name = std::string(subkey->GetName());
      name = dirname + "/" + name;
      MACH3LOG_INFO("{}", name);
      if (std::string(subkey->GetClassName()) != "TH1D"){continue;}
      else{MACH3LOG_WARN("continuing along my way for {}", dirname);}

      TH1D* blarb[4];
      MACH3LOG_INFO("Looking for {} from file {}", name.c_str(), fname1.Data());
      blarb[0] = static_cast<TH1D*>(infile->Get(name.c_str())->Clone());
      //KS: Some fixe params can go crazy
      if(TMath::IsNaN(blarb[0]->GetBinContent(1)) ) continue;

      RemoveFitter(blarb[0], "Fitter");
      blarb[0]->SetLineStyle(kSolid);
      blarb[0]->SetLineColor(kRed);
      blarb[0]->Draw();
      if( dirname == "AccProb") blarb[0]->GetYaxis()->SetRangeUser(0, 1.0);
      if( name == "AccProb/AcceptanceProbability" ) continue;
      if(infile2 != NULL)
      {
        blarb[1] = static_cast<TH1D*>(infile2->Get(name.c_str())->Clone());
        RemoveFitter(blarb[1], "Fitter");
        blarb[1]->SetLineStyle(kDashed);
        blarb[1]->SetLineColor(kBlue);
        blarb[1]->Draw("same");
      }
      if(infile3 != NULL)
      {
        blarb[2] = static_cast<TH1D*>(infile3->Get(name.c_str())->Clone());
        RemoveFitter(blarb[2], "Fitter");
        blarb[2]->SetLineStyle(kDotted );
        blarb[2]->SetLineColor(kGreen);
        blarb[2]->Draw("same");
      }
      if(infile4 != NULL)
      {
        blarb[3] = static_cast<TH1D*>(infile4->Get(name.c_str())->Clone());
        RemoveFitter(blarb[3], "Fitter");
        blarb[3]->SetLineStyle(kDashDotted );
        blarb[3]->SetLineColor(kOrange);
        blarb[3]->Draw("same");
      }

      c1->Print(Form("%s.pdf",dirname.c_str()), "pdf");
    }
    gDirectory->cd("..");
    c1->Print(Form("%s.pdf]",dirname.c_str()), "pdf");
  }

  infile->Close();
  if(infile2 != NULL)infile2->Close();
  if(infile3 != NULL)infile3->Close();
  if(infile4 != NULL)infile4->Close();
}

void PlotAutoCorr(TString fname1, TString fname2, TString fname3, TString fname4)
{
  TString fname[4];
  fname[0] = fname1; fname[1] = fname2; fname[2] = fname3; fname[3] = fname4;
  //Color_t PlotColor[4]={kRed, kBlue, kGreen, kOrange};

  TFile *infile[4];
  infile[0] = TFile::Open(fname[0].Data());
  //KS" We need to check number of files to loop over in very lazy way
  int Nfiles = 1;

  if(fname[1] != DUMMYFILE){ infile[1] = TFile::Open(fname[1].Data()); Nfiles++;}
  if(fname[2] != DUMMYFILE){ infile[2] = TFile::Open(fname[2].Data()); Nfiles++;}
  if(fname[3] != DUMMYFILE){ infile[3] = TFile::Open(fname[3].Data()); Nfiles++;}

  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 800,630);
  gStyle->SetOptStat(0); //Set 0 to disable statistic box
  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  c1->Print("Auto_Corr_PerFile.pdf[", "pdf");
  for(int ik = 0; ik < Nfiles; ik++)
  {
    TIter next(infile[ik]->GetListOfKeys());

    TKey* key;
    while ((key = static_cast<TKey*>(next())))
    {
      std::string dirname = std::string(key->GetName());

      //KS: Script We are only interested in auto corr
      if( (dirname != "Auto_corr")) continue;

      infile[ik]->cd(dirname.c_str());
      TIter nextsub(gDirectory->GetListOfKeys());

      TKey *subkey;
      bool FirstTime = true;
      while ((subkey = static_cast<TKey*>(nextsub())))
      {
        std::string name = std::string(subkey->GetName());
        name = dirname + "/" + name;

        if (std::string(subkey->GetClassName()) != "TH1D") continue;
        MACH3LOG_DEBUG("{}", name.c_str());
        TH1D* blarb = static_cast<TH1D*>(infile[ik]->Get(name.c_str())->Clone());
        //KS: Some fixe pramas can go crazy
        if(TMath::IsNaN(blarb->GetBinContent(1))) continue;
        //KS: This is unfortunately hardcoded, need to find better way to write this
        //blarb[0]->GetListOfFunctions()->ls();
        delete blarb->GetListOfFunctions()->FindObject("Fitter");

        double MinValue = GetMinimumInRange(blarb, 0, 24000);

        if(MinValue >= 0.80) blarb->SetLineColor(kRed);
        else if(MinValue >= 0.40 && MinValue < 0.80) blarb->SetLineColor(kOrange);
        else if(MinValue > 0.20 && MinValue < 0.40) blarb->SetLineColor(kYellow);
        else if(MinValue <= 0.20) blarb->SetLineColor(kGreen);
        blarb->GetXaxis()->UnZoom();

        if(FirstTime) blarb->SetTitle( Form("Auto_Corr_%s.pdf",fname[ik].Data()) );

        blarb->SetLineStyle(kDashed);

        if(FirstTime) blarb->Draw();
        if(!FirstTime) blarb->Draw("same");
        FirstTime = false;
      }
      gDirectory->cd("..");
    }
    c1->Print("Auto_Corr_PerFile.pdf", "pdf");
  }
  c1->Print("Auto_Corr_PerFile.pdf]", "pdf");
}

std::vector<TH1D*> GetAverageFileAC(TFile* infile, Color_t color = kOrange+7){
  /*
  Generates the average auto-correlation histogram from a given file.
  */
  TDirectoryFile *ac_dir = static_cast<TDirectoryFile*>(infile->Get("Auto_corr"));
  if(!ac_dir) {
    MACH3LOG_ERROR("No Auto_corr directory in file {}", infile->GetName());
    return {nullptr, nullptr, nullptr};
  }

  TKey *key;
  TIter next(ac_dir->GetListOfKeys());

  TH1D* auto_hist = nullptr;
  TH1D* min_hist = nullptr;
  TH1D* max_hist = nullptr;

  int n_pars = 0;
  while ((key = static_cast<TKey*>(next()))) {
    std::string name = std::string(key->GetName());
    if (std::string(key->GetClassName()) != "TH1D") continue;

    auto plot_hist = static_cast<TH1D*>(ac_dir->Get(name.c_str())->Clone());

    // Lets us get binning info
    if(!auto_hist){
      auto_hist = static_cast<TH1D*>(plot_hist->Clone("Average_Auto_Corr"));
      auto_hist->SetDirectory(0); // Prevent ROOT from deleting it
      min_hist = static_cast<TH1D*>(plot_hist->Clone("Average_Auto_Corr_Min_Max"));
      max_hist = static_cast<TH1D*>(plot_hist->Clone("Average_Auto_Corr_Min_Max"));
      auto_hist->Reset();
    }
    // Now we add the hists together
    auto_hist->Add(plot_hist);
    
    // Now now iterate to get min and max hists
    for(int i = 1; i <= plot_hist->GetNbinsX(); i++) {
      double bin_content = plot_hist->GetBinContent(i);
      if(bin_content < min_hist->GetBinContent(i)) {
        min_hist->SetBinContent(i, bin_content);
      }
      if(bin_content > max_hist->GetBinContent(i)) {
        max_hist->SetBinContent(i, bin_content);
      }
    }

    n_pars++;
  
  }

  auto_hist->Scale(1.0 / float(n_pars)); // Average the histogram

  // We now want to make things look pretty!
  auto error_hist = static_cast<TH1D*>(auto_hist->Clone("Average_Auto_Corr_Error"));
  error_hist->SetDirectory(0); // Prevent ROOT from deleting it

  auto_hist->SetLineColor(color);
  error_hist->SetFillColorAlpha(color, float(0.3));
  min_hist->SetFillColorAlpha(color, float(0.1));
  max_hist->SetFillColorAlpha(color, float(0.1));


  return {auto_hist, error_hist, min_hist, max_hist};
}

void DrawAverageAC(TFile* file, Color_t color = kOrange+7, TCanvas* canvas=nullptr, TLegend* legend = nullptr){
  auto ac_hists = GetAverageFileAC(file, color);

  if(ac_hists[0]==nullptr){
    MACH3LOG_ERROR("No auto-correlation histograms found in file {}", file->GetName());
    return;
  }
  
  canvas->cd();

  ac_hists[1]->SetTitle("Average Auto-Correlation");
  ac_hists[1]->GetXaxis()->SetTitle("Lag");
  ac_hists[1]->GetYaxis()->SetTitle("Auto-Correlation");


  // First we need to proces the minmnax bbands
  auto band = new TGraphAsymmErrors(ac_hists[2]->GetNbinsX());
  for(int i=1; i<=ac_hists[2]->GetNbinsX(); i++){
    double x = ac_hists[2]->GetBinCenter(i);
    double y_min = ac_hists[2]->GetBinContent(i);
    double y_max = ac_hists[3]->GetBinContent(i);

    band->SetPoint(i-1, x, (y_min + y_max) / 2.0);
    band->SetPointError(i-1, 0.0, 0.0, (y_max-y_min) / 2.0, (y_max-y_min) / 2.0);
  }

  if(legend){
    legend->AddEntry(ac_hists[0], TString(file->GetName()), "l");
  }

  // Now we draw!
  ac_hists[1]->Draw("SAME E3");
  band->Draw("SAME");
  ac_hists[0]->Draw("SAME HIST");

}


// TODO [HW]: Make this all a little more generic!
void PlotAverageAutoCorr(TString fname1, TString fname2, TString fname3, TString fname4)
{
  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 800,630);
  TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  gStyle->SetOptStat(0); //Set 0 to disable statistic box
  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  TFile *infile = TFile::Open(fname1.Data());
  DrawAverageAC(infile, kRed, c1, leg);

  TFile *infile2 = NULL;
  if(fname2 != DUMMYFILE){
    infile2 = TFile::Open(fname2.Data());
    DrawAverageAC(infile2, kBlue, c1, leg);
  };
  TFile *infile3 = NULL;
  if(fname3 != DUMMYFILE){
     infile3 = TFile::Open(fname3.Data());
     DrawAverageAC(infile3, kGreen, c1, leg);
  }
  TFile *infile4 = NULL;
  if(fname4 != DUMMYFILE){
    infile4 = TFile::Open(fname4.Data());
    DrawAverageAC(infile4, kOrange, c1, leg);
  }
    
  c1->Print("Average_Auto_Corr.pdf", "pdf");
  
}

int main(int argc, char *argv[]) {
  SetMaCh3LoggerFormat();
  if (argc < 2 || argc > 5)
  {
    MACH3LOG_ERROR("How to use: {} DiagMCMC_Output.root", argv[0]);
    MACH3LOG_ERROR("Up to 4 files");
    throw MaCh3Exception(__FILE__ , __LINE__);
  }

  if(argc == 2) {
    MakePlot(argv[1], DUMMYFILE, DUMMYFILE, DUMMYFILE);
    PlotAutoCorr(argv[1], DUMMYFILE, DUMMYFILE, DUMMYFILE);
    PlotAverageAutoCorr(argv[1], DUMMYFILE, DUMMYFILE, DUMMYFILE);
  } else if(argc == 3) {
    MakePlot(argv[1], argv[2], DUMMYFILE ,DUMMYFILE);
    PlotAutoCorr(argv[1], argv[2], DUMMYFILE, DUMMYFILE);
    PlotAverageAutoCorr(argv[1], argv[2], DUMMYFILE, DUMMYFILE);
  } else if(argc == 4) {
    MakePlot(argv[1], argv[2], argv[3], DUMMYFILE);
    PlotAutoCorr(argv[1], argv[2], argv[3], DUMMYFILE);
    PlotAverageAutoCorr(argv[1], argv[2], argv[3], DUMMYFILE);
  } else if(argc == 5) {
    MakePlot(argv[1], argv[2], argv[3], argv[4]);
    PlotAutoCorr(argv[1], argv[2], argv[3], argv[4]);
    PlotAverageAutoCorr(argv[1], argv[2], argv[3], argv[4]);
  }
  return 0;
}
