// MaCh3 includes
#include "Fitters/MCMCProcessor.h"
#include "Samples/HistogramUtils.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include <TSystem.h>

/// @file PlotMCMCDiag.cpp
/// @brief KS: This script is used to analyse output form DiagMCMC.
/// @warning This script support comparing up to 4 files, there is easy way to expand it up to five or six,
/// @todo this need serious refactor

TString DUMMYFILE = "DummyFile";
TString DUMMYNAME = "DummyName";

// Utilities

/// @brief KS: function which looks for minimum in given range
double GetMinimumInRange(TH1D *hist, double minRange, double maxRange)
{
  double MinVale = 1234567890.;
  TAxis *xaxis = hist->GetXaxis();
  for (int x = 1; x <= hist->GetNbinsX(); x++)
  {
    if (xaxis->GetBinLowEdge(x) > minRange && xaxis->GetBinUpEdge(x) < maxRange)
    {
      if (MinVale > hist->GetBinContent(x))
        MinVale = hist->GetBinContent(x);
    }
  }
  return MinVale;
}

/// @brief HW: Check if histogram is flat within a given tolerance.
// Utility functions
bool IsHistogramAllOnes(TH1D *hist, double tolerance = 0.001, int max_failures = 100)
{
  int failure_count = 0;

  for (int bin = 2; bin <= hist->GetNbinsX(); ++bin)
  {
    if (fabs(hist->GetBinContent(bin) - 1.0) > tolerance)
    {
      if (++failure_count > max_failures)
      {
        return false;
      }
    }
  }
  return true;
}

void MakePlot(TString fname1, TString flabel1, TString fname2, TString flabel2, TString fname3, TString flabel3, TString fname4, TString flabel4)
{
  TCanvas *c1 = new TCanvas("c1", " ", 0, 0, 800, 630);
  gStyle->SetOptStat(0); // Set 0 to disable statistic box
  // To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  TKey *key;
  TFile *infile = TFile::Open(fname1.Data());

  bool add_legend = false;

  TFile *infile2 = NULL;
  if (fname2 != DUMMYFILE)
  {
    infile2 = TFile::Open(fname2.Data());
    add_legend = true;
  }
  TFile *infile3 = NULL;
  if (fname3 != DUMMYFILE)
  {
    infile3 = TFile::Open(fname3.Data());
    add_legend = true;
  }
  TFile *infile4 = NULL;
  if (fname4 != DUMMYFILE)
  {
    infile4 = TFile::Open(fname4.Data());
    add_legend = true;
  }

  TIter next(infile->GetListOfKeys());
  while ((key = static_cast<TKey *>(next())))
  {
    std::string dirname = std::string(key->GetName());
    if (std::string(key->GetClassName()) != "TDirectoryFile")
      continue;
    // KS: Script will work with LogL and Batched_means, you can comment it if you are interested in it
    if ((dirname == "LogL") || (dirname == "Batched_means"))
      continue;
    // KS: Trace wo longer chains is super big, the way to avoid is to plot as png but I don't like png,
    // keep possibility to skip it
    // if( (dirname == "Trace") ) continue;
    infile->cd(dirname.c_str());
    TIter nextsub(gDirectory->GetListOfKeys());
    c1->Print(Form("%s_%s.pdf[", flabel1.Data(), dirname.c_str()), "pdf");
    TKey *subkey;
    while ((subkey = static_cast<TKey *>(nextsub())))
    {

      TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->SetFillColorAlpha(kWhite, float(0.7));
      std::string name = std::string(subkey->GetName());
      name = dirname + "/" + name;
      MACH3LOG_INFO("{}", name);
      if (std::string(subkey->GetClassName()) != "TH1D")
      {
        continue;
      }
      else
      {
        MACH3LOG_WARN("continuing along my way for {}", dirname);
      }

      TH1D *blarb[4];
      MACH3LOG_INFO("Looking for {} from file {}", name.c_str(), fname1.Data());
      blarb[0] = static_cast<TH1D *>(infile->Get(name.c_str())->Clone());
      // KS: Some fixe params can go crazy
      if (TMath::IsNaN(blarb[0]->GetBinContent(1)))
        continue;

      RemoveFitter(blarb[0], "Fitter");
      blarb[0]->SetLineStyle(kSolid);
      blarb[0]->SetLineColor(kRed);
      blarb[0]->Draw();

      leg->AddEntry(blarb[0], flabel1.Data(), "l");

      if (dirname == "AccProb")
        blarb[0]->GetYaxis()->SetRangeUser(0, 1.0);
      if (name == "AccProb/AcceptanceProbability")
        continue;
      if (infile2 != NULL)
      {
        blarb[1] = static_cast<TH1D *>(infile2->Get(name.c_str())->Clone());
        RemoveFitter(blarb[1], "Fitter");
        blarb[1]->SetLineStyle(kDashed);
        blarb[1]->SetLineColor(kBlue);
        blarb[1]->Draw("same");
        leg->AddEntry(blarb[1], flabel2.Data(), "l");
      }
      if (infile3 != NULL)
      {
        blarb[2] = static_cast<TH1D *>(infile3->Get(name.c_str())->Clone());
        RemoveFitter(blarb[2], "Fitter");
        blarb[2]->SetLineStyle(kDotted);
        blarb[2]->SetLineColor(kGreen);
        blarb[2]->Draw("same");
        leg->AddEntry(blarb[2], flabel3.Data(), "l");
      }
      if (infile4 != NULL)
      {
        blarb[3] = static_cast<TH1D *>(infile4->Get(name.c_str())->Clone());
        RemoveFitter(blarb[3], "Fitter");
        blarb[3]->SetLineStyle(kDashDotted);
        blarb[3]->SetLineColor(kOrange);
        blarb[3]->Draw("same");
        leg->AddEntry(blarb[3], flabel4.Data(), "l");
      }
      if (add_legend)
      {
        leg->Draw();
      }

      c1->Print(Form("%s_%s.pdf", flabel1.Data(), dirname.c_str()), "pdf");
      delete leg;
    }
    gDirectory->cd("..");
    c1->Print(Form("%s_%s.pdf]", flabel1.Data(), dirname.c_str()), "pdf");
  }

  infile->Close();
  if (infile2 != NULL)
    infile2->Close();
  if (infile3 != NULL)
    infile3->Close();
  if (infile4 != NULL)
    infile4->Close();
}

void PlotAutoCorr(TString fname1, TString flabel1, TString fname2, TString flabel2, TString fname3, TString flabel3, TString fname4, TString flabel4)
{
  TString fname[4];
  fname[0] = fname1;
  fname[1] = fname2;
  fname[2] = fname3;
  fname[3] = fname4;
  // Color_t PlotColor[4]={kRed, kBlue, kGreen, kOrange};
  std::vector<TString> flabel = {flabel1, flabel2, flabel3, flabel4};

  TFile *infile[4];
  infile[0] = TFile::Open(fname[0].Data());
  // KS" We need to check number of files to loop over in very lazy way
  int Nfiles = 1;

  if (fname[1] != DUMMYFILE)
  {
    infile[1] = TFile::Open(fname[1].Data());
    Nfiles++;
  }
  if (fname[2] != DUMMYFILE)
  {
    infile[2] = TFile::Open(fname[2].Data());
    Nfiles++;
  }
  if (fname[3] != DUMMYFILE)
  {
    infile[3] = TFile::Open(fname[3].Data());
    Nfiles++;
  }

  TCanvas *c1 = new TCanvas("c1", " ", 0, 0, 800, 630);
  gStyle->SetOptStat(0); // Set 0 to disable statistic box
  // To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  c1->Print("Auto_Corr_PerFile.pdf[", "pdf");
  for (int ik = 0; ik < Nfiles; ik++)
  {
    TIter next(infile[ik]->GetListOfKeys());

    TKey *key;
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

    while ((key = static_cast<TKey *>(next())))
    {
      std::string dirname = std::string(key->GetName());

      // KS: Script We are only interested in auto corr
      if ((dirname != "Auto_corr"))
        continue;

      infile[ik]->cd(dirname.c_str());
      TIter nextsub(gDirectory->GetListOfKeys());

      TKey *subkey;
      bool FirstTime = true;

      while ((subkey = static_cast<TKey *>(nextsub())))
      {
        std::string name = std::string(subkey->GetName());
        name = dirname + "/" + name;

        if (std::string(subkey->GetClassName()) != "TH1D")
          continue;
        MACH3LOG_DEBUG("{}", name.c_str());
        TH1D *blarb = static_cast<TH1D *>(infile[ik]->Get(name.c_str())->Clone());
        // KS: Some fixe pramas can go crazy
        if (TMath::IsNaN(blarb->GetBinContent(1)))
          continue;
        // KS: This is unfortunately hardcoded, need to find better way to write this
        // blarb[0]->GetListOfFunctions()->ls();
        delete blarb->GetListOfFunctions()->FindObject("Fitter");

        double MinValue = GetMinimumInRange(blarb, 0, 24000);

        if (MinValue >= 0.80)
          blarb->SetLineColor(kRed);
        else if (MinValue >= 0.40 && MinValue < 0.80)
          blarb->SetLineColor(kOrange);
        else if (MinValue > 0.20 && MinValue < 0.40)
          blarb->SetLineColor(kYellow);
        else if (MinValue <= 0.20)
          blarb->SetLineColor(kGreen);
        blarb->GetXaxis()->UnZoom();

        if (FirstTime)
          blarb->SetTitle(Form("Auto_Corr_%s.pdf", fname[ik].Data()));

        blarb->SetLineStyle(kDashed);

        if (FirstTime)
          blarb->Draw();

        TString integral = Form("Integral: %.2f", blarb->Integral());

        if (!FirstTime)
          blarb->Draw("same");
        FirstTime = false;
      }
      gDirectory->cd("..");
    }

    if (Nfiles > 1)
    {
      leg->Draw();
    }
    c1->Print("Auto_Corr_PerFile.pdf", "pdf");
    delete leg;
  }
  c1->Print("Auto_Corr_PerFile.pdf]", "pdf");
}

// HW: Utilities for average Auto Correlation plots

/// @brief HW: Create a band of minimum and maximum values from a histogram.
std::pair<TGraph *, TGraph *> CreateMinMaxBand(TH1D *hist, Color_t color)
{
  int nBins = hist->GetNbinsX();
  std::vector<double> x(nBins), ymin(nBins), ymax(nBins);

  for (int i = 0; i < nBins; ++i)
  {
    x[i] = hist->GetBinCenter(i + 1);
    ymin[i] = hist->GetBinContent(i + 1);
    ymax[i] = hist->GetBinContent(i + 1);
  }

  TGraph *minGraph = new TGraph(nBins, x.data(), ymin.data());
  TGraph *maxGraph = new TGraph(nBins, x.data(), ymax.data());

  minGraph->SetLineColor(color);
  maxGraph->SetLineColor(color);
  minGraph->SetFillColorAlpha(color, float(0.3));
  maxGraph->SetFillColorAlpha(color, float(0.3));

  return {minGraph, maxGraph};
}

TGraph *CalculateMinMaxBand(const std::vector<TH1D *> &histograms, Color_t color)
{
  if (histograms.empty())
  {
    throw std::invalid_argument("Empty histogram vector provided");
  }

  int nBins = histograms[0]->GetNbinsX();
  std::vector<double> x(nBins), ymin(nBins, 1e10), ymax(nBins, -1e10);

  for (int bin = 0; bin < nBins; bin++)
  {
    x[bin] = histograms[0]->GetBinCenter(bin + 1);

    for (const auto &hist : histograms)
    {
      double content = hist->GetBinContent(bin + 1);
      ymin[bin] = std::min(ymin[bin], content);
      ymax[bin] = std::max(ymax[bin], content);
    }
  }

  // Create a band using TGraphAsymmErrors for the shaded region
  TGraphAsymmErrors *band = new TGraphAsymmErrors(nBins);
  for (int i = 0; i < nBins; ++i)
  {
    band->SetPoint(i, x[i], (ymax[i] + ymin[i]) / 2.0);
    band->SetPointError(i, 0, 0, (ymax[i] - ymin[i]) / 2.0, (ymax[i] - ymin[i]) / 2.0);
  }

  band->SetFillColorAlpha(color, float(0.2));
  band->SetLineWidth(1);

  return band; // Return band and min graph (min graph can be used for legend)
}

std::pair<TH1D *, TH1D *> CalculateMinMaxHistograms(const std::vector<TH1D *> &histograms)
{
  if (histograms.empty())
  {
    throw std::invalid_argument("Empty histogram vector provided");
  }

  TH1D *min_hist = static_cast<TH1D *>(histograms[0]->Clone());
  TH1D *max_hist = static_cast<TH1D *>(histograms[0]->Clone());

  for (const auto &hist : histograms)
  {
    for (int bin = 1; bin <= min_hist->GetNbinsX(); ++bin)
    {
      double current_min = min_hist->GetBinContent(bin);
      double current_max = max_hist->GetBinContent(bin);
      double bin_content = hist->GetBinContent(bin);

      min_hist->SetBinContent(bin, std::min(current_min, bin_content));
      max_hist->SetBinContent(bin, std::max(current_max, bin_content));
    }
  }

  return {min_hist, max_hist};
}

// File processing functions
void ProcessAutoCorrelationDirectory(TDirectoryFile *autocor_dir,
                                     TH1D *&average_hist,
                                     int &parameter_count,
                                     std::vector<TH1D *> &histograms)
{
  TIter next(autocor_dir->GetListOfKeys());
  TKey *key;

  while ((key = dynamic_cast<TKey *>(next())))
  {
    TH1D *current_hist = nullptr;
    autocor_dir->GetObject(key->GetName(), current_hist);

    if (!current_hist ||
        current_hist->GetMaximum() <= 0 ||
        IsHistogramAllOnes(current_hist))
    {
      continue;
    }

    current_hist->SetDirectory(nullptr); // Detach from file
    histograms.push_back(current_hist);

    if (!average_hist)
    {
      average_hist = static_cast<TH1D *>(current_hist->Clone());
      average_hist->SetDirectory(nullptr);
    }
    else
    {
      average_hist->Add(current_hist);
    }
    parameter_count++;
  }
}

void ProcessDiagnosticFile(const TString &file_path,
                           TH1D *&average_hist,
                           int &parameter_count,
                           std::vector<TH1D *> &histograms)
{
  std::unique_ptr<TFile> input_file(TFile::Open(file_path));
  if (!input_file || input_file->IsZombie())
  {
    throw std::runtime_error("Could not open file: " + std::string(file_path.Data()));
  }

  TDirectoryFile *autocor_dir = nullptr;
  input_file->GetObject("Auto_corr", autocor_dir);

  if (!autocor_dir)
  {
    throw MaCh3Exception(__FILE__, __LINE__,
                         "Auto_corr directory not found in file: " + std::string(file_path.Data()));
  }

  ProcessAutoCorrelationDirectory(autocor_dir, average_hist, parameter_count, histograms);
}

TH1D *AutocorrProcessInputs(const TString &input_file, std::vector<TH1D *> &histograms)
{

  TH1D *average_hist = nullptr;
  int parameter_count = 0;

  try
  {
    ProcessDiagnosticFile(input_file, average_hist, parameter_count, histograms);
  }
  catch (const std::exception &e)
  {
    MACH3LOG_ERROR("Error processing file {} : {}", input_file, e.what());
  }

  if (average_hist && parameter_count > 0)
  {
    MACH3LOG_INFO("Processed {} parameters from {} files", parameter_count, input_file);
    average_hist->Scale(1.0 / parameter_count);
  }

  return average_hist;
}

void CompareAverageAC(const std::vector<std::vector<TH1D *>> &histograms,
                      const std::vector<TH1D *> &averages,
                      const std::vector<TString> &hist_labels,
                      const TString &output_name,
                      bool draw_min_max = true,
                      bool draw_all = false,
                      bool draw_errors = true)
{
  TCanvas *canvas = new TCanvas("AverageAC", "Average Auto Correlation", 800, 600);
  canvas->SetGrid();

  // Setup colour palette
  Int_t nb = 255.0;
  Double_t stops[8] = {0.0, 0.25, 0.5, 0.75};
  Double_t red[8] = {0.83203125, 0.796875, 0.0, 0.9375};
  Double_t green[8] = {0.3671875, 0.47265625, 0.4453125, 0.890625};
  Double_t blue[8] = {0.0, 0.65234375, 0.6953125, 0.2578125};
  TColor::CreateGradientColorTable(8, stops, red, green, blue, nb);

  TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
  leg->SetFillColorAlpha(kWhite, float(0.7));

  if (draw_min_max)
  {
    for (size_t i = 0; i < histograms.size(); ++i)
    {
      auto colour = static_cast<Color_t>(TColor::GetColorPalette(static_cast<Int_t>(i * nb / histograms.size())));

      {
        auto band = CalculateMinMaxBand(histograms[i], colour);
        band->SetLineWidth(0);
        band->SetTitle("Average Auto Correlation");
        band->GetXaxis()->SetTitle("lag");
        band->GetYaxis()->SetTitle("Autocorrelation Function");
        if (i == 0)
        {
          band->Draw("A3");
        }
        else
        {
          band->Draw("3 SAME");
        }
      }
    }
  }

  for (size_t i = 0; i < averages.size(); ++i)
  {
    auto colour = static_cast<Color_t>(TColor::GetColorPalette(static_cast<Int_t>(i * nb / histograms.size())));

    if (draw_errors)
    {
      auto error_hist = static_cast<TH1D *>(averages[i]->Clone());
      error_hist->SetFillColorAlpha(colour, float(0.3));
      error_hist->SetLineWidth(0);
      error_hist->SetTitle("Average Auto Correlation");
      error_hist->GetXaxis()->SetTitle("lag");
      error_hist->GetYaxis()->SetTitle("Autocorrelation Function");
      error_hist->Draw("E3 SAME");
    }
    if (draw_all)
    {
      for (const auto &hist : histograms[i])
      {
        hist->SetLineColorAlpha(colour, float(0.05));
        hist->SetLineWidth(1);
        hist->Draw("HIST SAME");
      }
    }

    averages[i]->SetLineColor(colour);
    averages[i]->SetLineWidth(2);
    averages[i]->SetTitle("Average Auto Correlation");
    averages[i]->GetXaxis()->SetTitle("lag");
    averages[i]->GetYaxis()->SetTitle("Autocorrelation Function");
    averages[i]->Draw("HIST SAME");

    TString int_str = TString::Format("Average integrated autocorrelation %f", averages[i]->Integral());

    leg->AddEntry(averages[i], hist_labels[i] + int_str, "l");
  }

  if (averages.size() > 1)
  {
    leg->Draw();
  }

  canvas->SaveAs(output_name + ".pdf");
  delete canvas;
}

void PlotAverageACMult(std::vector<TString> input_files,
                       std::vector<TString> hist_labels,
                       const TString &output_name,
                       bool draw_min_max = true)
// bool draw_all = false,
// bool draw_errors = true)
{
  // Process first folder
  std::vector<std::vector<TH1D *>> histograms;
  std::vector<TH1D *> averages;

  for (int i = 0; i < static_cast<int>(input_files.size()); i++)
  {
    TString folder = input_files[i];
    MACH3LOG_INFO("Folder : {}", folder);

    if (folder.IsNull() || folder == DUMMYFILE)
    {
      MACH3LOG_WARN("Skipping empty or dummy folder: {}", folder.Data());
      continue;
    }

    std::vector<TH1D *> histograms_i;
    averages.push_back((AutocorrProcessInputs(folder, histograms_i)));
    histograms.push_back(histograms_i);
  }

  CompareAverageAC(histograms, averages, hist_labels, output_name, draw_min_max); // draw_all, draw_errors);
}

int main(int argc, char *argv[])
{
  SetMaCh3LoggerFormat();
  if (argc != 3 && argc != 5 && argc != 7 && argc != 9)
  {
    MACH3LOG_ERROR("Wrong number of arguments ({}) provided", argc);
    MACH3LOG_ERROR("How to use: {} DiagMCMC_Output.root Plot Name", argv[0]);
    MACH3LOG_ERROR("Up to 4 files");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (argc == 3)
  {
    MakePlot(argv[1], argv[2], DUMMYFILE, DUMMYNAME, DUMMYFILE, DUMMYNAME, DUMMYFILE, DUMMYNAME);
    PlotAutoCorr(argv[1], argv[2], DUMMYFILE, DUMMYNAME, DUMMYFILE, DUMMYNAME, DUMMYFILE, DUMMYNAME);
    PlotAverageACMult({argv[1]}, {argv[2]}, Form("%s_Average_Auto_Corr", argv[2]), true);
  }
  else if (argc == 5)
  {
    MakePlot(argv[1], argv[2], argv[3], argv[4], DUMMYFILE, DUMMYNAME, DUMMYFILE, DUMMYNAME);
    PlotAutoCorr(argv[1], argv[2], argv[3], argv[4], DUMMYFILE, DUMMYNAME, DUMMYFILE, DUMMYNAME);
    PlotAverageACMult({argv[1], argv[3]}, {argv[2], argv[4]}, Form("%s_Average_Auto_Corr", argv[2]), true);
  }
  else if (argc == 7)
  {
    MakePlot(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], DUMMYFILE, DUMMYNAME);
    PlotAutoCorr(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], DUMMYFILE, DUMMYNAME);
    PlotAverageACMult({argv[1], argv[3], argv[5]}, {argv[2], argv[4], argv[6]}, Form("%s_Average_Auto_Corr", argv[2]), true);
  }
  else if (argc == 9)
  {
    MakePlot(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
    PlotAutoCorr(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
    PlotAverageACMult({argv[1], argv[3], argv[5], argv[7]}, {argv[2], argv[4], argv[6], argv[8]}, Form("%s_Average_Auto_Corr", argv[2]), true);
  }
  return 0;
}
