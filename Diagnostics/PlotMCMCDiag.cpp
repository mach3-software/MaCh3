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
/// @author Henry Wallace
/// @author Kamil Skwarczynski

// Utilities
/// @brief KS: function which looks for minimum in given range
double GetMinimumInRange(TH1D *hist, const double minRange, const double maxRange)
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
bool IsHistogramAllOnes(TH1D *hist, double tolerance = 0.001, int max_failures = 100)
{
  int failure_count = 0;

  for (int bin = 2; bin <= hist->GetNbinsX(); ++bin)
  {
    if (std::fabs(hist->GetBinContent(bin) - 1.0) > tolerance)
    {
      if (++failure_count > max_failures)
      {
        return false;
      }
    }
  }
  return true;
}

/// @brief function which loops over Diag MCMC output and prints everything to PDF
void MakeDiagPlot(std::vector<TString> input_files,
              std::vector<TString> hist_labels)
{
  auto c1 = std::make_unique<TCanvas>("c1", " ", 0, 0, 800, 630);
  gStyle->SetOptStat(0); // Set 0 to disable statistic box
  // To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  TKey *key = nullptr;
  std::vector<TFile*> infiles;
  bool add_legend = false;

  for (const auto& fname : input_files)
  {
    TFile* f = M3::Open(fname.Data(), "UPDATE", __FILE__, __LINE__);
    MACH3LOG_INFO("Adding file {}", fname.Data());
    infiles.push_back(f);
    add_legend = (input_files.size() > 1);
  }
  
  constexpr int NFiles_hard = 4;
  constexpr Color_t Colour[NFiles_hard] = {kRed, kBlue, kGreen, kOrange};
  constexpr Style_t Style[NFiles_hard] = {kSolid, kDashed, kDotted, kDashDotted};
  if(input_files.size() > NFiles_hard){
    MACH3LOG_ERROR("Passed {}, while I have nasty hardcoding for max {} filse", input_files.size(), NFiles_hard);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  TIter next(infiles[0]->GetListOfKeys());
  while ((key = static_cast<TKey *>(next())))
  {
    std::string dirname = std::string(key->GetName());
    if (std::string(key->GetClassName()) != "TDirectoryFile")
      continue;
    // KS: Script will work with LogL and Batched_means, you can comment it if you are interested in it
    if ((dirname == "LogL") || (dirname == "Batched_means") || (dirname == "PowerSpectrum"))
      continue;
    // KS: Trace wo longer chains is super big, the way to avoid is to plot as png but I don't like png,
    // keep possibility to skip it
    // if( (dirname == "Trace") ) continue;
    infiles[0]->cd(dirname.c_str());
    TIter nextsub(gDirectory->GetListOfKeys());
    c1->Print(Form("%s_%s.pdf[", hist_labels[0].Data(), dirname.c_str()), "pdf");
    TKey *subkey;
    while ((subkey = static_cast<TKey *>(nextsub())))
    {
      auto leg = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
      leg->SetFillColorAlpha(kWhite, 0.7f);
      std::string name = std::string(subkey->GetName());
      name = dirname + "/" + name;
      MACH3LOG_INFO("{}", name);
      if (std::string(subkey->GetClassName()) != "TH1D") {
        continue;
      } else {
        MACH3LOG_WARN("continuing along my way for {}", dirname);
      }
      if (name == "AccProb/AcceptanceProbability") continue;

      std::vector<std::unique_ptr<TH1D>> blarb(infiles.size());
      MACH3LOG_INFO("Looking for {} from file {}", name.c_str(), input_files[0].Data());
      blarb[0] = M3::Clone(static_cast<TH1D *>(infiles[0]->Get(name.c_str())));
      // KS: Some fixe params can go crazy
      if (std::isnan(blarb[0]->GetBinContent(1)))
        continue;

      RemoveFitter(blarb[0].get(), "Fitter");
      if (dirname == "AccProb")
        blarb[0]->GetYaxis()->SetRangeUser(0, 1.0);

      blarb[0]->SetLineStyle(Style[0]);
      blarb[0]->SetLineColor(Colour[0]);
      blarb[0]->Draw();
      leg->AddEntry(blarb[0].get(), hist_labels[0].Data(), "l");

      for (size_t i = 1; i < infiles.size(); ++i)
      {
        blarb[i] = M3::Clone(static_cast<TH1D *>(infiles[i]->Get(name.c_str())));

        blarb[i]->SetLineStyle(Style[i]);
        blarb[i]->SetLineColor(Colour[i]);

        blarb[i]->Draw(i == 0 ? "" : "same");
        leg->AddEntry(blarb[i].get(), hist_labels[i].Data(), "l");
      }
      if (add_legend)
      {
        leg->Draw();
      }
      c1->Print(Form("%s_%s.pdf", hist_labels[0].Data(), dirname.c_str()), "pdf");
    }
    gDirectory->cd("..");
    c1->Print(Form("%s_%s.pdf]", hist_labels[0].Data(), dirname.c_str()), "pdf");
  }

  for (size_t i = 0; i < infiles.size(); ++i)
  {
    infiles[i]->Close();
    delete infiles[i];
  }
}

/// @brief Plot autocorrelation for all params into single plot with color codding helping whether on average it is ok or not
void PlotAutoCorr(const std::vector<TString>& fname)
{
  std::vector<TFile*> infile(fname.size());
  int Nfiles = 0;
  for(size_t i = 0; i < fname.size(); i++)
  {
    infile[i] = M3::Open(fname[i].Data(), "open", __FILE__, __LINE__);
    Nfiles++;
  }

  auto c1 = std::make_unique<TCanvas>("c1", " ", 0, 0, 800, 630);
  gStyle->SetOptStat(0); // Set 0 to disable statistic box
  // To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  c1->Print("Auto_Corr_PerFile.pdf[", "pdf");
  for (int ik = 0; ik < Nfiles; ik++)
  {
    TIter next(infile[ik]->GetListOfKeys());

    TKey *key;
    // KS: Keep track of histos to delete them
    std::vector<std::unique_ptr<TH1D>> histos;
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
        histos.push_back(M3::Clone(static_cast<TH1D *>(infile[ik]->Get(name.c_str()))));
        // KS: Some fixe pramas can go crazy
        if (std::isnan(histos.back()->GetBinContent(1)))
          continue;
        // KS: This is unfortunately hardcoded, need to find better way to write this
        // histos.back()->GetListOfFunctions()->ls();
        RemoveFitter(histos.back().get(), "Fitter");
        // KS: Very arbitrary check if autocorrelation dropped below threshold in half of consider LagL
        double UpperEdgeBeforeLast = histos.back()->GetXaxis()->GetBinUpEdge(histos.back()->GetNbinsX()/2);
        double MinValue = GetMinimumInRange(histos.back().get(), 0, UpperEdgeBeforeLast);

        if (MinValue >= 0.80)
          histos.back()->SetLineColor(kRed);
        else if (MinValue >= 0.40 && MinValue < 0.80)
          histos.back()->SetLineColor(kOrange);
        else if (MinValue > 0.20 && MinValue < 0.40)
          histos.back()->SetLineColor(kYellow);
        else if (MinValue <= 0.20)
          histos.back()->SetLineColor(kGreen);
        histos.back()->GetXaxis()->UnZoom();

        if (FirstTime)
          histos.back()->SetTitle(Form("Auto_Corr_%s.pdf", fname[ik].Data()));

        histos.back()->SetLineStyle(kDashed);

        if (FirstTime)
          histos.back()->Draw();

        TString integral = Form("Integral: %.2f", histos.back()->Integral());

        if (!FirstTime)
          histos.back()->Draw("same");
        FirstTime = false;
      }
      gDirectory->cd("..");
    }
    c1->Print("Auto_Corr_PerFile.pdf", "pdf");
  }
  c1->Print("Auto_Corr_PerFile.pdf]", "pdf");
}

// HW: Utilities for average Auto Correlation plots

/// @brief HW: Create a band of minimum and maximum values from a histogram.
std::pair<std::unique_ptr<TGraph>, std::unique_ptr<TGraph>> CreateMinMaxBand(TH1D *hist, Color_t color)
{
  int nBins = hist->GetNbinsX();
  std::vector<double> x(nBins), ymin(nBins), ymax(nBins);

  for (int i = 0; i < nBins; ++i)
  {
    x[i] = hist->GetBinCenter(i + 1);
    ymin[i] = hist->GetBinContent(i + 1);
    ymax[i] = hist->GetBinContent(i + 1);
  }

  auto minGraph = std::make_unique<TGraph>(nBins, x.data(), ymin.data());
  auto maxGraph = std::make_unique<TGraph>(nBins, x.data(), ymax.data());

  minGraph->SetLineColor(color);
  maxGraph->SetLineColor(color);
  minGraph->SetFillColorAlpha(color, 0.3f);
  maxGraph->SetFillColorAlpha(color, 0.3f);

  return {std::move(minGraph), std::move(maxGraph)};
}

std::unique_ptr<TGraphAsymmErrors> CalculateMinMaxBand(const std::vector<std::unique_ptr<TH1D>> &histograms, Color_t color)
{
  if (histograms.empty()) {
    throw MaCh3Exception(__FILE__, __LINE__, "Empty histogram vector provided");
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
  auto band = std::make_unique<TGraphAsymmErrors>(nBins);
  for (int i = 0; i < nBins; ++i)
  {
    band->SetPoint(i, x[i], (ymax[i] + ymin[i]) / 2.0);
    band->SetPointError(i, 0, 0, (ymax[i] - ymin[i]) / 2.0, (ymax[i] - ymin[i]) / 2.0);
  }

  band->SetFillColorAlpha(color, 0.2f);
  band->SetLineWidth(1);

  return band; // Return band and min graph (min graph can be used for legend)
}

/// @brief Loop over directory get histograms into vector and add to averaged hist
void ProcessAutoCorrelationDirectory(TDirectoryFile *autocor_dir,
                                     std::unique_ptr<TH1D>& average_hist,
                                     int &parameter_count,
                                     std::vector<std::unique_ptr<TH1D>> &histograms)
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

    histograms.push_back(M3::Clone(current_hist));
    if (!average_hist) {
      average_hist = M3::Clone(current_hist);
    } else {
      average_hist->Add(current_hist);
    }
    parameter_count++;
  }
}

void ProcessDiagnosticFile(const TString &file_path,
                           std::unique_ptr<TH1D>& average_hist,
                           int &parameter_count,
                           std::vector<std::unique_ptr<TH1D>> &histograms)
{
  std::unique_ptr<TFile> input_file(TFile::Open(file_path));
  if (!input_file || input_file->IsZombie())
  {
    throw MaCh3Exception(__FILE__, __LINE__, "Could not open file: " + std::string(file_path.Data()));
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

std::unique_ptr<TH1D> AutocorrProcessInputs(const TString &input_file, std::vector<std::unique_ptr<TH1D>> &histograms)
{
  std::unique_ptr<TH1D> average_hist = nullptr;
  int parameter_count = 0;

  try {
    ProcessDiagnosticFile(input_file, average_hist, parameter_count, histograms);
  } catch (const std::exception &e) {
    MACH3LOG_ERROR("Error processing file {} : {}", input_file, e.what());
  }

  if (average_hist && parameter_count > 0)
  {
    MACH3LOG_INFO("Processed {} parameters from {} files", parameter_count, input_file);
    average_hist->Scale(1.0 / parameter_count);
  }

  return average_hist;
}

void CompareAverageAC(const std::vector<std::vector<std::unique_ptr<TH1D>>> &histograms,
                      const std::vector<std::unique_ptr<TH1D>> &averages,
                      const std::vector<TString> &hist_labels,
                      const TString &output_name,
                      bool draw_min_max = true,
                      bool draw_all = false,
                      bool draw_errors = true)
{
  auto canvas = std::make_unique<TCanvas>("AverageAC", "Average Auto Correlation", 800, 600);
  canvas->SetGrid();

  // Setup colour palette
  constexpr Int_t nb = 255.0;
  Double_t stops[8] = {0.0, 0.25, 0.5, 0.75};
  Double_t red[8] = {0.83203125, 0.796875, 0.0, 0.9375};
  Double_t green[8] = {0.3671875, 0.47265625, 0.4453125, 0.890625};
  Double_t blue[8] = {0.0, 0.65234375, 0.6953125, 0.2578125};
  TColor::CreateGradientColorTable(8, stops, red, green, blue, nb);

  auto leg = std::make_unique<TLegend>(0.5, 0.7, 0.9, 0.9);
  leg->SetFillColorAlpha(kWhite, 0.7f);
  if (draw_min_max)
  {
    for (size_t i = 0; i < histograms.size(); ++i)
    {
      auto colour = static_cast<Color_t>(TColor::GetColorPalette(static_cast<Int_t>(i * nb / histograms.size())));
      auto band = CalculateMinMaxBand(histograms[i], colour);
      band->SetLineWidth(0);
      band->SetTitle("Average Auto Correlation");
      band->GetXaxis()->SetTitle("lag");
      band->GetYaxis()->SetTitle("Autocorrelation Function");
      if (i == 0) {
        band->Draw("A3");
      } else {
        band->Draw("3 SAME");
      }
    }
  }
  std::vector<std::unique_ptr<TH1D>> error_hist(averages.size());
  for (size_t i = 0; i < averages.size(); ++i)
  {
    auto colour = static_cast<Color_t>(TColor::GetColorPalette(static_cast<Int_t>(i * nb / histograms.size())));

    if (draw_errors)
    {
      error_hist[i] = M3::Clone(averages[i].get());
      error_hist[i]->SetFillColorAlpha(colour, 0.3f);
      error_hist[i]->SetLineWidth(0);
      error_hist[i]->SetTitle("Average Auto Correlation");
      error_hist[i]->GetXaxis()->SetTitle("lag");
      error_hist[i]->GetYaxis()->SetTitle("Autocorrelation Function");
      error_hist[i]->Draw("E3 SAME");
    }
    if (draw_all)
    {
      for (const auto &hist : histograms[i])
      {
        hist->SetLineColorAlpha(colour, 0.05f);
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

    TString int_str = TString::Format("Average integrated autocorrelation %.2f", averages[i]->Integral());
    leg->AddEntry(averages[i].get(), hist_labels[i] + int_str, "l");
  }

  if (averages.size() > 1)
  {
    leg->Draw();
  }

  canvas->SaveAs(output_name + ".pdf");
}

/// @brief Calculate mean AC based on all parameters with error band. Great for comparing AC between different chains
void PlotAverageACMult(std::vector<TString> input_files,
                       std::vector<TString> hist_labels,
                       const TString &output_name,
                       bool draw_min_max = true,
                       bool draw_all = false,
                       bool draw_errors = true)
{
  // Process first folder
  std::vector<std::vector<std::unique_ptr<TH1D>>> histograms;
  std::vector<std::unique_ptr<TH1D>> averages;

  for (int i = 0; i < static_cast<int>(input_files.size()); i++)
  {
    TString folder = input_files[i];
    MACH3LOG_INFO("Folder : {}", folder);

    if (folder.IsNull())
    {
      MACH3LOG_WARN("Skipping empty or dummy folder: {}", folder.Data());
      continue;
    }

    std::vector<std::unique_ptr<TH1D>> histograms_i;
    averages.emplace_back((AutocorrProcessInputs(folder, histograms_i)));
    histograms.push_back(std::move(histograms_i));
  }

  CompareAverageAC(histograms, averages, hist_labels, output_name, draw_min_max, draw_all, draw_errors);
}

int main(int argc, char *argv[])
{
  SetMaCh3LoggerFormat();
  // Require at least one file-title pair
  if (argc < 3 || (argc - 1) % 2 != 0)
  {
    MACH3LOG_ERROR("Wrong number of arguments ({}) provided", argc);
    MACH3LOG_ERROR("Usage: {} file1 title1 [file2 title2 ... fileN titleN]", argv[0]);
    MACH3LOG_ERROR("Example: {} DiagMCMC_Output.root PlotName", argv[0]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  std::vector<TString> filename;
  std::vector<TString> title;

  std::vector<TString> filenames;
  std::vector<TString> titles;

  // Parse all pairs (filename, title)
  for (int i = 1; i < argc; i += 2)
  {
    filenames.emplace_back(argv[i]);
    titles.emplace_back(argv[i + 1]);
  }

  // Log how many pairs we found
  MACH3LOG_INFO("Processing {} file(s)", filenames.size());

  // Generate plots
  MakeDiagPlot(filenames, titles);
  PlotAutoCorr(filenames);
  PlotAverageACMult(filenames, titles, Form("%s_Average_Auto_Corr", argv[2]), true);
  return 0;
}
