//MaCh3 Includes
#include "PlottingUtils/PlottingUtils.h"
#include "PlottingUtils/PlottingManager.h"

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file MatrixPlotter.cpp
/// @ingroup MaCh3Plotting
/// @author Kamil Skwarczynski

/// @brief Grab large Matrix only only extract submatrix based on label naming
std::unique_ptr<TH2D> GetSubMatrix(TH2D *MatrixFull,
                                   const std::string& Title,
                                   const std::vector<std::string>& Params,
                                   const std::unique_ptr<MaCh3Plotting::PlottingManager>& man)
{
  std::vector<int> ParamIndex(Params.size(), M3::_BAD_INT_);

  for(size_t i = 0; i < Params.size(); i++)
  {
    for(int j = 0; j < MatrixFull->GetNbinsX(); j++)
    {
      if(MatrixFull->GetXaxis()->GetBinLabel(j+1) == Params[i])
      {
        ParamIndex[i] = j;
        break;
      }
    }
  }

  for(size_t i = 0; i < Params.size(); i++)
  {
    if(ParamIndex[i] == M3::_BAD_INT_)
      MACH3LOG_ERROR("Didn't find param {} in matrix within {} sub-block", Params[i], Title);
  }

  auto new_end = std::remove(ParamIndex.begin(), ParamIndex.end(), M3::_BAD_INT_);
  ParamIndex.erase(new_end, ParamIndex.end());

  auto Hist = std::make_unique<TH2D>(Title.c_str(), Title.c_str(), ParamIndex.size(), 0, ParamIndex.size(), ParamIndex.size(), 0, ParamIndex.size());
  Hist->SetDirectory(nullptr);
  Hist->GetZaxis()->SetTitle("Correlation");
  Hist->SetMinimum(-1.);
  Hist->SetMaximum(1.);

  for(size_t x = 0; x < ParamIndex.size(); x++)
  {
    for(size_t y = 0; y < ParamIndex.size(); y++)
    {
      Hist->SetBinContent(x+1, y+1, MatrixFull->GetBinContent(ParamIndex[x]+1, ParamIndex[y]+1));
    }

    std::string FancyLabel = man->style().prettifyParamName(MatrixFull->GetXaxis()->GetBinLabel(ParamIndex[x]+1));
    Hist->GetXaxis()->SetBinLabel(x+1, FancyLabel.c_str());
    Hist->GetYaxis()->SetBinLabel(x+1, FancyLabel.c_str());
  }
  return Hist;
}

void DynamicLabelSize(TH2D* Hist) {
  if (Hist->GetNbinsX() < 20) {
    Hist->SetMarkerSize(1.0);
    Hist->GetXaxis()->SetLabelSize(0.02);
    Hist->GetYaxis()->SetLabelSize(0.02);
  } else {
    Hist->SetMarkerSize(0.5);
    Hist->GetXaxis()->SetLabelSize(0.015);
    Hist->GetYaxis()->SetLabelSize(0.015);
  }
}

void SetupInfo(const std::string& Config, std::vector<std::string>& Title, std::vector<std::vector<std::string>>& Params)
{
  // Load the YAML file
  YAML::Node config = M3OpenConfig(Config);

  // Access the "MatrixPlotter" section
  YAML::Node settings = config["MatrixPlotter"];

  // Retrieve the list of Titles
  Title = settings["Titles"].as<std::vector<std::string>>();
  Params.resize(Title.size());

  // Retrieve parameters for each title (using the titles as keys)
  for(size_t it = 0; it < Title.size(); it++)
  {
    if (settings[Title[it]]) {
      Params[it] = settings[Title[it]].as<std::vector<std::string>>();
    } else {
      MACH3LOG_ERROR("Missing key in YAML for: {}", Title[it]);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }

  // Check if sizes match
  if(Title.size() != Params.size())
  {
    MACH3LOG_ERROR("Size doesn't match");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
}

void PlotMatrix(const std::unique_ptr<MaCh3Plotting::PlottingManager>& man, const std::string& Config, const std::string& File)
{
  // Open the ROOT file
  TFile *file = M3::Open(File, "UPDATE", __FILE__, __LINE__);
  TH2D *MatrixFull = nullptr;
  file->GetObject("Correlation_plot", MatrixFull);

  if (!MatrixFull) {
    MACH3LOG_ERROR("Error: Could not retrieve histogram");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  auto MatrixPlot = std::make_unique<TCanvas>("MatrixPlot", "MatrixPlot", 0, 0, 1024, 1024);
  MatrixPlot->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  MatrixPlot->SetTickx();
  MatrixPlot->SetTicky();
  MatrixPlot->SetBottomMargin(0.2);
  MatrixPlot->SetTopMargin(0.1);
  MatrixPlot->SetRightMargin(0.15);
  MatrixPlot->SetLeftMargin(0.15);
  gStyle->SetOptTitle(1);
  gStyle->SetPaintTextFormat("4.1f");

  // Make pretty Correlation colors (red to blue)
  constexpr int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;
  DynamicLabelSize(MatrixFull);
  MatrixFull->GetXaxis()->LabelsOption("v");
  MatrixPlot->Print("MatrixPlot.pdf[");
  MatrixFull->SetTitle("");
  MatrixFull->Draw("COLZ");
  MatrixPlot->Print("MatrixPlot.pdf");
  std::vector<std::string> Title;
  std::vector<std::vector<std::string>> Params;

  SetupInfo(Config, Title, Params);

  for(size_t it = 0; it < Title.size(); it++)
  {
    std::unique_ptr<TH2D> Hist = GetSubMatrix(MatrixFull, Title[it], Params[it], man);
    Hist->GetXaxis()->LabelsOption("v");
    DynamicLabelSize(Hist.get());

    if(Hist->GetNbinsX() < 50) {
      Hist->Draw("COLZ TEXT");
    } else {
      Hist->Draw("COLZ");
    }
    MatrixPlot->Print("MatrixPlot.pdf");
  }

  MatrixPlot->Print("MatrixPlot.pdf]");

  delete MatrixFull;
  file->Close();
  delete file;
}

void CompareMatrices(const std::unique_ptr<MaCh3Plotting::PlottingManager>& man,
                     const std::string& Config, const std::string& File1, const std::string& Title1,
                     const std::string& File2, const std::string& Title2)
{
  (void) man;
  // Open the ROOT file
  constexpr int NFiles = 2;
  TFile *file[NFiles];
  file[0] = M3::Open(File1, "UPDATE", __FILE__, __LINE__);
  file[1] = M3::Open(File2, "UPDATE", __FILE__, __LINE__);

  TH2D *MatrixFull[NFiles] = {nullptr};
  for(int i = 0; i < NFiles; i++) {
    file[i]->GetObject("Correlation_plot", MatrixFull[i]);
  }
  auto MatrixPlot = std::make_unique<TCanvas>("MatrixPlot", "MatrixPlot", 0, 0, 1024, 1024);
  MatrixPlot->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  MatrixPlot->SetTickx();
  MatrixPlot->SetTicky();
  MatrixPlot->SetBottomMargin(0.2);
  MatrixPlot->SetTopMargin(0.1);
  MatrixPlot->SetRightMargin(0.15);
  MatrixPlot->SetLeftMargin(0.15);
  gStyle->SetOptTitle(1);

  //KS: Fancy colors
  constexpr int NRGBs = 10;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.10, 0.25, 0.35, 0.50, 0.60, 0.65, 0.75, 0.90, 1.00 };
  Double_t red[NRGBs]   = { 0.50, 1.00, 1.00, 0.25, 0.00, 0.10, 0.50, 1.00, 0.75, 0.55 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00, 0.60, 0.90, 1.00, 0.75, 0.75 };
  Double_t blue[NRGBs]  = { 0.00, 0.25, 1.00, 1.00, 0.50, 0.60, 0.90, 1.00, 0.05, 0.05 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  MatrixPlot->Print("MatrixComparePlot.pdf[");

  std::vector<std::string> Title;
  std::vector<std::vector<std::string>> Params;

  SetupInfo(Config, Title, Params);

  for(size_t it = 0; it < Title.size(); it++)
  {
    std::unique_ptr<TH2D> Hist[2];
    for(int i = 0; i < NFiles; i++)
    {
      Hist[i] = GetSubMatrix(MatrixFull[i], Title[it], Params[it], man);
    }
    Hist[0]->GetZaxis()->SetTitle( (Title1 + "/" +  Title2).c_str());
    Hist[0]->GetXaxis()->LabelsOption("v");
    Hist[0]->Divide(Hist[1].get());

    DynamicLabelSize(Hist[0].get());

    Hist[0]->Draw("COLZ");
    MatrixPlot->Print("MatrixComparePlot.pdf");
  }
  MatrixPlot->Print("MatrixComparePlot.pdf]");

  for(int i = 0; i < NFiles; i++)
  {
    delete MatrixFull[i];
    file[i]->Close();
    delete file[i];
  }
}

int main(int argc, char *argv[]) 
{
  SetMaCh3LoggerFormat();

  auto man = std::make_unique<MaCh3Plotting::PlottingManager>();
  man->initialise();

  if (argc != 3 && argc != 6)
  {
    MACH3LOG_INFO("How to use: {} config.yaml MCMC_Processor_Output.root", argv[0]);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (argc == 3)
  {
    PlotMatrix(man, std::string(argv[1]), std::string(argv[2]));
  }

  if (argc == 6)
  {
    MACH3LOG_INFO("Comparing matrices");
    CompareMatrices(man, std::string(argv[1]), std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), std::string(argv[5]));
  }

  MACH3LOG_INFO("Finished plotting matrices");
  return 0;
}
