//MaCh3 Includes
#include "plottingUtils/plottingUtils.h"
#include "plottingUtils/plottingManager.h"

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

TH2D* GetSubMatrix(TH2D *MatrixFull, const std::string& Title, const std::vector<std::string>& Params)
{
  std::vector<int> ParamIndex(Params.size(), -999);

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
    if(ParamIndex[i] == -999 )
      MACH3LOG_ERROR("Didn't find param {} in matrix within {} sub-block", Params[i], Title);
  }

  auto new_end = std::remove(ParamIndex.begin(), ParamIndex.end(), -999);
  ParamIndex.erase(new_end, ParamIndex.end());

  TH2D* Hist = new TH2D(Title.c_str(), Title.c_str(), ParamIndex.size(), 0, ParamIndex.size(), ParamIndex.size(), 0, ParamIndex.size());
  Hist->GetZaxis()->SetTitle("Correlation");
  Hist->SetMinimum(-1);
  Hist->SetMaximum(1);
  Hist->GetXaxis()->SetLabelSize(0.015);
  Hist->GetYaxis()->SetLabelSize(0.015);

  for(size_t x = 0; x < ParamIndex.size(); x++)
  {
    for(size_t y = 0; y < ParamIndex.size(); y++)
    {
      Hist->SetBinContent(x+1, y+1, MatrixFull->GetBinContent(ParamIndex[x]+1, ParamIndex[y]+1));
    }
    Hist->GetXaxis()->SetBinLabel(x+1, MatrixFull->GetXaxis()->GetBinLabel(ParamIndex[x]+1));
    Hist->GetYaxis()->SetBinLabel(x+1, MatrixFull->GetXaxis()->GetBinLabel(ParamIndex[x]+1));
  }
  return Hist;
}

void SetupInfo(const std::string& Config, std::vector<std::string>& Title, std::vector<std::vector<std::string>>& Params)
{
  // Load the YAML file
  YAML::Node config = YAML::LoadFile(Config);

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

void PlotMatrix(std::string Config, std::string File)
{
  // Open the ROOT file
  TFile *file = TFile::Open(File.c_str());
  if (!file || file->IsZombie()) {
    MACH3LOG_ERROR("Error opening file");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

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
  const int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  MatrixPlot->Print("MatrixPlot.pdf[");

  MatrixFull->Draw("COLZ");
  MatrixPlot->Print("MatrixPlot.pdf");
  std::vector<std::string> Title;
  std::vector<std::vector<std::string>> Params;

  SetupInfo(Config, Title, Params);

  for(size_t it = 0; it < Title.size(); it++)
  {
    TH2D* Hist = GetSubMatrix(MatrixFull, Title[it], Params[it]);
    Hist->GetXaxis()->LabelsOption("v");
    Hist->Draw("COLZ TEXT");
    MatrixPlot->Print("MatrixPlot.pdf");
    delete Hist;
  }
  MatrixPlot->Print("MatrixPlot.pdf]");

  delete MatrixFull;
  file->Close();
  delete file;
}

void CompareMatrices(std::string Config, std::string File1, std::string Title1, std::string File2, std::string Title2)
{
  // Open the ROOT file
  const int NFiles = 2;
  TFile *file[NFiles];
  file[0] = TFile::Open(File1.c_str());
  file[1] = TFile::Open(File2.c_str());

  TH2D *MatrixFull[2] = {nullptr};
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

  //KS: Fancy colots
  const int NRGBs = 10;
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
    TH2D* Hist[2] = {nullptr};
    for(int i = 0; i < NFiles; i++)
    {
      Hist[i] = GetSubMatrix(MatrixFull[i], Title[it], Params[it]);
    }
    Hist[0]->GetZaxis()->SetTitle( (Title1 + "/" +  Title2).c_str());
    Hist[0]->GetXaxis()->LabelsOption("v");

    Hist[0]->Draw("COLZ");
    MatrixPlot->Print("MatrixComparePlot.pdf");

    for(int i = 0; i < NFiles; i++)
    {
      delete Hist[i];
    }
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

  if (argc != 3 && argc != 6)
  {
    MACH3LOG_INFO("How to use: {} MCMC_Processor_Output.root config", argv[0]);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (argc == 3)
  {
    PlotMatrix(std::string(argv[1]), std::string(argv[2]));
  }

  if (argc == 6)
  {
    MACH3LOG_INFO("Comparing matrices");
    CompareMatrices(std::string(argv[1]), std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), std::string(argv[5]));
  }

  MACH3LOG_INFO("Finished plotting matrices");
  return 0;
}

