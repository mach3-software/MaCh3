//MaCh3 Includes
#include "plottingUtils/plottingUtils.h"
#include "plottingUtils/plottingManager.h"

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file PredictivePlotting.cpp
/// @author Kamil Skwarczynski

std::vector<std::string> FindSamples(const std::string& File)
{
  TFile *file = M3::Open(File, "READ", __FILE__, __LINE__);
  TDirectoryFile *PredicitveDir = file->Get<TDirectoryFile>("Predictive");

  std::vector<std::string> SampleNames;
  //Get all entries in input file
  TIter next(PredicitveDir->GetListOfKeys());
  TKey *key = nullptr;

  // Loop through all entries
  while ((key = static_cast<TKey*>(next()))) {
    // get directory names, ignore flux
    auto classname = std::string(key->GetClassName());
    auto dirname = std::string(key->GetName());

    if (classname != "TDirectoryFile") continue;
    dirname = std::string(key->GetName());

    if(dirname == "Total") continue;

    SampleNames.push_back(dirname);
    MACH3LOG_DEBUG("Entering Sample {}", dirname);
  }

  file->Close();
  delete file;
  return SampleNames;
}

void OverlayPredicitve(const YAML::Node& Settings,
                       const std::vector<TFile*>& InputFiles,
                       const std::vector<std::string>& SampleNames,
                       const std::unique_ptr<TCanvas>& canv)
{
  canv->Print("Overlay_Predictive.pdf[", "pdf");

  MACH3LOG_INFO("Starting {}", __func__);
  canv->Clear();

  TPad* pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->AppendPad();
  TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.25);
  pad2->AppendPad();

  pad1->SetGrid();
  pad2->SetGrid();

  pad1->SetLeftMargin(canv->GetLeftMargin());
  pad1->SetRightMargin(canv->GetRightMargin());
  pad1->SetTopMargin(canv->GetTopMargin());
  pad1->SetBottomMargin(0);

  pad2->SetLeftMargin(canv->GetLeftMargin());
  pad2->SetRightMargin(canv->GetRightMargin());
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.28);

  auto PosteriorColor = Get<std::vector<Color_t >>(Settings["PosteriorColor"], __FILE__, __LINE__);
  auto Titles = Get<std::vector<std::string>>(Settings["FileTitle"], __FILE__, __LINE__);

  if(Titles.size() < InputFiles.size() || PosteriorColor.size() < InputFiles.size()){
    MACH3LOG_ERROR("Passed {} files, while only {} titles and {} colors", InputFiles.size(), Titles.size(), PosteriorColor.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for(size_t iSample = 0; iSample < SampleNames.size(); iSample++)
  {
    const int nFiles = static_cast<int>(InputFiles.size());
    TH1D* hist = InputFiles[0]->Get<TH1D>(("SampleFolder/data_" + SampleNames[iSample]).c_str());
    if(!hist) {
      MACH3LOG_WARN("Couldn't find hist for {}, most likely it is using 2D", SampleNames[iSample]);
      MACH3LOG_WARN("Currently only 1D, sorry");
      continue;
    }
    std::unique_ptr<TH1D> DataHist = M3::Clone(hist);
    DataHist->SetLineColor(kBlack);
    //KS: +1 for data, we want to get integral before scaling of the histogram
    std::vector<double> Integral(nFiles);
    Integral[nFiles] = DataHist->Integral();
    std::vector<std::unique_ptr<TH1D>> PredHist(nFiles);

    for(int iFile = 0; iFile < nFiles; iFile++)
    {
      InputFiles[iFile]->cd();
      PredHist[iFile] = M3::Clone(InputFiles[iFile]->Get<TH1D>(("Predictive/" + SampleNames[iSample] + "/" +
                                                                SampleNames[iSample] + "_mc_PostPred").c_str()));
      Integral[iFile] = PredHist[iFile]->Integral();
      PredHist[iFile]->SetTitle(SampleNames[iSample].c_str());
      PredHist[iFile]->SetLineColor(PosteriorColor[iFile]);
      PredHist[iFile]->SetMarkerColor(PosteriorColor[iFile]);
      PredHist[iFile]->SetFillColorAlpha(PosteriorColor[iFile], 0.35);
      PredHist[iFile]->SetFillStyle(1001);
      PredHist[iFile]->GetYaxis()->SetTitle("Events");
    }
    pad1->cd();

    PredHist[0]->Draw("p e2");
    for(int iFile = 1; iFile < nFiles; iFile++)
    {
      PredHist[iFile]->Draw("p e2 same");
    }
    DataHist->Draw("he same");

    auto legend = std::make_unique<TLegend>(0.50,0.52,0.90,0.88);
    legend->AddEntry(DataHist.get(), Form("Data, #int=%.0f", Integral[nFiles]),"le");
    for(int ig = 0; ig < nFiles; ig++ )
    {
      legend->AddEntry(PredHist[ig].get(), Form("%s, #int=%.2f", Titles[ig].c_str(), Integral[ig]), "lpf");
    }
    legend->SetLineStyle(0);
    legend->SetTextSize(0.03);
    legend->Draw();

    //// Now we do ratio
    pad2->cd();

    auto line = std::make_unique<TLine>(PredHist[0]->GetXaxis()->GetBinLowEdge(PredHist[0]->GetXaxis()->GetFirst()), 1.0, PredHist[0]->GetXaxis()->GetBinUpEdge(PredHist[0]->GetXaxis()->GetLast()), 1.0);

    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("");

    std::unique_ptr<TH1D> RatioPlotData = M3::Clone(DataHist.get());
    std::vector<std::unique_ptr<TH1D>> RatioPlot(nFiles);

    for(int ig = 0; ig < nFiles; ig++ )
    {
      RatioPlot[ig] = M3::Clone(DataHist.get());
      RatioPlot[ig]->SetLineColor(PosteriorColor[ig]);
      RatioPlot[ig]->SetMarkerColor(PosteriorColor[ig]);
      RatioPlot[ig]->SetFillColorAlpha(PosteriorColor[ig], 0.35);
      RatioPlot[ig]->SetFillStyle(1001);
      RatioPlot[ig]->GetYaxis()->SetTitle("Data/MC");
      RatioPlot[ig]->GetXaxis()->SetTitle(PredHist[0]->GetXaxis()->GetTitle());
      RatioPlot[ig]->SetBit(TH1D::kNoTitle);
      RatioPlot[ig]->GetXaxis()->SetTitleSize(0.12);
      RatioPlot[ig]->GetYaxis()->SetTitleOffset(0.4);
      RatioPlot[ig]->GetYaxis()->SetTitleSize(0.10);

      RatioPlot[ig]->GetXaxis()->SetLabelSize(0.10);
      RatioPlot[ig]->GetYaxis()->SetLabelSize(0.10);

      RatioPlot[ig]->Divide(PredHist[ig].get());
      PassErrorToRatioPlot(RatioPlot[ig].get(), PredHist[ig].get(), DataHist.get());
    }

    RatioPlotData->Divide(DataHist.get());
    PassErrorToRatioPlot(RatioPlotData.get(), DataHist.get(), DataHist.get());

    double maxz = -999;
    double minz = +999;
    for (int j = 0; j < nFiles; j++) {
      for (int i = 1; i < RatioPlot[0]->GetXaxis()->GetNbins(); i++)
      {
        maxz = std::max(maxz, RatioPlot[j]->GetBinContent(i));
        minz = std::min(minz, RatioPlot[j]->GetBinContent(i));
      }
    }
    maxz = maxz*1.001;
    minz = minz*1.001;

    if (std::fabs(1 - maxz) > std::fabs(1-minz))
      RatioPlot[0]->GetYaxis()->SetRangeUser(1-std::fabs(1-maxz),1+std::fabs(1-maxz));
    else
      RatioPlot[0]->GetYaxis()->SetRangeUser(1-std::fabs(1-minz),1+std::fabs(1-minz));

    RatioPlot[0]->Draw("p e2");
    for(int ig = 1; ig < nFiles; ig++ )
    {
      RatioPlot[ig]->Draw("p e2 same");
    }
    RatioPlotData->Draw("he same");

    canv->Print("Overlay_Predictive.pdf", "pdf");
  }

  canv->Print("Overlay_Predictive.pdf]", "pdf");
  delete pad1;
  delete pad2;
}

void PredictivePlotting(const std::string& ConfigName,
                        const std::vector<std::string>& FileNames)
{
  auto canvas = std::make_unique<TCanvas>("canv", "canv", 1080, 1080);
  // set the paper & margin sizes
  canvas->SetTopMargin(0.11);
  canvas->SetBottomMargin(0.16);
  canvas->SetRightMargin(0.075);
  canvas->SetLeftMargin(0.12);
  canvas->SetGrid();

  gStyle->SetOptStat(0);  //Set 0 to disable statystic box
  gStyle->SetPalette(51);
  gStyle->SetLegendBorderSize(0); //This option disables legends borders
  gStyle->SetFillStyle(0);

  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  auto Samples = FindSamples(FileNames[0]);
  std::vector<TFile*> InputFiles(FileNames.size());
  for(size_t i = 0; i < FileNames.size(); i++)
  {
    InputFiles[i] = M3::Open(FileNames[i], "READ", __FILE__, __LINE__);
  }

  // Load the YAML file
  YAML::Node Config = M3OpenConfig(ConfigName);
  // Access the "MatrixPlotter" section
  YAML::Node settings = Config["PredictivePlotting"];

  OverlayPredicitve(settings, InputFiles, Samples, canvas);

  for(size_t i = 0; i < FileNames.size(); i++)
  {
    InputFiles[i]->Close();
    delete InputFiles[i];
  }
}


int main(int argc, char **argv)
{
  SetMaCh3LoggerFormat();
  if (argc < 3)
  {
    MACH3LOG_ERROR("Need at least two arguments, {} <Config.Yaml> <Prior/Post_PredOutput.root>", argv[0]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  std::string ConfigName = std::string(argv[1]);
  // Collect all remaining arguments as file names
  std::vector<std::string> FileNames;
  for (int i = 2; i < argc; ++i) {
    FileNames.emplace_back(argv[i]);
  }

  PredictivePlotting(ConfigName, FileNames);

  return 0;
}
