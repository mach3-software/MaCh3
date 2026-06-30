//MaCh3 Includes
#include "PlottingUtils/PlottingUtils.h"
#include "PlottingUtils/PlottingManager.h"
#include <numeric>

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file PredictivePlotting.cpp
/// @ingroup MaCh3Plotting
/// @author Kamil Skwarczynski

/// @warning KS: keep raw pointer or ensure manual delete of PlotMan. If spdlog in automatically deleted before PlotMan then destructor has some spdlog and this could cause segfault
M3::Plotting::PlottingManager* PlotMan = nullptr;

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
    if(dirname == "BetaParameters") continue;
    if(dirname == "Correlations") continue;

    SampleNames.push_back(dirname);
    MACH3LOG_DEBUG("Entering Sample {}", dirname);
  }

  file->Close();
  delete file;
  return SampleNames;
}

std::vector<int> FindDimensions(const std::string& File, const std::vector<std::string>& Samples)
{
  TFile *file = M3::Open(File, "READ", __FILE__, __LINE__);
  TDirectoryFile *PredicitveDir = file->Get<TDirectoryFile>("Predictive");

  std::vector<int> SampleDimension;
  for (const auto& sample : Samples)
  {
    // Get directory for this sample
    TDirectoryFile* SampleDir = PredicitveDir->Get<TDirectoryFile>(sample.c_str());

    int Dimension = 0;

    while (true)
    {
      // Construct name Tutorial_mc_dimX
      std::string histName = fmt::format("{}_mc_dim{}", sample, Dimension);

      TH2D* hist = SampleDir->Get<TH2D>(histName.c_str());
      if (!hist) break;  // stop when next dimension does not exist

      Dimension++;
    }

    MACH3LOG_DEBUG("Sample '{}' has dimension {}", sample, Dimension);
    SampleDimension.push_back(Dimension);
  }

  file->Close();
  delete file;

  return SampleDimension;
}


std::vector<std::vector<std::string>> FindModes(const std::string& File,
                                                const std::vector<std::string>& SampleNames)
{
  TFile *file = M3::Open(File, "READ", __FILE__, __LINE__);
  TDirectoryFile *PredictiveDir = file->Get<TDirectoryFile>("Predictive");

  std::vector<std::vector<std::string>> ModeNames(SampleNames.size());

  for(size_t iSample = 0; iSample < SampleNames.size(); iSample++)
  {
    TDirectoryFile* SampleDir = PredictiveDir->Get<TDirectoryFile>(SampleNames[iSample].c_str());
    if(!SampleDir) continue;

    TDirectoryFile* ByModeDir = SampleDir->Get<TDirectoryFile>("ByMode");
    if(!ByModeDir) continue;

    // Loop over all keys in ByModeDir
    TIter next(ByModeDir->GetListOfKeys());
    TKey* key;

    while ((key = static_cast<TKey*>(next())))
    {
      TObject* obj = key->ReadObj();
      if (!obj->InheritsFrom("TH1")) continue;

      std::string histName = obj->GetName();

      // Example: sample_mode_dim0 → extract "mode"
      std::string prefix = SampleNames[iSample] + "_";
      std::string suffix = "_dim0";

      if (histName.find(prefix) == 0 &&
        histName.rfind(suffix) == histName.size() - suffix.size())
      {
        std::string Mode = histName.substr(prefix.size(),
                                           histName.size() - prefix.size() - suffix.size());
        MACH3LOG_DEBUG("Found mode '{}' for sample {}", Mode, SampleNames[iSample]);
        ModeNames[iSample].push_back(Mode);
      }
    }
  }

  file->Close();
  delete file;

  return ModeNames;
}


void PretifyHistogram(TH1* Hist, const std::string& SampleName) {
  Hist->SetTitle(PlotMan->style().prettifySampleName(SampleName).c_str());
  auto BinWidthScale = PlotMan->style().getBinWidthScale(Hist->GetXaxis()->GetTitle());
  auto PrettyX = PlotMan->style().prettifyKinematicName(Hist->GetXaxis()->GetTitle());
  Hist->GetXaxis()->SetTitle(PrettyX.c_str());
  Hist->GetYaxis()->SetTitle(fmt::format("Events/{:.0f}", BinWidthScale).c_str());
  M3::ScaleHistogram(Hist, BinWidthScale);
}

double GetPValue(const TH2D* hist)
{
  double pvalue = 0;
  std::string TempTile = hist->GetTitle();
  std::string temp = "=";

  std::string::size_type SizeType = TempTile.find(temp);
  TempTile.erase(0, SizeType+1);
  pvalue = atof(TempTile.c_str());
  return pvalue;
}

void PrintPosteriorPValue(const YAML::Node& Settings,
                          const std::vector<TFile*>& InputFiles,
                          const std::vector<std::string>& SampleNames)
{
  MACH3LOG_INFO("Starting {}", __func__);
  auto Titles = Get<std::vector<std::string>>(Settings["FileTitle"], __FILE__, __LINE__);
  std::vector<std::vector<double>> FlucDrawVec(InputFiles.size());
  // KS: Alternatively try "_drawfluc_draw"
  std::string FlucutationType = "_predfluc_draw";
  //KS: P-values per each sample
  std::cout<<"\\begin{table}[htb]"<<std::endl;
  std::cout<<"\\centering"<<std::endl;
  std::cout<<"\\begin{tabular}{ | l | ";

  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    std::cout<<"c | ";
  }

  std::cout<<"} \\hline"<<std::endl;
  std::cout<<"Sample ";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    std::cout<<"& \\multicolumn{1}{| c |}{" + Titles[f] +" p-value} ";
  }
  std::cout<<"\\\\"<<std::endl;
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    std::cout<<" & Fluctuation of Prediction ";
  }
  std::cout<<"\\\\ \\hline"<<std::endl;
  for(unsigned int i = 0; i < SampleNames.size(); i++)
  {
    std::cout<<SampleNames[i];
    for(unsigned int f = 0; f < InputFiles.size(); f++)
    {
      std::string TempString = "Predictive/" + SampleNames[i]+"/"+SampleNames[i] + FlucutationType;
      TH2D *hist2D = InputFiles[f]->Get<TH2D>(TempString.c_str());
      double FlucDraw = GetPValue(hist2D);
      std::cout<<" & "<<FlucDraw;
      FlucDrawVec[f].push_back(FlucDraw);
    }
    std::cout<<" \\\\"<<std::endl;
  }
  std::cout<<"Total ";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    TH2D *hFlucPred = InputFiles[f]->Get<TH2D>(("Predictive/Total/Total" + FlucutationType).c_str());
    double FlucDraw = GetPValue(hFlucPred);
    std::cout<<" & "<<FlucDraw;
  }
  std::cout<<" \\\\ \\hline"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{table}"<<std::endl;

  auto Threshold = GetFromManager<double>(Settings["Significance"], 0.05, __FILE__ , __LINE__);
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    MACH3LOG_INFO("Calculating Shape for file {}", Titles[f]);

    CheckBonferoniCorrectedpValue(SampleNames, FlucDrawVec[f], Threshold);
    MACH3LOG_INFO("Combined pvalue following Fisher method: {:.4f}", FisherCombinedPValue(FlucDrawVec[f]));
  }
}

void OverlayViolin(const YAML::Node& Settings,
                   const std::vector<TFile*>& InputFiles,
                   const std::vector<std::string>& SampleNames,
                   const std::vector<int>& SampleDimension,
                   const std::unique_ptr<TCanvas>& canv)
{
  MACH3LOG_INFO("Starting {}", __func__);
  canv->Clear();

  canv->SetTopMargin(0.10);
  canv->SetBottomMargin(0.12);
  canv->SetRightMargin(0.075);
  canv->SetLeftMargin(0.14);

  auto PosteriorColor = Get<std::vector<Color_t >>(Settings["PosteriorColor"], __FILE__, __LINE__);
  auto Titles = Get<std::vector<std::string>>(Settings["FileTitle"], __FILE__, __LINE__);
  const int nFiles = static_cast<int>(InputFiles.size());

  //KS: No idea why but ROOT changed treatment of violin in R6. If you have non uniform binning this will results in very hard to see violin plots.
  TCandle::SetScaledViolin(false);
  for(size_t iSample = 0; iSample < SampleNames.size(); iSample++)
  {
    for(int iDim = 0; iDim < SampleDimension[iSample]; iDim++)
    {
      std::vector<std::unique_ptr<TH2D>> ViolinHist(nFiles);
      for(int iFile = 0; iFile < nFiles; iFile++)
      {
        InputFiles[iFile]->cd();
        ViolinHist[iFile] = M3::Clone(InputFiles[iFile]->Get<TH2D>(("Predictive/" + SampleNames[iSample]
                                      + "/" + SampleNames[iSample] + "_mc_dim" + iDim).Data()));
        ViolinHist[iFile]->SetTitle(PlotMan->style().prettifySampleName(SampleNames[iSample]).c_str());
        ViolinHist[iFile]->SetLineColor(PosteriorColor[iFile]);
        ViolinHist[iFile]->SetMarkerColor(PosteriorColor[iFile]);
        ViolinHist[iFile]->SetFillColorAlpha(PosteriorColor[iFile], 0.35);
        ViolinHist[iFile]->SetFillStyle(1001);
        ViolinHist[iFile]->GetXaxis()->SetTitle(PlotMan->style().prettifyKinematicName(
                                                ViolinHist[iFile]->GetXaxis()->GetTitle()).c_str());
        ViolinHist[iFile]->GetYaxis()->SetTitle("Events");
      }

      ViolinHist[0]->Draw("violinX(03100300)");
      for(int iFile = 1; iFile < nFiles; iFile++) {
        ViolinHist[iFile]->Draw("violinX(03100300) same");
      }

      TLegend legend(0.50, 0.52, 0.90, 0.88);
      for(int ig = 0; ig < nFiles; ig++) {
        legend.AddEntry(ViolinHist[ig].get(), Form("%s", Titles[ig].c_str()), "lpf");
      }
      legend.SetLineStyle(0);
      legend.SetTextSize(0.03);
      legend.Draw();

      canv->Print("Overlay_Predictive.pdf", "pdf");
    }
  }
}

void OverlayPredicitve(const YAML::Node& Settings,
                       const std::vector<TFile*>& InputFiles,
                       const std::vector<std::string>& SampleNames,
                       const std::vector<int>& SampleDimension,
                       const std::unique_ptr<TCanvas>& canv)
{
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
    auto SampleName = SampleNames[iSample];
    const int nDims = (SampleDimension[iSample] == 2) ? 2 : 1;
    for(int iDim = 0; iDim < nDims; iDim++) {
      std::string DataLocation = "";
      if(nDims == 2) {
        DataLocation = "Predictive/" + SampleName + "/Data_" + SampleName + "_Dim" + std::to_string(iDim);
      } else {
        DataLocation = "SampleFolder/data_" + SampleName;
      }
      TH1D* hist = InputFiles[0]->Get<TH1D>((DataLocation).c_str());

      auto BinWidthScale = PlotMan->style().getBinWidthScale(hist->GetXaxis()->GetTitle());
      std::unique_ptr<TH1D> DataHist = M3::Clone(hist);
      auto DataPoissonErrors = PoissonGraphScaled(DataHist.get(), BinWidthScale);
      M3::ScaleHistogram(DataHist.get(), BinWidthScale);

      DataHist->SetLineColor(kBlack);
      DataPoissonErrors->SetLineColor(kBlack);
      //KS: +1 for data, we want to get integral before scaling of the histogram
      std::vector<double> Integral(nFiles+1);
      Integral[nFiles] = DataHist->Integral();
      std::vector<std::unique_ptr<TH1D>> PredHist(nFiles);

      for(int iFile = 0; iFile < nFiles; iFile++)
      {
        InputFiles[iFile]->cd();
        std::string HistLocation = "";
        if(nDims == 2) {
          HistLocation = "Predictive/" + SampleName + "/" + SampleName + "_mc_PostPred_dim" + std::to_string(iDim);
        } else {
          HistLocation = "Predictive/" + SampleName + "/" + SampleName + "_mc_PostPred";
        }
        PredHist[iFile] = M3::Clone(InputFiles[iFile]->Get<TH1D>((HistLocation).c_str()));
        Integral[iFile] = PredHist[iFile]->Integral();
        PredHist[iFile]->SetLineColor(PosteriorColor[iFile]);
        PredHist[iFile]->SetMarkerColor(PosteriorColor[iFile]);
        PredHist[iFile]->SetFillColorAlpha(PosteriorColor[iFile], 0.35);
        PredHist[iFile]->SetFillStyle(1001);
        PretifyHistogram(PredHist[iFile].get(), SampleName);
      }
      pad1->cd();

      PredHist[0]->Draw("p e2");
      for(int iFile = 1; iFile < nFiles; iFile++) {
        PredHist[iFile]->Draw("p e2 same");
      }
      DataPoissonErrors->Draw("p same");

      auto legend = std::make_unique<TLegend>(0.50,0.52,0.90,0.88);
      legend->AddEntry(DataPoissonErrors.get(), Form("Data, #int=%.0f", Integral[nFiles]),"le");
      for(int ig = 0; ig < nFiles; ig++ ) {
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
        auto PrettyX = PlotMan->style().prettifyKinematicName(PredHist[0]->GetXaxis()->GetTitle());
        RatioPlot[ig]->GetXaxis()->SetTitle(PrettyX.c_str());
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

      M3::Plotting::SetSymmetricRatioRange(RatioPlot);

      RatioPlot[0]->Draw("p e2");
      for(int ig = 1; ig < nFiles; ig++ ) {
        RatioPlot[ig]->Draw("p e2 same");
      }
      RatioPlotData->Draw("he same");

      canv->Print("Overlay_Predictive.pdf", "pdf");
    }
  }

  delete pad1;
  delete pad2;
}

void OverlayPredicitveByMode(const YAML::Node& Settings,
                       const std::vector<TFile*>& InputFiles,
                       const std::vector<std::string>& SampleNames,
                       const std::vector<int>& SampleDimension,
                       const std::vector<std::vector<std::string>>& Modes,
                       const std::unique_ptr<TCanvas>& canv)
{
  MACH3LOG_INFO("Starting {}", __func__);
  canv->cd();
  constexpr auto DefaultColor = kBlack;
  auto Titles = Get<std::vector<std::string>>(Settings["FileTitle"], __FILE__, __LINE__);
  auto RelevantModesName = Get<std::vector<std::string>>(Settings["RelevantModesName"], __FILE__, __LINE__);
  auto RelevantColors = Get<std::vector<Color_t>>(Settings["RelevantModesColors"], __FILE__, __LINE__);
  int nRelevantModes = static_cast<int>(RelevantModesName.size());
  const int nFiles = static_cast<int>(InputFiles.size());
  if(Titles.size() < InputFiles.size()){
    MACH3LOG_ERROR("Passed {} files, while only {} titles", InputFiles.size(), Titles.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if(RelevantModesName.size() != RelevantColors.size()) {
    MACH3LOG_ERROR("Colors ({}) doesn't match relevant modes {}", RelevantColors.size(), RelevantModesName.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for(int iFile = 0; iFile < nFiles; iFile++ )
  {
    for(size_t iSample = 0; iSample < SampleNames.size(); iSample++)
    {
      auto SampleName = SampleNames[iSample];
      for(int iDim = 0; iDim < SampleDimension[iSample]; iDim++)
      {
        const int nDims = (SampleDimension[iSample] == 2) ? 2 : 1;
        std::string HistLocation = "";
        if(nDims == 2) {
          HistLocation = "Predictive/" + SampleName + "/" + SampleName + "_mc_PostPred_dim" + std::to_string(iDim);
        } else {
          HistLocation = "Predictive/" + SampleName + "/" + SampleName + "_mc_PostPred";
        }
        std::unique_ptr<TH1D> Sample_MC_Full = M3::Clone(InputFiles[iFile]->Get<TH1D>((HistLocation).c_str()));
        Sample_MC_Full->SetLineColor(kOrange);
        Sample_MC_Full->SetLineWidth(2);
        Sample_MC_Full->SetMarkerColor(kOrange);
        PretifyHistogram(Sample_MC_Full.get(), SampleName);

        std::string DataLocation = "";
        std::unique_ptr<TH1D> Sample_Data;
        if(nDims == 2) {
         DataLocation = "Predictive/" + SampleName + "/Data_" + SampleName + "_Dim" + std::to_string(iDim);
        } else if(nDims == 1) {
          DataLocation = "SampleFolder/data_" + SampleName;
        }
        if(DataLocation != "") {
          Sample_Data = M3::Clone(InputFiles[iFile]->Get<TH1D>((DataLocation).c_str()));
          Sample_Data->SetLineColor(kBlack);
          Sample_Data->SetLineWidth(2);
          Sample_Data->SetMarkerColor(kBlack);
          PretifyHistogram(Sample_Data.get(), SampleName);
        }
        int nModes = static_cast<int>(Modes[iSample].size());
        // Simple map to keep track which mode is relevant and which will be added to "Other"
        std::vector<bool> isRelevantMode(nModes, false);
        std::vector<Color_t > ColorMap(nModes, DefaultColor);
        for(int iMode = 0; iMode < nModes; iMode++) {
          for(int iRelevant = 0; iRelevant < nRelevantModes; iRelevant++) {
            if(Modes[iSample][iMode] == RelevantModesName[iRelevant]) {
              isRelevantMode[iMode] = true;
              ColorMap[iMode] = RelevantColors[iRelevant];
            }
          }
        }
        auto Sample_Stack = std::make_unique<THStack>(SampleName.c_str(), SampleName.c_str());
        // This will hold values for "Other" modes
        std::unique_ptr<TH1D> Sample_MC_Other;
        // KS: Store histogram for each mode
        std::vector<std::unique_ptr<TH1D>> Sample_MC(nModes);
        std::vector<double> Integrals(nModes, 0.);
        for(int iMode = 0; iMode < nModes; iMode++)
        {
          std::string FileLocaction = "Predictive/" + SampleName + "/ByMode/" + SampleName
                                    + "_" + Modes[iSample][iMode] + "_dim" + std::to_string(iDim);
          auto SpectraByMode = InputFiles[iFile]->Get<TH2D>((FileLocaction).c_str());
          if(SpectraByMode == nullptr){
            MACH3LOG_ERROR("Something went wrong and didn't find histogram: {}", FileLocaction);
            throw MaCh3Exception(__FILE__, __LINE__);
          }
          Sample_MC[iMode] = MakeSummaryFromSpectra(SpectraByMode, SpectraByMode->GetTitle());
          Integrals[iMode] = Sample_MC[iMode]->Integral();
          PretifyHistogram(Sample_MC[iMode].get(), SampleName);

          if(Sample_MC_Other == nullptr) {
            Sample_MC_Other = M3::Clone(Sample_MC[iMode].get());
            Sample_MC_Other->Reset();
            Sample_MC_Other->SetFillColor(DefaultColor);
            Sample_MC_Other->SetLineColor(DefaultColor);
          }
          if(!isRelevantMode[iMode]) {
            Sample_MC_Other->Add(Sample_MC[iMode].get());
          }
          Sample_MC[iMode]->SetFillColor(ColorMap[iMode]);
          Sample_MC[iMode]->SetLineColor(ColorMap[iMode]);
        } // end loop over modes
        Sample_Stack->Add(Sample_MC_Other.get());
        // KS: We do this other way around as we want to have most relevant modes first
        for(int iMode = nModes-1; iMode >= 0; iMode--) {
          if(isRelevantMode[iMode]) Sample_Stack->Add( Sample_MC[iMode].get() );
        }
        Sample_Stack->Draw("hist");
        Sample_MC_Full->Draw("SAME he");
        if(Sample_Data) Sample_Data->Draw("SAME pe");
        canv->cd();
        Sample_Stack->GetXaxis();
        Sample_Stack->SetTitle(Sample_MC_Other->GetTitle());
        Sample_Stack->GetXaxis()->SetTitle(Sample_MC_Other->GetXaxis()->GetTitle());
        Sample_Stack->GetYaxis()->SetTitle(Sample_MC_Other->GetYaxis()->GetTitle());

        double FullIntegral = std::accumulate(Integrals.begin(), Integrals.end(), 0.0);
        double OtherIntegral = 0.;
        TLegend legend(0.50,0.52,0.85,0.88);
        if(Sample_Data) legend.AddEntry(Sample_Data.get(), "Data","ple");
        legend.AddEntry(Sample_MC_Full.get(), Titles[iFile].c_str(),"fple");
        for(int iMode = 0; iMode < nModes; iMode++) {
          if(isRelevantMode[iMode]) {
            std::string Label = Form("%s: %.1f%%", Modes[iSample][iMode].c_str(), Integrals[iMode]/FullIntegral*100);
            legend.AddEntry(Sample_MC[iMode].get(), Label.c_str() ,"lf");
          } else{
            OtherIntegral += Integrals[iMode]/FullIntegral;
          }
        }
        legend.AddEntry(Sample_MC_Other.get(), Form("Other: %.1f%%", OtherIntegral*100), "lf");
        legend.SetTextSize(0.03);
        legend.Draw();

        canv->Print("Overlay_Predictive.pdf", "pdf");
      } // end loop over dimensions
    } // end loop over samples
  } // end loop over files
}


/// @brief KS: Get mean and error from gaussian fit to event distribution
void GetMeanError(TH1D* hist, double &Mean, double &Error){
  TF1 *Gauss = hist->GetFunction("Fit"); //This name is hardcoded be careful
  //KS: Get mean and error from Gauss
  Mean = Gauss->GetParameter(1);
  Error = Gauss->GetParameter(2);

  //KS: Get mean and error from HPD
  //Mean = hist->GetMean();
  //Error = hpost->GetRMS();
}

/// @brief KS Print event rates in Latex like table
void PrintPosteriorEventRates(const std::vector<TFile*>& InputFiles,
                              const std::vector<std::string>& SampleNames) {
  MACH3LOG_INFO("Starting {}", __func__);
  MACH3LOG_INFO("");

  double mean, error;
  //KS: We now prepare to make tables for TN etc.
  std::cout<<"\\begin{table}[htb]"<<std::endl;
  std::cout<<"\\centering"<<std::endl;
  std::cout<<"\\begin{tabular}{ | l |";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    std::cout<<" c |";
  }
  std::cout<<"} \\hline"<<std::endl;
  std::cout<<"Sample ";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    std::cout<<"& Event Rates  ";
  }
  std::cout<<"\\\\ \\hline"<<std::endl;
  for(unsigned int i = 0; i < SampleNames.size(); i++)
  {
    std::cout<<SampleNames[i];
    std::string TempString = "Predictive/" + SampleNames[i]+"/"+SampleNames[i]+"_sum";
    for(unsigned int f = 0; f < InputFiles.size(); f++)
    {
      TH1D *hist = static_cast<TH1D*>(InputFiles[f]->Get(TempString.c_str()));
      GetMeanError(hist, mean, error);
      std::cout<<" & "<<mean<<" $\\pm$ "<<error;
    }
    std::cout<<" \\\\"<<std::endl;
  }
  std::cout<<"Total";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    TH1D *histTot = static_cast<TH1D*>(InputFiles[f]->Get("Predictive/Total/Total_sum"));
    GetMeanError(histTot, mean, error);
    std::cout<<" & "<<mean<<" $\\pm$ "<<error;
  }
  std::cout<<" \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{table}"<<std::endl;
  MACH3LOG_INFO("");
}

/// @brief KS: Print Fractional Uncertainties into Latex table format
void PrintPosteriorFractionalUncertainties(const std::vector<TFile*>& InputFiles,
                                           const std::vector<std::string>& SampleNames) {
  MACH3LOG_INFO("Starting {}", __func__);
  MACH3LOG_INFO("");
  double mean, error;

  //KS: Fractional uncertainties on the prior and posterior predictive event rates.
  std::cout<<"\\begin{table}[htb]"<<std::endl;
  std::cout<<"\\centering"<<std::endl;
  std::cout<<"\\begin{tabular}{ | l |";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    std::cout<<" c |";
  }
  std::cout<<"} \\hline"<<std::endl;

  std::cout<<"Sample ";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    std::cout<<"& $\\delta N / N (\\%)$";
  }
  std::cout<<"\\\\ \\hline"<<std::endl;

  for(unsigned int i = 0; i < SampleNames.size(); i++)
  {
    std::cout<<SampleNames[i];
    std::string TempString = "Predictive/" + SampleNames[i]+"/"+SampleNames[i]+"_sum";
    for(unsigned int f = 0; f < InputFiles.size(); f++)
    {
      TH1D *hist = static_cast<TH1D*>(InputFiles[f]->Get(TempString.c_str()));
      GetMeanError(hist, mean, error);
      std::cout<<" & "<<error/mean*100;
    }
    std::cout<<" \\\\"<<std::endl;
  }
  std::cout<<"Total";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    TH1D *histTotal = static_cast<TH1D*>(InputFiles[f]->Get("Predictive/Total/Total_sum"));
    GetMeanError(histTotal, mean, error);
    std::cout<<" & "<<error/mean*100;
  }
  std::cout<<"\\\\ \\hline"<<std::endl;
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{table}"<<std::endl;
  MACH3LOG_INFO("");
}

double GetLLH(TH1* hist)
{
  std::string TempTile = hist->GetTitle();
  std::string temp = "=";

  std::string::size_type SizeType = TempTile.find(temp);
  TempTile.erase(0, SizeType+1);
  double llh = atof(TempTile.c_str());
  return llh;
}

/// @brief KS Print Predictive LLH into Latex table format
void PrintPredictiveLLH(const std::vector<TFile*>& InputFiles,
                        const std::vector<std::string>& SampleNames) {
  MACH3LOG_INFO("Starting {}", __func__);
  MACH3LOG_INFO("");

  std::vector<double> Total(InputFiles.size());
  //KS: We now prepare to make tables for TN etc.
  std::cout<<"\\begin{table}[htb]"<<std::endl;
  std::cout<<"\\centering"<<std::endl;
  std::cout<<"\\begin{tabular}{ | l |";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    Total[f] = 0.;
    std::cout<<" c |";
  }
  std::cout<<"} \\hline"<<std::endl;
  std::cout<<"Sample ";
  for(unsigned int f = 0; f < InputFiles.size(); f++)
  {
    std::cout<<"& 2#log#mathcal{L}_{stat}  ";
  }
  std::cout<<"\\\\ \\hline"<<std::endl;
  for(unsigned int i = 0; i < SampleNames.size(); i++)
  {
    std::cout<<SampleNames[i];
    std::string TempString = "Predictive/" + SampleNames[i]+"/"+SampleNames[i]+"_mc_PostPred";
    for(unsigned int f = 0; f < InputFiles.size(); f++)
    {
      TH1 *hist = static_cast<TH1*>(InputFiles[f]->Get(TempString.c_str()));

      double llh = GetLLH(hist);
      std::cout<<" & "<<llh;
      Total[f] += llh;
    }
    std::cout<<" \\\\"<<std::endl;
  }
  std::cout<<"Total";
  for(unsigned int f = 0; f < InputFiles.size(); f++) {
    std::cout<<" & "<<Total[f];
  }
  std::cout<<" \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{table}"<<std::endl;
  std::cout<<" "<<std::endl;
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

  gStyle->SetOptStat(0);  //Set 0 to disable statistic box
  gStyle->SetPalette(51);
  gStyle->SetLegendBorderSize(0); //This option disables legends borders
  gStyle->SetFillStyle(0);

  //To avoid TCanvas::Print> messages
  gErrorIgnoreLevel = kWarning;

  auto Samples = FindSamples(FileNames[0]);
  auto Dimensions = FindDimensions(FileNames[0], Samples);
  auto Modes = FindModes(FileNames[0], Samples);

  std::vector<TFile*> InputFiles(FileNames.size());
  for(size_t i = 0; i < FileNames.size(); i++) {
    InputFiles[i] = M3::Open(FileNames[i], "READ", __FILE__, __LINE__);
  }

  // Load the YAML file
  YAML::Node Config = M3OpenConfig(ConfigName);
  // Access the "MatrixPlotter" section
  YAML::Node settings = Config["PredictivePlotting"];
  canvas->Print("Overlay_Predictive.pdf[", "pdf");

  // Make overlay of 1D hists
  OverlayPredicitve(settings, InputFiles, Samples, Dimensions, canvas);
  // Make overlay of violin plots
  OverlayViolin(settings, InputFiles, Samples, Dimensions, canvas);
  // Make By Mode post pred
  if(Modes[0].size() != 0) OverlayPredicitveByMode(settings, InputFiles, Samples, Dimensions, Modes, canvas);
  // Get PValue per sample
  PrintPosteriorPValue(settings, InputFiles, Samples);
  // KS: Print Fractional Uncertainties into Latex table format
  PrintPosteriorEventRates(InputFiles, Samples);
  // KS: Print Fractional Uncertainties into Latex table format
  PrintPosteriorFractionalUncertainties(InputFiles, Samples);
  // KS: Print Predictive LLH into Latex table format
  PrintPredictiveLLH(InputFiles, Samples);
  canvas->Print("Overlay_Predictive.pdf]", "pdf");

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

  PlotMan = new M3::Plotting::PlottingManager();
  PlotMan->initialise();

  PredictivePlotting(ConfigName, FileNames);

  if(PlotMan) delete PlotMan;
  return 0;
}
