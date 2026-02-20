#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "PlottingUtils/PlottingUtils.h"
#include "PlottingUtils/PlottingManager.h"

_MaCh3_Safe_Include_Start_ //{
#include "TROOT.h"
#include "TGaxis.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TCandle.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"
_MaCh3_Safe_Include_End_ //}

/// @file GetPostfitParamPlots.cpp
/// This script generates post-fit parameter plots. The central postfit value is
/// taken as the Highest Posterior Density (HPD), but can be easily changed to
/// another method such as Gaussian. Be cautious as parameter names and the number
/// of parameters per plot are currently hardcoded.
///
/// @details
/// Usage:
/// ```
/// ./GetPostfitParamPlots ProcessMCMC_Output1.root <ProcessMCMC_Output2.root> <ProcessMCMC_Output3.root>
/// ```
///
/// @ingroup MaCh3Plotting
///
/// @author Clarence Wret
/// @author Will Parker
/// @author Kamil Skwarczynski
/// @author Ewan Miller

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

MaCh3Plotting::PlottingManager *PlotMan;
TH1D *Prefit;

int NDParameters;
int NDParametersStartingPos;

std::vector<int> NDSamplesBins;
std::vector<std::string> NDSamplesNames;
int nBins;
TCanvas *canv;

std::string SaveName;

TPad *p3;
TPad *p4;

// KS: Color for 0 - prefit, 1 postfit, 2 another postfit, 3 you know the drill
constexpr Color_t PlotColor[] = {kRed, kBlack, kBlue, kGreen};
std::string plotType;

void copyParToBlockHist(const int localBin, const std::string& paramName, TH1D* blockHist,
                        const std::string& type, const int fileId, const bool setLabels = true){
  // Set the values in the sub-histograms
  MACH3LOG_DEBUG("copying data from at local bin {}: for parameter {}", localBin, paramName);
  MACH3LOG_DEBUG("  Fitter specific name: {}", PlotMan->input().translateName(fileId, MaCh3Plotting::kPostFit, paramName));
  MACH3LOG_DEBUG("  value: {:.4f}", PlotMan->input().getPostFitValue(fileId, paramName, type));
  MACH3LOG_DEBUG("  error: {:.4f}", PlotMan->input().getPostFitError(fileId, paramName, type));

  blockHist->SetBinContent(localBin +1, PlotMan->input().getPostFitValue(fileId, paramName, type));
  blockHist->SetBinError(localBin +1, PlotMan->input().getPostFitError(fileId, paramName, type));

  if(setLabels){
    blockHist->GetXaxis()->SetBinLabel(localBin +1, paramName.c_str());
    blockHist->GetXaxis()->LabelsOption("v");
  }
}

inline void InitializePads(TCanvas* canvas, TPad*& pad3, TPad*& pad4) {
  // Initialize TPad p3
  pad3 = new TPad("Top", "Top", 0.0, 0.4, 1.0, 1.0);
  pad3->SetLeftMargin(canvas->GetLeftMargin());
  pad3->SetRightMargin(canvas->GetRightMargin());
  pad3->SetTopMargin(canvas->GetTopMargin());
  pad3->SetBottomMargin(0);
  pad3->SetGrid();

  // Initialize TPad p4
  pad4 = new TPad("Bottom", "Bottom", 0.0, 0.0, 1.0, 0.4);
  pad4->SetLeftMargin(canvas->GetLeftMargin());
  pad4->SetRightMargin(canvas->GetRightMargin());
  pad4->SetTopMargin(0);
  pad4->SetBottomMargin(0.75);
  pad4->SetGrid();
}

void CopyViolinToBlock(TH2D* FullViolin, TH2D* ReducedViolin, const std::vector<std::string>& ParamNames) {
  for(unsigned int i = 0; i < ParamNames.size(); i++)
  {
    int ParamBinId = M3::_BAD_INT_;
    for (int ix = 0; ix < FullViolin->GetXaxis()->GetNbins(); ++ix) {
      if(FullViolin->GetXaxis()->GetBinLabel(ix+1) == ParamNames[i])
      {
        ParamBinId = ix+1;
        break;
      }
    }
    if(ParamBinId == M3::_BAD_INT_) {
      MACH3LOG_WARN("Didn't find param {}", ParamNames[i]);
      continue;
    }
    //KS Fill content of reduced violin
    for (int iy = 0; iy < FullViolin->GetYaxis()->GetNbins(); ++iy) {
      ReducedViolin->SetBinContent(i+1, iy+1, FullViolin->GetBinContent(ParamBinId, iy+1));
      ReducedViolin->GetXaxis()->SetBinLabel(i+1, ParamNames[i].c_str());
    }
  }
  ReducedViolin->SetFillColor(FullViolin->GetFillColor());
  ReducedViolin->SetFillColorAlpha(FullViolin->GetMarkerColor(), 0.35);
  ReducedViolin->SetLineColor(FullViolin->GetMarkerColor());

  ReducedViolin->SetMarkerColor(FullViolin->GetMarkerColor());
  ReducedViolin->SetMarkerStyle(FullViolin->GetMarkerStyle());
  ReducedViolin->SetMarkerSize(FullViolin->GetMarkerSize());

  ReducedViolin->GetYaxis()->SetTitleOffset(FullViolin->GetTitleOffset());
  ReducedViolin->GetYaxis()->SetTitle(FullViolin->GetYaxis()->GetTitle());
  ReducedViolin->GetXaxis()->LabelsOption("v");
}

void PrettifyTitles(TH1D *Hist) {
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i)
  {
    std::string title = Hist->GetXaxis()->GetBinLabel(i+1);
    title = PlotMan->style().prettifyParamName(title);
    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
  }
}

void PrettifyTitles(TH2D *Hist) {
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i) 
  {
    std::string title = Hist->GetXaxis()->GetBinLabel(i+1);

    title = PlotMan->style().prettifyParamName(title);
    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
  }
}

bool ReadSettings(const std::shared_ptr<TFile>& File1)
{
  MACH3LOG_DEBUG("Reading settings for file {}", File1->GetName());
  #ifdef DEBUG
  File1->ls();
  #endif
  TTree *Settings = (File1->Get<TTree>("Settings"));

  // can't find settings tree :(
  if (!Settings) return false;

  MACH3LOG_DEBUG("Got settings tree");
  MaCh3Utils::Print(Settings);

  Settings->SetBranchAddress("NDParameters", &NDParameters);
  Settings->SetBranchAddress("NDParametersStartingPos", &NDParametersStartingPos);

  std::vector<int> *NDSamples_Bins = 0;
  std::vector<std::string> *NDSamples_Names = 0;
  Settings->SetBranchAddress("NDSamplesNames", &NDSamples_Names);
  Settings->SetBranchAddress("NDSamplesBins", &NDSamples_Bins);

  Settings->GetEntry(0);

  NDSamplesNames = *NDSamples_Names;
  NDSamplesBins = *NDSamples_Bins;

  MACH3LOG_DEBUG("Read settings tree successfully");
  return true;
}

std::unique_ptr<TH1D> makeRatio(TH1D *PrefitCopy, TH1D *PostfitCopy, bool setAxes) {
  // set up the ratio hist
  std::unique_ptr<TH1D> Ratio = M3::Clone(PrefitCopy);
  Ratio->GetYaxis()->SetTitle("(x_{Post}-#mu_{Prior})/#sigma_{Prior}");
  Ratio->SetMinimum(-3.7);
  Ratio->SetMaximum(3.7);    

  for (int j = 0; j < Ratio->GetXaxis()->GetNbins(); ++j) 
  {
    if ( PrefitCopy->GetBinError(j+1) > 1.e-5 )
    {
      Ratio->SetBinContent(j+1, (PostfitCopy->GetBinContent(j+1)-PrefitCopy->GetBinContent(j+1))/PrefitCopy->GetBinError(j+1));

      double up = (PostfitCopy->GetBinContent(j+1)+PostfitCopy->GetBinError(j+1)-PrefitCopy->GetBinContent(j+1))/PrefitCopy->GetBinError(j+1);
      double down = (PostfitCopy->GetBinContent(j+1)-PostfitCopy->GetBinError(j+1)-PrefitCopy->GetBinContent(j+1))/PrefitCopy->GetBinError(j+1);

      double maximum = up-Ratio->GetBinContent(j+1);
      double minimum = Ratio->GetBinContent(j+1)-down;

      Ratio->SetBinError(j+1, std::max(maximum, minimum));
    }
    //KS: Most likely flat prior
    else {
      Ratio->SetBinContent(j+1, (PostfitCopy->GetBinContent(j+1)-PrefitCopy->GetBinContent(j+1)));

      double up = (PostfitCopy->GetBinContent(j+1)+PostfitCopy->GetBinError(j+1)-PrefitCopy->GetBinContent(j+1));
      double down = (PostfitCopy->GetBinContent(j+1)-PostfitCopy->GetBinError(j+1)-PrefitCopy->GetBinContent(j+1));

      double maximum = up-Ratio->GetBinContent(j+1);
      double minimum = Ratio->GetBinContent(j+1)-down;

      Ratio->SetBinError(j+1, std::max(maximum, minimum));
    }
  } //end loop over parameters

  if(setAxes){
    Ratio->SetFillStyle(0);
    Ratio->SetFillColor(0);

    Ratio->SetLineColor(PostfitCopy->GetLineColor());
    if (Ratio->GetLineColor() == 0) Ratio->SetLineColor(kBlack);
    Ratio->SetMarkerColor(PostfitCopy->GetMarkerColor());

    Ratio->SetLineWidth(PlotMan->getOption<int>("plotLineWidth"));
    Ratio->SetTitle("");

    Ratio->SetMarkerSize(2);
    Ratio->SetMarkerStyle(20);

    Ratio->GetYaxis()->SetTitleSize(25);
    Ratio->GetYaxis()->SetTitleFont(43);
    Ratio->GetYaxis()->SetTitleOffset(2.0);
    Ratio->GetYaxis()->SetLabelFont(43);
    Ratio->GetYaxis()->SetLabelSize(25);
    Ratio->GetYaxis()->CenterTitle();
    Ratio->GetYaxis()->SetNdivisions(5,2,0);

    Ratio->GetXaxis()->SetTitleSize(25);
    Ratio->GetXaxis()->SetTitleFont(43);
    Ratio->GetXaxis()->SetTitleOffset(4.0);
    Ratio->GetXaxis()->SetLabelOffset(0.025);
    Ratio->GetXaxis()->SetLabelFont(43);
    Ratio->GetXaxis()->SetLabelSize(25);
  }

  return Ratio;
}

void DrawPlots(TCanvas *plotCanv, TH1D* PrefitCopy, const std::vector<std::unique_ptr<TH1D>>& PostfitVec, TPad *mainPad, TPad *ratioPad) {
  // Draw!
  plotCanv->cd();
  mainPad->Draw();
  mainPad->cd();
  PrefitCopy->GetYaxis()->SetTitle("Parameter Value");
  
  PrefitCopy->GetYaxis()->SetLabelSize(0.);
  PrefitCopy->GetYaxis()->SetTitleSize(0.05);
  PrefitCopy->GetYaxis()->SetTitleOffset(1.3);
  PrefitCopy->Draw("e2");

  for (int fileId = 0; fileId < static_cast<int>(PostfitVec.size()); fileId++) {
    TH1D *postFitHist = PostfitVec[fileId].get();
    
    postFitHist->SetMarkerColor(TColor::GetColorPalette(fileId));
    postFitHist->SetLineColor(TColor::GetColorPalette(fileId));
    postFitHist->SetMarkerStyle(7);
    postFitHist->SetLineStyle(1+fileId);
    postFitHist->SetLineWidth(PlotMan->getOption<int>("plotLineWidth"));

    postFitHist->Draw("e1, same");
  }
  
  plotCanv->Update();
  auto axis = std::make_unique<TGaxis>(PrefitCopy->GetXaxis()->GetBinLowEdge(PrefitCopy->GetXaxis()->GetFirst()), gPad->GetUymin()+0.01,
                            PrefitCopy->GetXaxis()->GetBinLowEdge(PrefitCopy->GetXaxis()->GetFirst()), gPad->GetUymax(),
                            gPad->GetUymin()+0.01, gPad->GetUymax(), 510, "");
  axis->SetLabelFont(43);
  axis->SetLabelSize(25);
  axis->Draw();

  plotCanv->cd();
  ratioPad->Draw();
  ratioPad->cd();
  
  std::vector<std::unique_ptr<TH1D>> ratioHists;
  // save pointers to these so we can delete them once we are done
  ratioHists.push_back(makeRatio(PrefitCopy, PostfitVec[0].get(), true));

  ratioHists[0]->Draw("p");
  for(int postFitIdx = 1; postFitIdx < static_cast<int>(PostfitVec.size()); postFitIdx++){
    ratioHists.push_back(makeRatio(PrefitCopy, PostfitVec[postFitIdx].get(), true));
    
    ratioHists[postFitIdx]->SetMarkerColor(TColor::GetColorPalette(postFitIdx));
    ratioHists[postFitIdx]->SetLineColor(TColor::GetColorPalette(postFitIdx));
    ratioHists[postFitIdx]->SetMarkerStyle(7);
    ratioHists[postFitIdx]->SetLineStyle(1+postFitIdx);
    ratioHists[postFitIdx]->SetLineWidth(PlotMan->getOption<int>("plotLineWidth"));

    ratioHists[postFitIdx]->Draw("p same");
  }

  // draw lines across the plot at +-1 and 0
  TLine line(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), 0.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), 0.0);
  TLine line2(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), 1.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), 1.0);
  TLine line3(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), -1.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), -1.0);
  
  line.SetLineColor(kRed);
  line.SetLineStyle(kDashed);
  line.SetLineWidth(PlotMan->getOption<int>("refLineWidth"));
  line2.SetLineColor(kRed);
  line2.SetLineStyle(kDashed);
  line2.SetLineWidth(PlotMan->getOption<int>("refLineWidth"));
  line3.SetLineColor(kRed);
  line3.SetLineStyle(kDashed);
  line3.SetLineWidth(PlotMan->getOption<int>("refLineWidth"));

  line.Draw("same");
  line2.Draw("same");
  line3.Draw("same");

  plotCanv->Print((SaveName).c_str());
}

void MakeParameterPlots()
{  
  // get the names of the blocks of parameters to group together
  std::vector<std::string> const blockNames = PlotMan->getOption<std::vector<std::string>>("paramGroups");
  const int nPlots = static_cast<int>(blockNames.size());

  for (int i = 0; i < nPlots; i++)
  {
    // get the configuration for this parameter
    std::string blockName = blockNames[i];
    YAML::Node paramBlock = PlotMan->getOption(blockName);
    auto blockTitle       = paramBlock[0].as<std::string>();
    auto blockLimits      = paramBlock[1].as<std::vector<double>>();
    auto blockContents    = paramBlock[2].as<std::vector<std::string>>();

    // get num of params in the block
    const int nParams = static_cast<int>(blockContents.size());

    // set some plot things
    auto blockHist_prefit = std::make_unique<TH1D>(blockName.c_str(), blockTitle.c_str(),
                                                  nParams, 0.0, static_cast<double>(nParams));
    
    PlotMan->style().setTH1Style(blockHist_prefit.get(), PlotMan->getOption<std::string>("prefitHistStyle"));

    // set the errors for the prefit block hist
    for(int localBin=0; localBin < nParams; localBin ++){
      // the "local" bin is the params index within the group of parameters
      std::string paramName = blockContents[localBin];
      copyParToBlockHist(localBin, paramName, blockHist_prefit.get(), "Prior", 0);
    }

    // now set for the postfit blocks for all files
    std::vector <std::unique_ptr<TH1D>> blockHist_postfit_Vec;
    for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++) {
      auto blockHist_postfit = std::make_unique<TH1D>((blockName + PlotMan->getFileName(fileId)).c_str(),
                                                      blockTitle.c_str(), nParams, 0.0, static_cast<double>(nParams));

      // loop through all the parameters in this block and set the contents in the blocks TH1
      for(int localBin=0; localBin < nParams; localBin ++){
        // the "local" bin is the params index within the group of parameters
        std::string paramName = blockContents[localBin];
        copyParToBlockHist(localBin, paramName, blockHist_postfit.get(), "", fileId);
      }
      blockHist_postfit_Vec.push_back(std::move(blockHist_postfit));
    } 
    // set the y axis limits we got from config
    blockHist_prefit->GetYaxis()->SetRangeUser(blockLimits[0], blockLimits[1]); 

    // Do some fancy replacements
    PrettifyTitles(blockHist_prefit.get());

    DrawPlots(canv, blockHist_prefit.get(), blockHist_postfit_Vec, p3, p4);
  }
}

void MakeFluxPlots()
{
  // these for non named params where we don't need as much space
  auto p1 = std::make_unique<TPad>("p1", "p1", 0.0, 0.3, 1.0, 1.0);
  auto p2 = std::make_unique<TPad>("p2", "p2", 0.0, 0.0, 1.0, 0.3);
  p1->SetLeftMargin(canv->GetLeftMargin());
  p1->SetRightMargin(canv->GetRightMargin());
  p1->SetTopMargin(canv->GetTopMargin());
  p1->SetBottomMargin(0);
  p1->SetGrid();
  p2->SetLeftMargin(canv->GetLeftMargin());
  p2->SetRightMargin(canv->GetRightMargin());
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.25);
  p2->SetGrid();

  p1->SetLogx(true);
  p2->SetLogx(true);
  
  // get the names of the blocks of parameters to group together
  std::vector<std::string> const fluxBlockNames = PlotMan->getOption<std::vector<std::string>>("fluxGroups");
  auto const fluxBinningTable = PlotMan->getOption("FluxBinning");

  const int FluxPlots = static_cast<int>(fluxBlockNames.size());

  for (int i = 0; i < FluxPlots; i++) 
  {
    // get the configuration for this block
    std::string fluxBlockName = fluxBlockNames[i];
    YAML::Node paramBlock           = PlotMan->getOption(fluxBlockName);
    std::string blockTitle          = paramBlock[0].as<std::string>();
    std::vector<double> blockLimits = paramBlock[1].as<std::vector<double>>();
    std::string blockBinningName    = paramBlock[2].as<std::string>();
    std::vector<int> blockContents  = paramBlock[3].as<std::vector<int>>();

    // get the binning for this block of flux params
    std::vector<double> binning    = fluxBinningTable[blockBinningName].as<std::vector<double>>();

    // get num of params in the block
    int nParams = blockContents[1] - blockContents[0] +1;
    // check for sanity
    if(nParams <= 0 || blockContents.size() > 2){
      MACH3LOG_CRITICAL("Invalid flux parameter block endpoints specified for {}", fluxBlockName);
      MACH3LOG_CRITICAL("  Should have the form [<low index>, <up index>]");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    if (nParams != static_cast<int>(binning.size()) - 1) {
      MACH3LOG_CRITICAL("Binning provided for flux param block {} does not match the number of parameters specified for the block", fluxBlockName);
      MACH3LOG_CRITICAL("  Provided {} parameters but {} bins", nParams, binning.size() -1);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    auto blockHist_prefit = std::make_unique<TH1D>(fluxBlockName.c_str(), blockTitle.c_str(), nParams, binning.data());
    blockHist_prefit->GetYaxis()->SetTitle("Parameter Variation");
    blockHist_prefit->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    blockHist_prefit->GetXaxis()->SetTitleOffset(blockHist_prefit->GetXaxis()->GetTitleOffset()*1.2);
    PlotMan->style().setTH1Style(blockHist_prefit.get(), PlotMan->getOption<std::string>("prefitHistStyle"));
    // set the errors for the prefit block hist
    for(int fluxParId = blockContents[0]; fluxParId <= blockContents[1]; fluxParId++){
        int localBin = fluxParId - blockContents[0];
        std::string paramName = "b_" + std::to_string(fluxParId);
        copyParToBlockHist(localBin, paramName, blockHist_prefit.get(), "Prior", 0, false);
    }

    // now set for the postfit blocks for all files
    std::vector <std::unique_ptr<TH1D>> blockHist_postfit_Vec;
    for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++){
      auto blockHist_postfit = std::make_unique<TH1D>(fluxBlockName.c_str(), blockTitle.c_str(), nParams, binning.data());

      for(int fluxParId = blockContents[0]; fluxParId <= blockContents[1]; fluxParId++){
        int localBin = fluxParId - blockContents[0];
        std::string paramName = "b_" + std::to_string(fluxParId);

        copyParToBlockHist(localBin, paramName, blockHist_postfit.get(), "", fileId, false);
      }
      blockHist_postfit_Vec.push_back(std::move(blockHist_postfit));
    }
    // set the y axis limits we got from config
    blockHist_prefit->GetYaxis()->SetRangeUser(blockLimits[0], blockLimits[1]); 
    DrawPlots(canv, blockHist_prefit.get(), blockHist_postfit_Vec, p1.get(), p2.get());
  }

  canv->cd();
  canv->SetLogx(false);
  canv->SetBottomMargin(canv->GetBottomMargin()*1.7);
}

/// @warning This is legacy functions and will become deprecated
void MakeNDDetPlots()
{
  MACH3LOG_INFO("ND detector parameters: {}", NDParameters);
  Prefit->GetYaxis()->SetTitleOffset(Prefit->GetYaxis()->GetTitleOffset()*1.2);

  TPad* pTop = nullptr;
  TPad* pDown = nullptr;
  InitializePads(canv, pTop, pDown);

  int NDbinCounter = NDParametersStartingPos;
  int Start = NDbinCounter;
  
  MACH3LOG_INFO("Running on {} samples", NDSamplesNames.size());

  for (unsigned int i = 0; i < NDSamplesNames.size(); ++i)
  {
    MACH3LOG_DEBUG("--- On sample {}", NDSamplesNames[i]); 
    NDbinCounter += NDSamplesBins[i];

    std::vector<std::unique_ptr<TH1D>> PostfitNDDetHistVec(PlotMan->getNFiles());
    TH1D *PreFitNDDetHist = PlotMan->input().getFile(0).file->Get<TH1D>(Form("param_%s_prefit", NDSamplesNames[i].c_str()));
    PlotMan->style().setTH1Style(PreFitNDDetHist, PlotMan->getOption<std::string>("prefitHistStyle"));

    std::string temp = NDSamplesNames[i].c_str();
    while (temp.find("_") != std::string::npos) {
       temp.replace(temp.find("_"), 1, std::string(" "));
     }
    PreFitNDDetHist->SetTitle(temp.c_str());
    PreFitNDDetHist->GetXaxis()->SetRangeUser(Start, NDbinCounter);
    
    MACH3LOG_DEBUG("  Start bin: {} :: End bin: {}", Start, NDbinCounter);
    // set the x range for the postfits
    for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++){
      PostfitNDDetHistVec[fileId] = M3::Clone(PlotMan->input().getFile(fileId).file->Get<TH1D>(Form("param_%s_%s", NDSamplesNames[i].c_str(), plotType.c_str())));
    }

    //KS: We dont' need name for every nd param
    for(int j = 0; j < NDSamplesBins[i]; ++j)
    {
      bool ProductOfTen = false;
      if(j % 10) ProductOfTen = true;
      if(j != 0 && ProductOfTen) PreFitNDDetHist->GetXaxis()->SetBinLabel(Start+j+1, " ");
      else{
        PreFitNDDetHist->GetXaxis()->SetBinLabel(Start+j+1, Form("Det Variation Bin %i", Start+j));
      }
    }

    PreFitNDDetHist->GetYaxis()->SetRangeUser(PlotMan->getOption<double>("detParYRange_low"), PlotMan->getOption<double>("detParYRange_high"));
    
    Start += NDSamplesBins[i];

    DrawPlots(canv, PreFitNDDetHist, PostfitNDDetHistVec, pTop, pDown);
    canv->Update();
  }
  delete pTop;
  delete pDown;
}

void MakeRidgePlots()
{  
  gStyle->SetPalette(51);

  auto blankCanv = std::make_unique<TCanvas>("blankCanv", "blankCanv", 2048, 2048);
  blankCanv->SaveAs("RidgePlots.pdf[");

  // get the names of the blocks of parameters to group together
  const auto blockNames = PlotMan->getOption<std::vector<std::string>>("paramGroups");
  const int nPlots = static_cast<int>(blockNames.size());

  constexpr double padTopMargin = 0.9;
  constexpr double padBottomMargin = 0.1;
  constexpr double padOverlap = 0.9;
  constexpr double ridgeLineWidth = 1.0;
  for (int i = 0; i < nPlots; i++)
  {
    // get the configuration for this parameter
    std::string blockName = blockNames[i];
    auto const &paramBlock = PlotMan->getOption(blockName);
    auto blockTitle    = paramBlock[0].as<std::string>();
    auto blockLimits   = paramBlock[1].as<std::vector<double>>();
    auto blockContents = paramBlock[2].as<std::vector<std::string>>();

    // the directory of histograms
    TDirectoryFile *posteriorDir = PlotMan->input().getFile(0).file->Get<TDirectoryFile>("Post_1d_hists");

    // get num of params in the block
    int nParams = static_cast<int>(blockContents.size());

    if (nParams == 1) {
      MACH3LOG_WARN("{} doesn't work for single param", __func__);
      continue;
    }
    auto ridgeCanv = std::make_unique<TCanvas>("RidgePlotCanv", "RidgePlotCanv", 2048, 2048);
    ridgeCanv->Divide(1,1+nParams, 0.01, 0.0);

    auto title = std::make_unique<TLatex>();
    title->SetTextAlign(21);
    title->SetTextSize(0.03);
    title->DrawLatex(0.5, 0.95, blockTitle.c_str());

    auto label = std::make_unique<TLatex>();
    label->SetTextAlign(31);
    label->SetTextSize(0.02);
    
    auto line = std::make_unique<TLine>();
    line->SetLineColor(kBlack);
    line->SetLineWidth(ridgeLineWidth);

    // use this to set the limits and also to plot the x axis and grid
    auto axisPlot = std::make_unique<TH1D>("axis plot", "", 1, blockLimits[0], blockLimits[1]);

    std::vector<std::unique_ptr<TH1D>> axisPlot_holder(nParams);
    std::vector<std::unique_ptr<TPad>> graph_holder(nParams);
    for(int parId = 0; parId < nParams; parId++) {
      std::string paramName = blockContents[parId];

      TH1D *posteriorDist = nullptr;
      // get the list of objects in the directory
      TIter next(posteriorDir->GetListOfKeys());
      while (TKey* key = static_cast<TKey*>(next())) {
        // check if the end of the param name matches with the MaCh3 name, do this so we exclude things like nds_ at the start of the name
        std::string str(key->GetTitle());
        std::string name = PlotMan->input().translateName(0, MaCh3Plotting::kPostFit, paramName);
        uint pos = str.find(name);
        bool foundPar = (pos == str.length() - name.length());

        MACH3LOG_TRACE("Looking for {} in {}", name, str);
        if(foundPar){
          MACH3LOG_TRACE("Found it");
          posteriorDist = posteriorDir->Get<TH1D>(key->GetName());
        }
      }

      if(posteriorDist == nullptr){
        MACH3LOG_WARN("Couldn't find parameter {} when making ridgeline plots", paramName);
        MACH3LOG_WARN("It could be fixed param");
        continue;
      }

      // EM: do some funky scaling so that we always get evenly spaced pads in the range [bottomMargin, TopMargin] with the specified overlap
      double padAnchor = padBottomMargin + (static_cast<double>(nParams - parId - 1) /
                                  static_cast<double>(nParams - 1)) * (padTopMargin - padBottomMargin);
      double padWidth = (padTopMargin - padBottomMargin) / static_cast<double>(nParams);
      double norm = (padTopMargin - padBottomMargin);

      double padTop = padWidth * (1.0 + padOverlap) * (padTopMargin - padAnchor) / norm + padAnchor;
      double padBottom = padAnchor - padWidth * (1.0 + padOverlap) * (padAnchor - padBottomMargin) / norm;

      auto pad = std::make_unique<TPad>(paramName.c_str(), "", 0.3, padBottom, 0.9, padTop, -1, 0, -1);
      ridgeCanv->cd();

      pad->SetBottomMargin(0.0);
      pad->SetTopMargin(0.0);
      pad->SetLeftMargin(0.0);
      pad->SetRightMargin(0.0);

      pad->Draw();
      pad->cd();
      pad->SetFillStyle(4000);

      gPad->SetFrameFillStyle(4000);
      posteriorDist->GetFunction("Gauss")->SetBit(TF1::kNotDraw);
      posteriorDist->SetTitle("");
      posteriorDist->SetLineWidth(ridgeLineWidth);
      
      auto axisPlot_tmp = M3::Clone(axisPlot.get(), Form("AxisPlot_%s", paramName.c_str()));
      axisPlot_tmp->Draw("A");
      posteriorDist->Draw("H SAME");
    
      axisPlot_tmp->GetYaxis()->SetRangeUser(0.0, 0.7 *posteriorDist->GetMaximum());
      posteriorDist->SetLineColor(kWhite);
      posteriorDist->SetFillColorAlpha(TColor::GetColorPalette(floor(static_cast<float>(parId) *
                                          TColor::GetNumberOfColors() / static_cast<float>(nParams))), 0.85);
      
      posteriorDist->GetXaxis()->SetRangeUser(blockLimits[0], blockLimits[1]);
      posteriorDist->GetYaxis()->SetTitle(paramName.c_str());
      
      //EM: Have to get the frame drawn by the histogram and then set it to be transparent... it took me an hour to get rid of this one line on a plot
      gPad->Modified(); gPad->Update();
      TFrame *frame = gPad->GetFrame();
      frame->SetLineColorAlpha(0, 0.0);
      
      ridgeCanv->cd();
      label->DrawLatexNDC(0.29, padBottom + 0.005, PlotMan->style().prettifyParamName(paramName).c_str());
      line->DrawLine(0.1, padBottom, 0.9, padBottom);

      axisPlot_holder[parId] = std::move(axisPlot_tmp);
      graph_holder[parId] = std::move(pad);
    }

    ridgeCanv->cd();
    ridgeCanv->SetGrid(1,1);
    auto axisPad = std::make_unique<TPad>("AxisPad", "", 0.3, 0.0, 0.9, 1.0, -1, 0, -1);
    axisPad->SetLeftMargin(0.0);
    axisPad->SetRightMargin(0.0);
    axisPad->Draw();
    axisPad->cd();
    axisPad->SetGrid(1,1);
    axisPad->SetFrameFillStyle(4000);
    
    axisPlot->GetXaxis()->SetTickSize(0.01);
    axisPlot->GetXaxis()->SetTitle("Parameter Variation");
    axisPlot->GetYaxis()->SetLabelOffset(9999);
    axisPlot->GetYaxis()->SetLabelSize(0);
    axisPlot->GetYaxis()->SetTickSize(0);
    axisPlot->GetYaxis()->SetAxisColor(0,0.0);
    axisPlot->Draw("AXIS");
    axisPlot->Draw("AXIG SAME");

    axisPlot->SetFillStyle(4000);
    axisPad->SetFillStyle(4000);
    
    axisPad->SetGrid(1,1);
    gPad->Modified(); gPad->Update();
    gPad->SetFrameFillStyle(4000);

    gPad->Modified(); gPad->Update();
    TFrame *frame = gPad->GetFrame();
    frame->SetLineColorAlpha(0, 0.0);

    ridgeCanv->SaveAs("RidgePlots.pdf");
  }
  blankCanv->SaveAs("RidgePlots.pdf]");
}

void GetPostfitParamPlots()
{
  SaveName = PlotMan->getOutputName();
  
  //KS: By default we take HPD values, but by changing "plotType" you can use for example Gauss
  plotType = "HPD";
  //plotType = "gaus"; 
    
  MACH3LOG_INFO("Plotting {} errors", plotType);
  
  // if we have one MaCh3 nd file then we can get settings from it 
  bool plotNDDet = false;
  for (size_t fileId = 0; fileId < PlotMan->input().getNInputFiles(); fileId++) {
    if(!ReadSettings(PlotMan->input().getFile(0).file)) {
      MACH3LOG_INFO("at least one file provided does not have 'settings' tree indicating it is not MaCh3 ND file");
      MACH3LOG_INFO("  sadly this means I cannot plot ND Det parameters as this is only supported for MaCh3 ND files for now... sorry :(");
      plotNDDet = false;
    }
  }

  canv = new TCanvas("canv", "canv", 1024, 1024);
  //gStyle->SetPalette(51);
  gStyle->SetOptStat(0); //Set 0 to disable statistic box
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.12);
  canv->SetTopMargin(0.08);
  canv->SetRightMargin(0.04);

  canv->Print((SaveName+"[").c_str());

  // these for named parameters where we need a nice big gap at the botto to fit the names
  InitializePads(canv, p3, p4);

  // Make a Legend page
  auto leg = std::make_unique<TLegend>(0.0, 0.0, 1.0, 1.0);
  // make a dummy TH1 to set out legend
  Prefit = new TH1D();
  PlotMan->style().setTH1Style(Prefit, PlotMan->getOption<std::string>("prefitHistStyle"));
  leg->AddEntry(Prefit, "Prior", "lpf");

  for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++){
    TH1D *postFitHist_tmp = new TH1D();
    postFitHist_tmp->SetBit(kCanDelete);
    
    postFitHist_tmp->SetMarkerColor(TColor::GetColorPalette(fileId));
    postFitHist_tmp->SetLineColor(TColor::GetColorPalette(fileId));
    postFitHist_tmp->SetMarkerStyle(7);
    postFitHist_tmp->SetLineStyle(1+fileId);
    postFitHist_tmp->SetLineWidth(PlotMan->getOption<int>("plotLineWidth"));
    leg->AddEntry(postFitHist_tmp, PlotMan->getFileLabel(fileId).c_str(), "lpf");
  }

  canv->cd();
  canv->Clear();
  leg->Draw();
  canv->Print((SaveName).c_str());

  MakeParameterPlots();

  MakeFluxPlots();

  //KS: By default we don't run ProcessMCMC with PlotDet as this take some time, in case we did let's make fancy plots
  if(plotNDDet & (NDParameters > 0)) MakeNDDetPlots();

  canv->Print((SaveName+"]").c_str());

  MakeRidgePlots();

  delete canv;
  delete Prefit;
}

std::unique_ptr<TGraphAsymmErrors> MakeTGraphAsymmErrors(const std::shared_ptr<TFile>& File,  std::vector<int> Index = {})
{
  int GraphBins = Index.size() == 0 ? nBins : Index.size();
  std::vector<double> x(GraphBins);
  std::vector<double> y(GraphBins);
  std::vector<double> exl(GraphBins);
  std::vector<double> eyl(GraphBins);
  std::vector<double> exh(GraphBins);
  std::vector<double> eyh(GraphBins);

  TH1D* PostHist = static_cast<TH1D*>(File->Get( ("param_xsec_"+plotType).c_str() ));

  auto Errors_HPD_Positive = static_cast<TVectorD*>(File->Get( "Errors_HPD_Positive" ));
  auto Errors_HPD_Negative = static_cast<TVectorD*>(File->Get( "Errors_HPD_Negative" ));
  //KS: I am tempted to multithread this...
  for(int i = 0; i < GraphBins; ++i)
  {
    int Counter = Index.size() == 0 ? i : Index[i];
    //KS: We are extracting value from three object each having different numbering scheme, I have checked carefully so this is correct please don't change all these +1 +0.5 etc. it just work...
    x[i] = i + 0.5;
    y[i] = PostHist->GetBinContent(Counter+1);

    //KS: We don't want x axis errors as they are confusing in Violin plot
    exh[i] = 0.00001;
    exl[i] = 0.00001;
    eyh[i] = (*Errors_HPD_Positive)(Counter);
    eyl[i] = (*Errors_HPD_Negative)(Counter);
  }
  auto PostGraph = std::make_unique<TGraphAsymmErrors>(GraphBins, x.data(), y.data(), exl.data(), exh.data(), eyl.data(), eyh.data());
  PostGraph->SetTitle("");

  return PostGraph;
}
   
/// @brief KS: Make fancy violin plots
void GetViolinPlots()
{
  //KS: Should be in some config... either way it control whether you plot symmetric or asymmetric error bars
  bool PlotAssym = true;
    
  //KS: No idea why but ROOT changed treatment of violin in R6. If you have non uniform binning this will results in very hard to see violin plots.
  TCandle::SetScaledViolin(false);

  std::string OutputName = "";
  for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++){
    MACH3LOG_INFO("File {}: {} ", fileId, PlotMan->getFileName(fileId));
    OutputName += PlotMan->getFileName(fileId);
    OutputName = OutputName.substr(0, OutputName.find(".root"));
  }
  MACH3LOG_INFO("Making Violin Plot");
  OutputName += "_Violin";
  if(PlotAssym) OutputName += "_Assym";

  auto canvas = std::make_unique<TCanvas>("canv", "canv", 1024, 1024);
  canvas->SetGrid();
  gStyle->SetOptStat(0);
  //KS: Remove errors on X axis as they are confusing in violin type of plot
  if(!PlotAssym) gStyle->SetErrorX(0.0001);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetBottomMargin(0.25);
  canvas->SetTopMargin(0.08);
  canvas->SetRightMargin(0.03);
  canvas->SetLeftMargin(0.10);
  canvas->Print((OutputName+".pdf[").c_str());
  canvas->SetGrid();

  if(PlotMan->input().getFile(0).file->Get<TH2D>( "param_violin_prior" ) == nullptr)
  {
    MACH3LOG_WARN("Couldn't find violin plot, make sure method from MCMCProcessor is being called");
    return;
  }
  std::unique_ptr<TH2D> ViolinPre = M3::Clone(PlotMan->input().getFile(0).file->Get<TH2D>( "param_violin_prior" ));
  // Do some fancy replacements
  ViolinPre->SetFillColor(kRed);
  ViolinPre->SetFillColorAlpha(kRed, 0.35);
  ViolinPre->SetMarkerColor(kRed);
  ViolinPre->SetMarkerStyle(20);
  ViolinPre->SetMarkerSize(0.5);

  ViolinPre->GetYaxis()->SetTitleOffset(1.3);
  ViolinPre->GetYaxis()->SetTitle("Parameter Value");
  ViolinPre->GetXaxis()->LabelsOption("v");

  std::unique_ptr<TH1D> Postfit = M3::Clone(PlotMan->input().getFile(0).file->Get<TH1D>( ("param_xsec_"+plotType).c_str() ));
  Postfit->SetMarkerColor(kRed);
  Postfit->SetLineColor(kRed);
  Postfit->SetMarkerStyle(7);

  std::vector<std::unique_ptr<TH2D>> Violin(PlotMan->getNFiles());
  for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++) {
    Violin[fileId] = M3::Clone(PlotMan->input().getFile(fileId).file->Get<TH2D>( "param_violin" ));
    if(Violin[fileId] == nullptr)
    {
      MACH3LOG_ERROR("Couldn't find violin plot, make sure method from MCMCProcessor is being called");
      return;
    }
    //KS: I know hardcoded but we can figure out later...
    if(fileId == 0){
      Violin[fileId]->SetFillColor(kBlue);
      Violin[fileId]->SetFillColorAlpha(kBlue, 0.35);
      Violin[fileId]->SetMarkerColor(kBlue);
      Violin[fileId]->SetMarkerStyle(20);
      Violin[fileId]->SetMarkerSize(0.5);
    } else if (fileId == 1) {
      Violin[fileId]->SetMarkerColor(kGreen);
      Violin[fileId]->SetLineColor(kGreen);
      Violin[fileId]->SetFillColor(kGreen);
      Violin[fileId]->SetFillColorAlpha(kGreen, 0.35);
    } else if (fileId == 2) {
      Violin[fileId]->SetMarkerColor(kMagenta);
      Violin[fileId]->SetLineColor(kMagenta);
      Violin[fileId]->SetFillColor(kMagenta);
      Violin[fileId]->SetFillColorAlpha(kMagenta, 0.35);
    } else {
      MACH3LOG_ERROR("Too many file, not implemented...");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }

  std::unique_ptr<TGraphAsymmErrors> PostGraphAll = MakeTGraphAsymmErrors(PlotMan->input().getFile(0).file);
  // Make a Legend page
  auto leg = std::make_unique<TLegend>(0.0, 0.0, 1.0, 1.0);
  if (ViolinPre != nullptr) leg->AddEntry(ViolinPre.get(), "Prior", "lpf");
  for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++) {
    leg->AddEntry(Violin[fileId].get(), PlotMan->getFileLabel(fileId).c_str(), "lpf");
  }
  if(PlotAssym) leg->AddEntry(PostGraphAll.get(), "HPD Assym", "lp");
  else          leg->AddEntry(Postfit.get(), "HPD", "lpf");

  canvas->cd();
  canvas->Clear();
  leg->Draw();
  canvas->Print((OutputName+".pdf").c_str());
  
  // get the names of the blocks of parameters to group together
  std::vector<std::string> const blockNames = PlotMan->getOption<std::vector<std::string>>("paramGroups");
  const int nPlots = static_cast<int>(blockNames.size());

  for (int i = 0; i < nPlots; i++)
  {
    // get the configuration for this parameter
    std::string blockName = blockNames[i];
    YAML::Node paramBlock = PlotMan->getOption(blockName);
    std::string blockTitle                 = paramBlock[0].as<std::string>();
    std::vector<double>      blockLimits   = paramBlock[1].as<std::vector<double>>();
    std::vector<std::string> blockContents = paramBlock[2].as<std::vector<std::string>>();

    // get num of params in the block
    const int nParams = static_cast<int>(blockContents.size());

    // set some plot things
    auto blockHist_prefit = std::make_unique<TH2D>((blockName + "_Prefit").c_str(), blockTitle.c_str(), nParams, 0.0, static_cast<double>(nParams),
                                      ViolinPre->GetYaxis()->GetNbins(), ViolinPre->GetYaxis()->GetXmin(), ViolinPre->GetYaxis()->GetXmax());
    CopyViolinToBlock(ViolinPre.get(), blockHist_prefit.get(), blockContents);
    // set the y axis limits we got from config
    blockHist_prefit->GetYaxis()->SetRangeUser(blockLimits[0], blockLimits[1]);

    std::vector<std::unique_ptr<TH2D>> blockHist(PlotMan->getNFiles());
    for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++) {
      blockHist[fileId] = std::make_unique<TH2D>((blockTitle + "Violin" + fileId).Data(), (blockTitle + "Violin" + fileId).Data(),
                                                   nParams, 0.0, static_cast<double>(nParams), Violin[fileId]->GetYaxis()->GetNbins(),
                                                   Violin[fileId]->GetYaxis()->GetXmin(), Violin[fileId]->GetYaxis()->GetXmax());
      CopyViolinToBlock(Violin[fileId].get(), blockHist[fileId].get(), blockContents);
      blockHist[fileId]->GetYaxis()->SetRangeUser(blockLimits[fileId], blockLimits[1]);
    }
    // Do some fancy replacements
    PrettifyTitles(blockHist_prefit.get());

    // set some plot things
    auto blockHist_Best = std::make_unique<TH1D>(blockName.c_str(), blockTitle.c_str(),
                                                 nParams, 0.0, static_cast<double>(nParams));
    // set the errors for the prefit block hist
    for(int localBin=0; localBin < nParams; localBin ++){
      // the "local" bin is the params index within the group of parameters
      std::string paramName = blockContents[localBin];
      copyParToBlockHist(localBin, paramName, blockHist_Best.get(), "", 0);
    }

    std::vector<int> Index;
    for(unsigned int is = 0; is < blockContents.size(); is++) {
      int ParamBinId = M3::_BAD_INT_;
      for (int ix = 0; ix < ViolinPre->GetXaxis()->GetNbins(); ++ix) {
        if(ViolinPre->GetXaxis()->GetBinLabel(ix+1) == blockContents[is]) {
          ParamBinId = ix;
          break;
        }
      }
      Index.push_back(ParamBinId);
    }
    std::unique_ptr<TGraphAsymmErrors> PostGraph = MakeTGraphAsymmErrors(PlotMan->input().getFile(0).file, Index);
    PostGraph->SetMarkerColor(kBlack);
    PostGraph->SetLineColor(kBlack);
    PostGraph->SetMarkerStyle(7);
    PostGraph->SetLineWidth(2);
    PostGraph->SetLineStyle(kSolid);

    blockHist_prefit->Draw("violinX(03100300)");
    for(unsigned int fileId = 0; fileId < PlotMan->getNFiles(); fileId++) {
      blockHist[fileId]->Draw("violinX(03100300) SAME");
    }

    if(PlotAssym) PostGraph->Draw("P SAME");
    else Postfit->Draw("SAME");
    canvas->Print((OutputName+".pdf").c_str());
  }
  canvas->Print((OutputName+".pdf]").c_str());
}

/// @brief KS: Make comparison of 2D Posteriors
void Get2DComparison(const std::string& FileName1, const std::string& FileName2)
{
  auto canvas = std::make_unique<TCanvas>("canvas", "canvas", 0, 0, 1024, 1024);
  canvas->SetBottomMargin(0.1f);
  canvas->SetTopMargin(0.05f);
  canvas->SetRightMargin(0.03f);
  canvas->SetLeftMargin(0.15f);

  // Open the two ROOT files
  TFile* File1 = M3::Open(FileName1, "READ", __FILE__, __LINE__);
  TFile* File2 = M3::Open(FileName2, "READ", __FILE__, __LINE__);

  // Get the Post_2d_hists directory from both files
  TDirectory* Dir1 = File1->Get<TDirectory>("Post_2d_hists");
  TDirectory* Dir2 = File2->Get<TDirectory>("Post_2d_hists");

  if (!Dir1 || !Dir2) {
    MACH3LOG_WARN("Post_2d_hists directory not found in one or both files while running {}.", __func__);
    File1->Close();
    delete File1;
    File2->Close();
    delete File2;
    return;
  }

  // Get the list of keys in the first directory
  TIter next1(Dir1->GetListOfKeys());
  TKey* key1 = nullptr;

  // Prepare the output PDF filename
  std::string SaveName2D = "2DComparison_" + FileName1 + "_" + FileName2;
  SaveName2D = SaveName2D.substr(0, SaveName2D.find(".root"));
  SaveName2D = SaveName2D + ".pdf";

  canvas->Print((SaveName2D+"[").c_str());
  // Loop over keys in the first directory
  while ((key1 = static_cast<TKey*>(next1()))) {
    TString histName = key1->GetName();

    // Check if the key is a TH2D
    if (TString(key1->GetClassName()) == "TH2D") {
      TH2D* hist1 = static_cast<TH2D*>(key1->ReadObj());

      // Try to get the histogram with the same name from the second directory
      TH2D* hist2 = static_cast<TH2D*>(Dir2->Get(histName));

      if (hist2) {
        hist1->SetTitle("");
        hist1->SetTitle("");

        // Prettify axis titles
        std::string Xtitle = PlotMan->style().prettifyParamName(hist1->GetXaxis()->GetTitle());
        std::string Ytitle = PlotMan->style().prettifyParamName(hist1->GetYaxis()->GetTitle());

        // Adjust the axis ranges of hist1 to include both histograms
        double xmin = std::min(hist1->GetXaxis()->GetXmin(), hist2->GetXaxis()->GetXmin());
        double xmax = std::max(hist1->GetXaxis()->GetXmax(), hist2->GetXaxis()->GetXmax());
        double ymin = std::min(hist1->GetYaxis()->GetXmin(), hist2->GetYaxis()->GetXmin());
        double ymax = std::max(hist1->GetYaxis()->GetXmax(), hist2->GetYaxis()->GetXmax());

        hist1->GetXaxis()->SetRangeUser(xmin, xmax);
        hist1->GetYaxis()->SetRangeUser(ymin, ymax);

        hist1->GetXaxis()->SetTitle(Xtitle.c_str());
        hist1->GetYaxis()->SetTitle(Ytitle.c_str());

        hist1->SetLineColor(kBlue);
        hist1->SetLineStyle(kSolid);
        hist1->SetLineWidth(2);

        hist2->SetLineColor(kRed);
        hist2->SetLineStyle(kDashed);
        hist2->SetLineWidth(2);

        hist1->Draw("CONT3");
        hist2->Draw("CONT3 SAME");

        auto Legend = std::make_unique<TLegend>(0.20, 0.7, 0.4, 0.92);
        Legend->AddEntry(hist1, PlotMan->getFileLabel(0).c_str(), "l");
        Legend->AddEntry(hist2, PlotMan->getFileLabel(1).c_str(), "l");
        Legend->SetTextSize(0.03);
        Legend->SetLineColor(0);
        Legend->SetLineStyle(0);
        Legend->SetFillColor(0);
        Legend->SetFillStyle(0);
        Legend->SetBorderSize(0);
        Legend->Draw("SAME");
        canvas->Print((SaveName2D).c_str());
      }
    }
  }
  canvas->Print((SaveName2D+"]").c_str());

  File1->Close();
  delete File1;
  File2->Close();
  delete File2;
}

int main(int argc, char *argv[]) 
{
  SetMaCh3LoggerFormat();
  // Avoid Info in <TCanvas::Print>
  gErrorIgnoreLevel = kWarning;

  PlotMan = new MaCh3Plotting::PlottingManager();
  PlotMan->parseInputs(argc, argv);
  #ifdef DEBUG
  PlotMan->input().getFile(0).file->ls();
  #endif
  PlotMan->setExec("GetPostfitParamPlots");

  PlotMan->style().setPalette(PlotMan->getOption<std::string>("colorPalette"));

  GetPostfitParamPlots();
  GetViolinPlots();

  if (PlotMan->input().getNInputFiles() == 2)
  {
    std::string filename1 = PlotMan->getFileName(0);
    std::string filename2 = PlotMan->getFileName(1);
    Get2DComparison(filename1, filename2);
  }

  delete PlotMan;
  return 0;
}
