#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
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
#pragma GCC diagnostic pop

#include "plottingUtils/plottingUtils.h"
#include "plottingUtils/plottingManager.h"

/// @file GetPostfitParamPlots
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
/// @author Clarence Wret
/// @author Will Parker
/// @author Kamil Skwarczynski
/// @author Ewan Miller

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

MaCh3Plotting::PlottingManager *man;
TH1D *Prefit;

TH2D *Violin;
TH2D *Violin2;
TH2D *ViolinPre;

int NDParameters;
int NDParametersStartingPos;

int FDParameters;
int FDParametersStartingPos;

int OscParameters;
int OscParametersStartingPos;

int CrossSectionParameters;
int XsecStartingPos;

int FluxParameters;

std::vector<int> NDSamplesBins;
std::vector<std::string> NDSamplesNames;
int nBins;
TCanvas *canv;

std::string SaveName;

TPad *p3;
TPad *p4;

// KS: Color for 0 - prefit, 1 postfit, 2 another postfit, 3 you know the drill
Color_t PlotColor[] = {kRed, kBlack, kBlue, kGreen};
std::string plotType;

std::vector <TH1D *> PostfitHistVec;

void copyParToBlockHist(int localBin, std::string paramName, TH1D*blockHist, std::string type, int fileId, bool setLabels = true){
  // Set the values in the sub-histograms
  MACH3LOG_DEBUG("copying data from at local bin {}: for parameter {}", localBin, paramName);
  MACH3LOG_DEBUG("  Fitter specific name: {}", man->input().translateName(fileId, MaCh3Plotting::kPostFit, paramName));
  MACH3LOG_DEBUG("  value: {}", man->input().getPostFitValue(fileId, paramName, type));
  MACH3LOG_DEBUG("  error: {}", man->input().getPostFitError(fileId, paramName, type));

  blockHist->SetBinContent(localBin +1, man->input().getPostFitValue(fileId, paramName, type));
  blockHist->SetBinError(localBin +1, man->input().getPostFitError(fileId, paramName, type));

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
    int ParamBinId = -999;
    for (int ix = 0; ix < FullViolin->GetXaxis()->GetNbins(); ++ix) {
      if(FullViolin->GetXaxis()->GetBinLabel(ix+1) == ParamNames[i])
      {
        ParamBinId = ix+1;
        break;
      }
    }
    if(ParamBinId == -999) {
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
    title = man->style().prettifyParamName(title);
    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
  }
}

void PrettifyTitles(TH2D *Hist) {
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i) 
  {
    std::string title = Hist->GetXaxis()->GetBinLabel(i+1);

    title = man->style().prettifyParamName(title);
    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
  }
}

void ReadSettings(std::shared_ptr<TFile> File1)
{
  MACH3LOG_DEBUG("Reading settings for file {}", File1->GetName());
  #ifdef DEBUG
  File1->ls();
  #endif
  TTree *Settings = (File1->Get<TTree>("Settings"));
  MACH3LOG_DEBUG("Got settings tree");
  Settings->Print();

  Settings->SetBranchAddress("CrossSectionParameters", &CrossSectionParameters);
  MACH3LOG_DEBUG("XSec params: {}", CrossSectionParameters);
  Settings->SetBranchAddress("FluxParameters", &FluxParameters);
  Settings->SetBranchAddress("CrossSectionParametersStartingPos", &XsecStartingPos);
  Settings->SetBranchAddress("NDParameters", &NDParameters);
  Settings->SetBranchAddress("NDParametersStartingPos", &NDParametersStartingPos);
  Settings->SetBranchAddress("FDParameters", &FDParameters);
  Settings->SetBranchAddress("FDParametersStartingPos", &FDParametersStartingPos);
  Settings->SetBranchAddress("OscParameters", &OscParameters);
  Settings->SetBranchAddress("OscParametersStartingPos", &OscParametersStartingPos);

  std::vector<int> *NDSamples_Bins = 0;
  std::vector<std::string> *NDSamples_Names = 0;
  Settings->SetBranchAddress("NDSamplesNames", &NDSamples_Names);
  Settings->SetBranchAddress("NDSamplesBins", &NDSamples_Bins);

  Settings->GetEntry(0);

  NDSamplesNames = *NDSamples_Names;
  NDSamplesBins = *NDSamples_Bins;

  MACH3LOG_DEBUG("Read successfully");
}

inline TH1D* makeRatio(TH1D *PrefitCopy, TH1D *PostfitCopy, bool setAxes){
  // set up the ratio hist
  TH1D* Ratio = static_cast<TH1D*>(PrefitCopy->Clone());
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

      Ratio->SetLineWidth(man->getOption<int>("plotLineWidth"));
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

inline void DrawPlots(TCanvas *plotCanv, TH1D* PrefitCopy, const std::vector<TH1D*>& PostfitVec, TPad *mainPad, TPad *ratioPad) {
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
    TH1D *postFitHist = PostfitVec[fileId];
    
    postFitHist->SetMarkerColor(TColor::GetColorPalette(fileId));
    postFitHist->SetLineColor(TColor::GetColorPalette(fileId));
    postFitHist->SetMarkerStyle(7);
    postFitHist->SetLineStyle(1+fileId);
    postFitHist->SetLineWidth(man->getOption<int>("plotLineWidth"));

    postFitHist->Draw("e1, same");
  }
  
  plotCanv->Update();
  TGaxis *axis = new TGaxis(PrefitCopy->GetXaxis()->GetBinLowEdge(PrefitCopy->GetXaxis()->GetFirst()), gPad->GetUymin()+0.01,
                            PrefitCopy->GetXaxis()->GetBinLowEdge(PrefitCopy->GetXaxis()->GetFirst()), gPad->GetUymax(),
                            gPad->GetUymin()+0.01, gPad->GetUymax(), 510, "");
  axis->SetLabelFont(43);
  axis->SetLabelSize(25);
  axis->Draw();

  plotCanv->cd();
  ratioPad->Draw();
  ratioPad->cd();
  
  std::vector<TH1D*> ratioHists;

  // save pointers to these so we can delete them once we are done
  ratioHists.push_back(makeRatio(PrefitCopy, PostfitVec[0], true));

  ratioHists[0]->Draw("p");
  for(int postFitIdx = 1; postFitIdx < static_cast<int>(PostfitVec.size()); postFitIdx++){
    ratioHists.push_back(makeRatio(PrefitCopy, PostfitVec[postFitIdx], true));
    
    ratioHists[postFitIdx]->SetMarkerColor(TColor::GetColorPalette(postFitIdx));
    ratioHists[postFitIdx]->SetLineColor(TColor::GetColorPalette(postFitIdx));
    ratioHists[postFitIdx]->SetMarkerStyle(7);
    ratioHists[postFitIdx]->SetLineStyle(1+postFitIdx);
    ratioHists[postFitIdx]->SetLineWidth(man->getOption<int>("plotLineWidth"));

    ratioHists[postFitIdx]->Draw("p same");
  }

  // draw lines across the plot at +-1 and 0
  TLine line(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), 0.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), 0.0);
  TLine line2(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), 1.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), 1.0);
  TLine line3(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), -1.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), -1.0);
  
  line.SetLineColor(kRed);
  line.SetLineStyle(kDashed);
  line.SetLineWidth(man->getOption<int>("refLineWidth"));
  line2.SetLineColor(kRed);
  line2.SetLineStyle(kDashed);
  line2.SetLineWidth(man->getOption<int>("refLineWidth"));
  line3.SetLineColor(kRed);
  line3.SetLineStyle(kDashed);
  line3.SetLineWidth(man->getOption<int>("refLineWidth"));

  line.Draw("same");
  line2.Draw("same");
  line3.Draw("same");

  plotCanv->Print((SaveName).c_str());
  
  ratioHists.clear();
  delete axis;
}

void MakeXsecPlots()
{  
  // get the names of the blocks of parameters to group together
  std::vector<std::string> const blockNames = man->getOption<std::vector<std::string>>("paramGroups");
  const int XsecPlots = static_cast<int>(blockNames.size());

  for (int i = 0; i < XsecPlots; i++) 
  {
    // get the configuration for this parameter
    std::string blockName = blockNames[i];
    YAML::Node paramBlock = man->getOption(blockName);
    std::string blockTitle                 = paramBlock[0].as<std::string>();
    std::vector<double>      blockLimits   = paramBlock[1].as<std::vector<double>>();
    std::vector<std::string> blockContents = paramBlock[2].as<std::vector<std::string>>();

    // get num of params in the block
    const int nParams = static_cast<int>(blockContents.size());

    // set some plot things
    TH1D *blockHist_prefit = new TH1D(blockName.c_str(),
                                      blockTitle.c_str(), nParams, 0.0, static_cast<double>(nParams));
    
    man->style().setTH1Style(blockHist_prefit, man->getOption<std::string>("prefitHistStyle"));

    // set the errors for the prefit block hist
    for(int localBin=0; localBin < nParams; localBin ++){
      // the "local" bin is the params index within the group of parameters
      std::string paramName = blockContents[localBin];
      copyParToBlockHist(localBin, paramName, blockHist_prefit, "Prior", 0);
    }

    // now set for the postfit blocks for all files
    std::vector <TH1D *> blockHist_postfit_Vec;
    for(unsigned int fileId = 0; fileId < man->getNFiles(); fileId++){

      TH1D *blockHist_postfit = new TH1D((blockName + man->getFileName(fileId)).c_str(),
                                         blockTitle.c_str(), nParams, 0.0, static_cast<double>(nParams));

      // loop through all the parameters in this block and set the contents in the blocks TH1
      for(int localBin=0; localBin < nParams; localBin ++){
        // the "local" bin is the params index within the group of parameters
        std::string paramName = blockContents[localBin];
        copyParToBlockHist(localBin, paramName, blockHist_postfit, "", fileId);
      }

      blockHist_postfit_Vec.push_back(blockHist_postfit);
    } 

    // set the y axis limits we got from config
    blockHist_prefit->GetYaxis()->SetRangeUser(blockLimits[0], blockLimits[1]); 

    // Do some fancy replacements
    PrettifyTitles(blockHist_prefit);

    DrawPlots(canv, blockHist_prefit, blockHist_postfit_Vec, p3, p4);
    delete blockHist_prefit;
    for (auto hist : blockHist_postfit_Vec) {
      delete hist;
    }
    blockHist_postfit_Vec.clear();
  }
}

void MakeFluxPlots()
{
  // these for non named params where we dont need as much space
  TPad* p1 = new TPad("p1", "p1", 0.0, 0.3, 1.0, 1.0);
  TPad* p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 0.3);
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
  std::vector<std::string> const fluxBlockNames = man->getOption<std::vector<std::string>>("fluxGroups");
  auto const fluxBinningTable = man->getOption("FluxBinning");

  const int FluxPlots = static_cast<int>(fluxBlockNames.size());

  for (int i = 0; i < FluxPlots; i++) 
  {
    // get the configuration for this block
    std::string fluxBlockName = fluxBlockNames[i];
    YAML::Node paramBlock           = man->getOption(fluxBlockName);
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

    TH1D *blockHist_prefit = new TH1D(fluxBlockName.c_str(), blockTitle.c_str(), nParams, binning.data());
    blockHist_prefit->GetYaxis()->SetTitle("Parameter Variation");
    blockHist_prefit->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    blockHist_prefit->GetXaxis()->SetTitleOffset(blockHist_prefit->GetXaxis()->GetTitleOffset()*1.2);
    man->style().setTH1Style(blockHist_prefit, man->getOption<std::string>("prefitHistStyle"));
    // set the errors for the prefit block hist
    for(int fluxParId = blockContents[0]; fluxParId <= blockContents[1]; fluxParId++){
        int localBin = fluxParId - blockContents[0];
        std::string paramName = "b_" + std::to_string(fluxParId);
        copyParToBlockHist(localBin, paramName, blockHist_prefit, "Prior", 0, false);
    }

    // now set for the postfit blocks for all files
    std::vector <TH1D *> blockHist_postfit_Vec;
    for(unsigned int fileId = 0; fileId < man->getNFiles(); fileId++){
      TH1D *blockHist_postfit = new TH1D(fluxBlockName.c_str(), blockTitle.c_str(), nParams, binning.data());

      for(int fluxParId = blockContents[0]; fluxParId <= blockContents[1]; fluxParId++){
        int localBin = fluxParId - blockContents[0];
        std::string paramName = "b_" + std::to_string(fluxParId);

        copyParToBlockHist(localBin, paramName, blockHist_postfit, "", fileId, false);
      }

      blockHist_postfit_Vec.push_back(blockHist_postfit);
    }
    
    // set the y axis limits we got from config
    blockHist_prefit->GetYaxis()->SetRangeUser(blockLimits[0], blockLimits[1]); 

    DrawPlots(canv, blockHist_prefit, blockHist_postfit_Vec, p1, p2);
    delete blockHist_prefit;
    for (auto hist : blockHist_postfit_Vec) {
      delete hist;
    }
    blockHist_postfit_Vec.clear();
  }

  canv->cd();
  canv->SetLogx(false);
  canv->SetBottomMargin(canv->GetBottomMargin()*1.7);
  delete p1;
  delete p2;
}

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

    std::vector<TH1D*> PostfitNDDetHistVec(man->getNFiles());
    TH1D *PreFitNDDetHist = man->input().getFile(0).file->Get<TH1D>(Form("param_%s_prefit", NDSamplesNames[i].c_str()));
    man->style().setTH1Style(PreFitNDDetHist, man->getOption<std::string>("prefitHistStyle"));

    std::string temp = NDSamplesNames[i].c_str();
    while (temp.find("_") != std::string::npos) {
       temp.replace(temp.find("_"), 1, std::string(" "));
     }
    PreFitNDDetHist->SetTitle(temp.c_str());
    PreFitNDDetHist->GetXaxis()->SetRangeUser(Start, NDbinCounter);
    
    MACH3LOG_DEBUG("  Start bin: {} :: End bin: {}", Start, NDbinCounter);
    // set the x range for the postfits
    for(unsigned int fileId = 0; fileId < man->getNFiles(); fileId++){
      PostfitNDDetHistVec[fileId] = man->input().getFile(fileId).file->Get<TH1D>(Form("param_%s_%s", NDSamplesNames[i].c_str(), plotType.c_str()));
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

    PreFitNDDetHist->GetYaxis()->SetRangeUser(man->getOption<double>("detParYRange_low"), man->getOption<double>("detParYRange_high"));
    
    Start += NDSamplesBins[i];

    DrawPlots(canv, PreFitNDDetHist, PostfitNDDetHistVec, pTop, pDown);
    canv->Update();
  }
  delete pTop;
  delete pDown;
}

void MakeFDDetPlots()
{
  Prefit->GetYaxis()->SetTitleOffset(Prefit->GetYaxis()->GetTitleOffset()*1.2);
  
  // these for named parameters where we need a nice big gap at the bottom to fit the names
  TPad* pTop = nullptr;
  TPad* pDown = nullptr;
  InitializePads(canv, pTop, pDown);

  int FDbinCounter = FDParametersStartingPos;
  int Start = FDbinCounter;
  //KS: WARNING This is hardcoded
  double FDSamplesBins[7] = {12,18,30,36,44,57,58};
  std::string FDSamplesNames[7] = {"FHC 1Re", "FHC 1R#mu", "RHC 1Re","RHC 1R#mu","FHC 1Re 1 d.e.","FHC MR#mu 1 or 2 d.e.", "Momentum Scale"};
  for (unsigned int i = 0; i < 7; ++i)
  {
    FDbinCounter = FDParametersStartingPos + FDSamplesBins[i];

    canv->cd();
    Prefit->SetTitle(FDSamplesNames[i].c_str());
    Prefit->GetXaxis()->SetRangeUser(Start, FDbinCounter);
    //KS: We don't' need name for every FD param
    for(int j = 0; j < FDSamplesBins[i]; ++j)
    {
      //bool ProductOfTen = false;
      //if(j % 10) ProductOfTen = true;
      //if(j != 0 && ProductOfTen) Prefit->GetXaxis()->SetBinLabel(Start+j+1, " ");
    }
    Prefit->GetYaxis()->SetRangeUser(man->getOption<double>("detParYRange_low"), man->getOption<double>("detParYRange_high"));
    Prefit->Draw("e2");
    
    Start = FDParametersStartingPos + FDSamplesBins[i];

    DrawPlots(canv, Prefit, PostfitHistVec, pTop, pDown);
    canv->Update();
   }
   delete pTop;
   delete pDown;
}

void MakeOscPlots()
{
  // these for named parameters where we need a nice big gap at the bottom to fit the names
  TPad* pTop = nullptr;
  TPad* pDown = nullptr;
  InitializePads(canv, pTop, pDown);

  Prefit->GetYaxis()->SetTitleOffset(Prefit->GetYaxis()->GetTitleOffset()*1.2);
  
  canv->cd();
  Prefit->SetTitle("Oscillation Parameters");
  Prefit->GetXaxis()->SetRangeUser(OscParametersStartingPos, OscParametersStartingPos+OscParameters);

  Prefit->GetYaxis()->SetRangeUser(man->getOption<double>("oscParYRange_low"), man->getOption<double>("oscParYRange_high"));
  Prefit->Draw("e2");

  DrawPlots(canv, Prefit, PostfitHistVec, pTop, pDown);
  delete pTop;
  delete pDown;
}

void MakeXsecRidgePlots()
{  
  gStyle->SetPalette(51);

  TCanvas *blankCanv = new TCanvas ("blankCanv", "blankCanv", 2048, 2048);
  blankCanv->SaveAs("RidgePlots.pdf[");

  // get the names of the blocks of parameters to group together
  std::vector<std::string> const blockNames = man->getOption<std::vector<std::string>>("paramGroups");
  const int XsecPlots = static_cast<int>(blockNames.size());

  double padTopMargin = 0.9;
  double padBottomMargin = 0.1;
  double padOverlap = 0.9;
  double ridgeLineWidth = 1.0;

  for (int i = 0; i < XsecPlots; i++) 
  {
    // get the configuration for this parameter
    std::string blockName = blockNames[i];
    auto const &paramBlock = man->getOption(blockName);
    auto blockTitle    = paramBlock[0].as<std::string>();
    auto blockLimits   = paramBlock[1].as<std::vector<double>>();
    auto blockContents = paramBlock[2].as<std::vector<std::string>>();

    // the directory of histograms
    TDirectoryFile *posteriorDir = man->input().getFile(0).file->Get<TDirectoryFile>("Post_1d_hists");

    // get num of params in the block
    int nParams = static_cast<int>(blockContents.size());
    TCanvas *ridgeCanv = new TCanvas ("RidgePlotCanv", "RidgePlotCanv", 2048, 2048);
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
    TH1D *axisPlot = new TH1D("axis plot", "", 1, blockLimits[0], blockLimits[1]);

    for(int parId = 0; parId < nParams; parId++) {
      std::string paramName = blockContents[parId];

      TH1D *posteriorDist = nullptr;
      // get the list of objects in the directory
      TIter next(posteriorDir->GetListOfKeys());
      while (TKey* key = static_cast<TKey*>(next())) {
        // check if the end of the param name matches with the MaCh3 name, do this so we exclude things like nds_ at the start of the name
        std::string str(key->GetTitle());
        std::string name = man->input().translateName(0, MaCh3Plotting::kPostFit, paramName);
        uint pos = str.find(name);
        bool foundPar = (pos == str.length() - name.length());

        if(foundPar){
          posteriorDist = posteriorDir->Get<TH1D>(key->GetName());
        }
      }

      if(posteriorDist == nullptr){
        MACH3LOG_WARN("Couldn't find parameter {} when making ridgeline plots", paramName);
        continue;
      }
      
      // EM: do some funky scaling so that we always get evenly spaced pads in the range [bottomMargin, TopMargin] with the specified overlap
      double padAnchor = padBottomMargin + (static_cast<double>(nParams - parId - 1) /
                                  static_cast<double>(nParams - 1)) * (padTopMargin - padBottomMargin);
      double padWidth = (padTopMargin - padBottomMargin) / static_cast<double>(nParams);
      double norm = (padTopMargin - padBottomMargin);

      double padTop = padWidth * (1.0 + padOverlap) * (padTopMargin - padAnchor) / norm + padAnchor;
      double padBottom = padAnchor - padWidth * (1.0 + padOverlap) * (padAnchor - padBottomMargin) / norm;

      TPad *pad = new TPad(paramName.c_str(), "", 0.3, padBottom, 0.9, padTop, -1, 0, -1); 
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
      
      TH1D* axisPlot_tmp = static_cast<TH1D*>(axisPlot->Clone(Form("AxisPlot_%s", paramName.c_str())));
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
      label->DrawLatexNDC(0.29, padBottom + 0.005, man->style().prettifyParamName(paramName).c_str());
      line->DrawLine(0.1, padBottom, 0.9, padBottom);
    }

    ridgeCanv->cd();
    ridgeCanv->SetGrid(1,1);
    TPad *axisPad = new TPad("AxisPad", "", 0.3, 0.0, 0.9, 1.0, -1, 0, -1);
    axisPad->SetLeftMargin(0.0);
    axisPad->SetRightMargin(0.0);
    axisPad->Draw();
    axisPad->cd();
    axisPad->SetGrid(1,1);
    axisPad->SetFrameFillStyle(4000);
    
    axisPlot->GetXaxis()->SetTickSize(0.01);
    axisPlot->GetXaxis()->SetTitle("");
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
    delete ridgeCanv;
    delete axisPlot;
  }
  blankCanv->SaveAs("RidgePlots.pdf]");
  delete blankCanv;
}

void GetPostfitParamPlots()
{
  SaveName = man->getOutputName();
  
  //KS: By default we take HPD values, but by changing "plotType" you can use for example Gauss
  plotType = "HPD";
  //plotType = "gaus"; 
    
  MACH3LOG_INFO("Plotting {} errors", plotType);
  
  ReadSettings(man->input().getFile(0).file);

  canv = new TCanvas("canv", "canv", 1024, 1024);
  //gStyle->SetPalette(51);
  gStyle->SetOptStat(0); //Set 0 to disable statystic box
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
  man->style().setTH1Style(Prefit, man->getOption<std::string>("prefitHistStyle"));
  leg->AddEntry(Prefit, "Prior", "lpf");

  for(unsigned int fileId = 0; fileId < man->getNFiles(); fileId++){
    TH1D *postFitHist_tmp = new TH1D();
    postFitHist_tmp->SetBit(kCanDelete);
    
    postFitHist_tmp->SetMarkerColor(TColor::GetColorPalette(fileId));
    postFitHist_tmp->SetLineColor(TColor::GetColorPalette(fileId));
    postFitHist_tmp->SetMarkerStyle(7);
    postFitHist_tmp->SetLineStyle(1+fileId);
    postFitHist_tmp->SetLineWidth(man->getOption<int>("plotLineWidth"));
    leg->AddEntry(postFitHist_tmp, man->getFileLabel(fileId).c_str(), "lpf");
  }
  /// @todo this is temporary hack
  for(unsigned int fileId = 0; fileId < man->getNFiles(); fileId++)
  {
    TH1D *Hist = static_cast<TH1D*>((man->input().getFile(fileId).file->Get( ("param_xsec_"+plotType).c_str()))->Clone());

    Hist->SetMarkerColor(TColor::GetColorPalette(fileId));
    Hist->SetLineColor(TColor::GetColorPalette(fileId));
    Hist->SetMarkerStyle(7);
    Hist->SetLineStyle(1+fileId);
    Hist->SetLineWidth(man->getOption<int>("plotLineWidth"));
    PostfitHistVec.push_back(Hist);
  }

  canv->cd();
  canv->Clear();
  leg->Draw();
  canv->Print((SaveName).c_str());

  MakeXsecPlots();

  MakeFluxPlots();

  //KS: By default we don't run ProcessMCMC with PlotDet as this take some time, in case we did let's make fancy plots
  if(NDParameters > 0) MakeNDDetPlots();
  
  //KS: Same as above but for FD parameters,
  if(FDParameters > 0) MakeFDDetPlots();
  
  if(OscParameters > 0) MakeOscPlots();

  canv->Print((SaveName+"]").c_str());

  MakeXsecRidgePlots();

  delete canv;
  delete Prefit;
}

inline TGraphAsymmErrors* MakeTGraphAsymmErrors(std::shared_ptr<TFile> File,  std::vector<int> Index = {})
{
  int GraphBins = Index.size() == 0 ? nBins : Index.size();
  std::vector<double> x(GraphBins);
  std::vector<double> y(GraphBins);
  std::vector<double> exl(GraphBins);
  std::vector<double> eyl(GraphBins);
  std::vector<double> exh(GraphBins);
  std::vector<double> eyh(GraphBins);

  TH1D* PostHist = static_cast<TH1D*>(File->Get( ("param_xsec_"+plotType).c_str() ));

  TVectorD* Errors_HPD_Positive = static_cast<TVectorD*>(File->Get( "Errors_HPD_Positive" ));
  TVectorD* Errors_HPD_Negative = static_cast<TVectorD*>(File->Get( "Errors_HPD_Negative" ));
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
  TGraphAsymmErrors* PostGraph = new TGraphAsymmErrors(GraphBins, x.data(), y.data(), exl.data(), exh.data(), eyl.data(), eyh.data());
  PostGraph->SetTitle("");

  return PostGraph;
}
   
//KS: Make fancy violin plots
void GetViolinPlots(std::string FileName1 = "", std::string FileName2 = "")
{
  //KS: Should be in some config... either way it control whether you plot symmetric or asymmetric error bars
  bool PlotAssym = true;
    
  //KS: No idea why but ROOT changed treatment of violin in R6. If you have non uniform binning this will results in very hard to see violin plots.
  TCandle::SetScaledViolin(false);

  MACH3LOG_INFO("Making Violin Plot");
  if (!FileName1.empty()) MACH3LOG_INFO("File 1: {} ", FileName1);
  if (!FileName2.empty()) MACH3LOG_INFO("File 2: {}", FileName2);
    
  SaveName = FileName1;
  SaveName = SaveName.substr(0, SaveName.find(".root"));
  if(FileName2 != "") SaveName += FileName2;
  if(FileName2 != "") SaveName = SaveName.substr(0, SaveName.find(".root"));
  SaveName += "_Violin";
  if(PlotAssym) SaveName += "_Assym";
   
  std::shared_ptr<TFile> File1 = std::make_shared<TFile>(FileName1.c_str());
  std::shared_ptr<TFile> File2 = nullptr;
  if(FileName2 != "") File2 = std::make_shared<TFile>(FileName2.c_str());
  
  canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetGrid();
  gStyle->SetOptStat(0);
  //KS: Remove errors on X axis as they are confusing in violin type of plot
  if(!PlotAssym) gStyle->SetErrorX(0.0001);
  canv->SetTickx();
  canv->SetTicky();
  canv->SetBottomMargin(0.25);
  canv->SetTopMargin(0.08);
  canv->SetRightMargin(0.03);
  canv->SetLeftMargin(0.10);
  canv->Print((SaveName+".pdf[").c_str());
  canv->SetGrid();

  Violin = nullptr;
  ViolinPre = File1->Get<TH2D>( "param_violin_prior" );
  Violin = File1->Get<TH2D>( "param_violin" );
  if(Violin == nullptr)
  {
    MACH3LOG_ERROR("Couldn't find violin plot, make sure method from MCMCProcessor is being called");
    return;
  }
  // Do some fancy replacements
  ViolinPre->SetFillColor(kRed);
  ViolinPre->SetFillColorAlpha(kRed, 0.35);
  ViolinPre->SetMarkerColor(kRed);
  ViolinPre->SetMarkerStyle(20);
  ViolinPre->SetMarkerSize(0.5);

  ViolinPre->GetYaxis()->SetTitleOffset(1.3);
  ViolinPre->GetYaxis()->SetTitle("Parameter Value");
  ViolinPre->GetXaxis()->LabelsOption("v");

  Violin->SetFillColor(kBlue);
  Violin->SetFillColorAlpha(kBlue, 0.35);
  Violin->SetMarkerColor(kBlue);
  Violin->SetMarkerStyle(20);
  Violin->SetMarkerSize(0.5);
    
  TH1D* Postfit = File1->Get<TH1D>( ("param_xsec_"+plotType).c_str() );
  Postfit->SetMarkerColor(kRed);
  Postfit->SetLineColor(kRed);
  Postfit->SetMarkerStyle(7);
  
 if(File2 != nullptr)
  {
    Violin2 = File2->Get<TH2D>( "param_violin" );
    Violin2->SetMarkerColor(kGreen);
    Violin2->SetLineColor(kGreen);
    Violin2->SetFillColor(kGreen);
    Violin2->SetFillColorAlpha(kGreen, 0.35);
  }
  
  TGraphAsymmErrors* PostGraphAll = MakeTGraphAsymmErrors(File1);

  // Make a Legend page
  auto leg = std::make_unique<TLegend>(0.0, 0.0, 1.0, 1.0);
  if (ViolinPre != nullptr) leg->AddEntry(ViolinPre, "Prior", "lpf");
  if (Violin != nullptr)    leg->AddEntry(Violin, "Posterior", "lpf");
  if (Violin2 != nullptr)   leg->AddEntry(Violin2, "Second Violin", "lpf");
  if(PlotAssym)             leg->AddEntry(PostGraphAll, "HPD Assym", "lp");
  else                      leg->AddEntry(Postfit, "HPD", "lpf");

  canv->cd();
  canv->Clear();
  leg->Draw();
  canv->Print((SaveName+".pdf").c_str());
  
  // get the names of the blocks of parameters to group together
  std::vector<std::string> const blockNames = man->getOption<std::vector<std::string>>("paramGroups");
  const int XsecPlots = static_cast<int>(blockNames.size());

  for (int i = 0; i < XsecPlots; i++)
  {
    // get the configuration for this parameter
    std::string blockName = blockNames[i];
    YAML::Node paramBlock = man->getOption(blockName);
    std::string blockTitle                 = paramBlock[0].as<std::string>();
    std::vector<double>      blockLimits   = paramBlock[1].as<std::vector<double>>();
    std::vector<std::string> blockContents = paramBlock[2].as<std::vector<std::string>>();

    // get num of params in the block
    const int nParams = static_cast<int>(blockContents.size());

    // set some plot things
    auto blockHist_prefit = std::make_unique<TH2D>(blockName.c_str(), blockTitle.c_str(), nParams, 0.0, static_cast<double>(nParams),
                                      ViolinPre->GetYaxis()->GetNbins(), ViolinPre->GetYaxis()->GetXmin(), ViolinPre->GetYaxis()->GetXmax());
    CopyViolinToBlock(ViolinPre, blockHist_prefit.get(), blockContents);
    // set the y axis limits we got from config
    blockHist_prefit->GetYaxis()->SetRangeUser(blockLimits[0], blockLimits[1]);
    auto blockHist_Violin1 = std::make_unique<TH2D>((blockTitle + "Violin1").c_str(), (blockTitle + "Violin1").c_str(), nParams, 0.0, static_cast<double>(nParams),
                                       Violin->GetYaxis()->GetNbins(), Violin->GetYaxis()->GetXmin(), Violin->GetYaxis()->GetXmax());
    CopyViolinToBlock(Violin, blockHist_Violin1.get(), blockContents);
    blockHist_Violin1->GetYaxis()->SetRangeUser(blockLimits[0], blockLimits[1]);
    std::unique_ptr<TH2D> blockHist_Violin2 = nullptr;
    if(Violin2 != nullptr) {
      blockHist_Violin2 = std::make_unique<TH2D>((blockTitle + "Violin2").c_str(), (blockTitle + "Violin2").c_str(), nParams, 0.0, static_cast<double>(nParams),
                                   Violin2->GetYaxis()->GetNbins(), Violin2->GetYaxis()->GetXmin(), Violin2->GetYaxis()->GetXmax());
      CopyViolinToBlock(Violin2, blockHist_Violin2.get(), blockContents);
      blockHist_Violin2->GetYaxis()->SetRangeUser(blockLimits[0], blockLimits[1]);
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
      int ParamBinId = -999;
      for (int ix = 0; ix < ViolinPre->GetXaxis()->GetNbins(); ++ix) {
        if(ViolinPre->GetXaxis()->GetBinLabel(ix+1) == blockContents[is]) {
          ParamBinId = ix;
          break;
        }
      }
      Index.push_back(ParamBinId);
    }
    TGraphAsymmErrors* PostGraph = MakeTGraphAsymmErrors(File1, Index);
    PostGraph->SetMarkerColor(kBlack);
    PostGraph->SetLineColor(kBlack);
    PostGraph->SetMarkerStyle(7);
    PostGraph->SetLineWidth(2);
    PostGraph->SetLineStyle(kSolid);

    blockHist_prefit->Draw("violinX(03100300)");
    blockHist_Violin1->Draw("violinX(03100300) SAME");
    if(blockHist_Violin2 != nullptr) {
      blockHist_Violin2->Draw("violinX(03100300) SAME");
    }
    if(PlotAssym) PostGraph->Draw("P SAME");
    else Postfit->Draw("SAME");
    canv->Print((SaveName+".pdf").c_str());
    delete PostGraph;
  }

  canv->Print((SaveName+".pdf]").c_str());
  delete canv;
  delete ViolinPre;
  delete Violin;
  delete PostGraphAll;
  if(Violin2 != nullptr) delete Violin2;
  delete Postfit;
  File1->Close();
  if(File2 != nullptr) {
    File2->Close();
  }
}

int main(int argc, char *argv[]) 
{
  SetMaCh3LoggerFormat();
  // Avoid Info in <TCanvas::Print>
  gErrorIgnoreLevel = kWarning;

  man = new MaCh3Plotting::PlottingManager();
  man->parseInputs(argc, argv);
  #ifdef DEBUG
  std::cout << std::endl << std::endl << "====================" << std::endl;
  man->input().getFile(0).file->ls();
  #endif
  man->setExec("GetPostfitParamPlots");

  man->style().setPalette(man->getOption<std::string>("colorPalette"));

  GetPostfitParamPlots();

  if (argc != 2 && argc != 3 && argc !=4)
  {
    MACH3LOG_CRITICAL("Invalid command line options specified");
    MACH3LOG_CRITICAL("How to use: {} <MCMC_Processor_Output>.root", argv[0]);
    MACH3LOG_CRITICAL("You can add up to 3 different files");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (argc == 2)
  {
    std::string filename = argv[1];
    GetViolinPlots(filename);
  }
  else if (argc == 3)
  {
    std::string filename1 = argv[1];
    std::string filename2 = argv[2];
    GetViolinPlots(filename1, filename2);
  }
  else if (argc == 4)
  {
    std::string filename1 = argv[1];
    std::string filename2 = argv[2];
    std::string filename3 = argv[3];
    //KS: Violin plot currently not supported by three file version although it should be super easy to adapt
  }

  delete man;
  return 0;
}
