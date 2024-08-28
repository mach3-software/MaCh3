#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

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

#include "plottingUtils/plottingUtils.h"
#include "plottingUtils/plottingManager.h"

//How to use ./GetPostfitParamPlots ProcessMCMC_Output1.root <ProcessMCMC_Output2.root> <ProcessMCMC_Output3.root>
//Originally written by Clarence, with some changes by Will, updated by Kamil, made generic plotter by Ewan
//Central postfit value taken is the Highest Posterior Density but can be changed easily to Gaussian etc. Watch out for parameter names and number of parameters per plot being quite hardcoded

//g++ `root-config --cflags` -g -std=c++11 -o GetPostfitParamPlots GetPostfitParamPlots.cpp -I`root-config --incdir` `root-config --glibs --libs`

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
TPad *p1;
TPad *p2;

TPad *p3;
TPad *p4;

// KS: Color for 0 - prefit, 1 postfit, 2 another postfit, 3 you know the drill
Color_t PlotColor[] = {kRed, kBlack, kBlue, kGreen};
std::string plotType;

std::vector <TH1D *> PostfitHistVec;

void copyParToBlockHist(int localBin, std::string paramName, TH1D*blockHist, std::string type, int fileId, bool setLabels = true){
  // Set the values in the sub-histograms
  MACH3LOG_DEBUG("copyin data from at local bin {}: for parameter {}", localBin, paramName);
  MACH3LOG_DEBUG("  Fitter specific name: {}", man->input().TranslateName(fileId, MaCh3Plotting::kPostFit, paramName));
  MACH3LOG_DEBUG("  value: {}", man->input().GetPostFitValue(fileId, paramName, type));
  MACH3LOG_DEBUG("  error: {}", man->input().GetPostFitError(fileId, paramName, type));

  blockHist->SetBinContent(localBin +1, man->input().GetPostFitValue(fileId, paramName, type));
  blockHist->SetBinError(localBin +1, man->input().GetPostFitError(fileId, paramName, type));

  if(setLabels){
    blockHist->GetXaxis()->SetBinLabel(localBin +1, paramName.c_str());
    blockHist->GetXaxis()->LabelsOption("v");
  }
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
  File1->ls();
  TTree *Settings = (TTree*)(File1->Get("Settings"));
  MACH3LOG_DEBUG("Got settings tree");
  Settings->Print();

  std::cout << CrossSectionParameters << std::endl;

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
  TH1D *Ratio = (TH1D*)PrefitCopy->Clone();
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

inline void DrawPlots(TCanvas *plotCanv, TH1D* PrefitCopy, std::vector<TH1D *>PostfitVec, TPad *mainPad, TPad *ratioPad) {
  // Draw!
  plotCanv->cd();
  mainPad->Draw();
  mainPad->cd();
  PrefitCopy->GetYaxis()->SetTitle("Parameter Value");
  
  PrefitCopy->GetYaxis()->SetLabelSize(0.);
  PrefitCopy->GetYaxis()->SetTitleSize(0.05);
  PrefitCopy->GetYaxis()->SetTitleOffset(1.3);
  PrefitCopy->Draw("e2");

  for(int fileId=0; fileId < (int)PostfitVec.size(); fileId++){
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
  for(int postFitIdx = 1; postFitIdx < (int)PostfitVec.size(); postFitIdx++){
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
  const int XsecPlots = (int)blockNames.size();

  for (int i = 0; i < XsecPlots; i++) 
  {
    // get the configuration for this parameter
    std::string blockName = blockNames[i];
    YAML::Node paramBlock = man->getOption(blockName);
    std::string blockTitle                 = paramBlock[0].as<std::string>();
    std::vector<double>      blockLimits   = paramBlock[1].as<std::vector<double>>();
    std::vector<std::string> blockContents = paramBlock[2].as<std::vector<std::string>>();

    // get num of params in the block
    int nParams = (int)blockContents.size();

    // set some plot things
    TH1D *blockHist_prefit = new TH1D(blockName.c_str(), blockTitle.c_str(), nParams, 0.0, (double)nParams);
    
    man->style().setTH1Style(blockHist_prefit, man->getOption<std::string>("prefitHistStyle"));

    // set the errors for the prefit block hist
    for(int localBin=0; localBin < nParams; localBin ++){
      // the "local" bin is the params index within the group of parameters
      std::string paramName = blockContents[localBin];
      copyParToBlockHist(localBin, paramName, blockHist_prefit, "Prior", 0);
    }

    // now set for the postfit blocks for all files
    std::vector <TH1D *> blockHist_postfit_Vec;
    for(int fileId = 0; fileId < man->getNFiles(); fileId++){

      TH1D *blockHist_postfit = new TH1D((blockName + man->getFileName(fileId)).c_str(), blockTitle.c_str(), nParams, 0.0, (double)nParams);

      // loop throught all the parameters in this block and set the contents in the blocks TH1
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
    blockHist_postfit_Vec.clear();
   }
}

void MakeFluxPlots()
{
  p1->SetLogx(true);
  p2->SetLogx(true);
  
  // get the names of the blocks of parameters to group together
  std::vector<std::string> const fluxBlockNames = man->getOption<std::vector<std::string>>("fluxGroups");
  auto const fluxBinningTable = man->getOption("FluxBinning");

  const int FluxPlots = (int)fluxBlockNames.size();

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
    if(nParams != (int)binning.size() -1){
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
    for(int fileId = 0; fileId < man->getNFiles(); fileId++){
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
    blockHist_postfit_Vec.clear();
  }

  canv->cd();
  canv->SetLogx(false);
  canv->SetBottomMargin(canv->GetBottomMargin()*1.7);
}

void MakeNDDetPlots()
{
  Prefit->GetYaxis()->SetTitleOffset(Prefit->GetYaxis()->GetTitleOffset()*1.2);

  int NDbinCounter = NDParametersStartingPos;
  int Start = NDbinCounter;
  
  MACH3LOG_INFO("Running on {} samples", NDSamplesNames.size());

  for (unsigned int i = 0; i < NDSamplesNames.size(); ++i)
  {
    MACH3LOG_DEBUG("--- On sample {}", NDSamplesNames[i]); 
    NDbinCounter += NDSamplesBins[i];

    std::vector<TH1D*> PostfitNDDetHistVec(man->getNFiles());
    TH1D *PreFitNDDetHist = (TH1D*)man->input().GetFile(0).file->Get(Form("param_%s_prefit", NDSamplesNames[i].c_str()));
    man->style().setTH1Style(PreFitNDDetHist, man->getOption<std::string>("prefitHistStyle"));

    std::string temp = NDSamplesNames[i].c_str();
    while (temp.find("_") != std::string::npos) {
       temp.replace(temp.find("_"), 1, std::string(" "));
     }
    PreFitNDDetHist->SetTitle(temp.c_str());
    PreFitNDDetHist->GetXaxis()->SetRangeUser(Start, NDbinCounter);
    
    MACH3LOG_DEBUG("  Start bin: {} :: End bin: {}", Start, NDbinCounter);
    // set the x range for the postfits
    for(int fileId = 0; fileId < man->getNFiles(); fileId++){
      PostfitNDDetHistVec[fileId] = (TH1D*)man->input().GetFile(fileId).file->Get(Form("param_%s_%s", NDSamplesNames[i].c_str(), plotType.c_str()));
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

    DrawPlots(canv, PreFitNDDetHist, PostfitNDDetHistVec, p3, p4);
    canv->Update();
  }
}

void MakeFDDetPlots()
{
  Prefit->GetYaxis()->SetTitleOffset(Prefit->GetYaxis()->GetTitleOffset()*1.2);
  
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
    //KS: We dont' need name for every FD param
    for(int j = 0; j < FDSamplesBins[i]; ++j)
    {
      //bool ProductOfTen = false;
      //if(j % 10) ProductOfTen = true;
      //if(j != 0 && ProductOfTen) Prefit->GetXaxis()->SetBinLabel(Start+j+1, " ");
    }
    Prefit->GetYaxis()->SetRangeUser(man->getOption<double>("detParYRange_low"), man->getOption<double>("detParYRange_high"));
    Prefit->Draw("e2");
    
    Start = FDParametersStartingPos + FDSamplesBins[i];

    DrawPlots(canv, Prefit, PostfitHistVec, p3, p4);
    canv->Update();
   }
}

void MakeOscPlots()
{
  Prefit->GetYaxis()->SetTitleOffset(Prefit->GetYaxis()->GetTitleOffset()*1.2);
  
  canv->cd();
  Prefit->SetTitle("Oscilation Parameters");
  Prefit->GetXaxis()->SetRangeUser(OscParametersStartingPos, OscParametersStartingPos+OscParameters);

    Prefit->GetYaxis()->SetRangeUser(man->getOption<double>("oscParYRange_low"), man->getOption<double>("oscParYRange_high"));
  Prefit->Draw("e2");

  DrawPlots(canv, Prefit, PostfitHistVec, p3, p4);
}


void MakeXsecRidgePlots()
{  
  
  
  gStyle->SetPalette(51);

  TCanvas *blankCanv = new TCanvas ("blankCanv", "blankCanv", 2048, 2048);
  blankCanv->SaveAs("RidgePlots.pdf[");

  // get the names of the blocks of parameters to group together
  std::vector<std::string> const blockNames = man->getOption<std::vector<std::string>>("paramGroups");
  const int XsecPlots = (int)blockNames.size();

  double padTopMargin = 0.9;
  double padBottomMargin = 0.1;
  double padOverlap = 0.9;
  double ridgeLineWidth = 1.0;

  for (int i = 0; i < XsecPlots; i++) 
  {
    
    // get the configuration for this parameter
    std::string blockName = blockNames[i];
    auto const &paramBlock = man->getOption(blockName);
    std::string blockTitle                 = paramBlock[0].as<std::string>();
    std::vector<double>      blockLimits   = paramBlock[1].as<std::vector<double>>();
    std::vector<std::string> blockContents = paramBlock[2].as<std::vector<std::string>>();

    // the directory of histograms
    TDirectoryFile *posteriorDir = (TDirectoryFile *)man->input().GetFile(0).file->Get("Post");

    // get num of params in the block
    int nParams = (int)blockContents.size();
    TCanvas *ridgeCanv = new TCanvas ("RidgePlotCanv", "RidgePlotCanv", 2048, 2048);
    ridgeCanv->Divide(1,1+nParams, 0.01, 0.0);

    TLatex *title = new TLatex();
    title->SetTextAlign(21);
    title->SetTextSize(0.03);
    title->DrawLatex(0.5, 0.95, blockTitle.c_str());

    TLatex *label = new TLatex();
    label->SetTextAlign(31);
    label->SetTextSize(0.02);
    
    TLine *line = new TLine();
    line->SetLineColor(kBlack);
    line->SetLineWidth(ridgeLineWidth);

    // use this to set the limits and also to plot the x axis and grid
    TH1D *axisPlot = new TH1D("axis plot", "", 1, blockLimits[0], blockLimits[1]);

    for(int parId=0; parId < nParams; parId++){
      std::string paramName = blockContents[parId];

      TCanvas *posteriorDistCanv = NULL;
      TH1D *posteriorDist = NULL;

      // get the list of objects in the directory
      TIter next(posteriorDir->GetListOfKeys());
      while(TKey *key = (TKey*) next()){ 
        // check if the end of the param name matches with the MaCh3 name, do this so we exclude things like nds_ at the start of the name
        std::string str(key->GetTitle());
        std::string name = man->input().TranslateName(0, MaCh3Plotting::kPostFit, paramName);
        uint pos = str.find(name);
        bool foundPar = (pos == str.length() - name.length());

        if(foundPar){
          posteriorDistCanv = (TCanvas*)posteriorDir->Get(key->GetName());
          posteriorDist = (TH1D*)posteriorDistCanv->GetPrimitive(key->GetName());
        }
      }

      if(posteriorDist == NULL){
        MACH3LOG_WARN("Couldnt find parameter {} when making ridgeline plots", paramName);
        continue;
      }
      
      // EM: do some funky scaling so that we always get evenly spaced pads in the range [bottomMargin, TopMargin] with the specified overlap
      double padAnchor = padBottomMargin + ((float)(nParams - parId -1) / (float)(nParams -1)) * (padTopMargin - padBottomMargin);
      double padWidth = ((padTopMargin - padBottomMargin) / (float)(nParams));
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
      
      TH1D *axisPlot_tmp = (TH1D*)axisPlot->Clone(Form("AxisPlot_%s", paramName.c_str()));
      axisPlot_tmp->Draw("A");
      posteriorDist->Draw("H SAME");
    
      axisPlot_tmp->GetYaxis()->SetRangeUser(0.0, 0.7 *posteriorDist->GetMaximum());
      posteriorDist->SetLineColor(kWhite);
      posteriorDist->SetFillColorAlpha(TColor::GetColorPalette(floor((float)parId * TColor::GetNumberOfColors()/ (float)nParams)), 0.85);
      
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
  }
  
  blankCanv->SaveAs("RidgePlots.pdf]");
}

void GetPostfitParamPlots()
{
  SaveName = man->getOutputName();
  
  //KS: By default we take HPD values, but by changing "plotType" you can use for example Gauss
  plotType = "HPD";
  //plotType = "gaus"; 
    
  MACH3LOG_INFO("Plotting {} errors", plotType);
  
  ReadSettings(man->input().GetFile(0).file);

  canv = new TCanvas("canv", "canv", 1024, 1024);
  //gStyle->SetPalette(51);
  gStyle->SetOptStat(0); //Set 0 to disable statystic box
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.12);
  canv->SetTopMargin(0.08);
  canv->SetRightMargin(0.04);

  canv->Print((SaveName+"[").c_str());

  // these for non named params where we dont need as much space
  p1 = new TPad("p1", "p1", 0.0, 0.3, 1.0, 1.0);
  p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 0.3);
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

  // these for named parameters where we need a nice big gap at the botto to fit the names
  p3 = new TPad("p3", "p3", 0.0, 0.4, 1.0, 1.0);
  p4 = new TPad("p4", "p4", 0.0, 0.0, 1.0, 0.4);
  p3->SetLeftMargin(canv->GetLeftMargin());
  p3->SetRightMargin(canv->GetRightMargin());
  p3->SetTopMargin(canv->GetTopMargin());
  p3->SetBottomMargin(0);
  p3->SetGrid();
  p4->SetLeftMargin(canv->GetLeftMargin());
  p4->SetRightMargin(canv->GetRightMargin());
  p4->SetTopMargin(0);
  p4->SetBottomMargin(0.75);
  p4->SetGrid();

  // Make a Legend page
  TLegend *leg = new TLegend(0.0, 0.0, 1.0, 1.0);

  // make a dummy TH1 to set out legend
  Prefit = new TH1D();
  man->style().setTH1Style(Prefit, man->getOption<std::string>("prefitHistStyle"));
  leg->AddEntry(Prefit, "Prior", "lpf");

  for(int fileId = 0; fileId < man->getNFiles(); fileId++){
    TH1D *postFitHist_tmp = new TH1D();
    postFitHist_tmp->SetBit(kCanDelete);
    
    postFitHist_tmp->SetMarkerColor(TColor::GetColorPalette(fileId));
    postFitHist_tmp->SetLineColor(TColor::GetColorPalette(fileId));
    postFitHist_tmp->SetMarkerStyle(7);
    postFitHist_tmp->SetLineStyle(1+fileId);
    postFitHist_tmp->SetLineWidth(man->getOption<int>("plotLineWidth"));
    leg->AddEntry(postFitHist_tmp, man->getFileLabel(fileId).c_str(), "lpf");
  }
  canv->cd();
  canv->Clear();
  leg->Draw();
  canv->Print((SaveName).c_str());
  delete leg;

  MakeXsecPlots();

  MakeFluxPlots();

  //KS: By default we don't run ProcessMCMC with PlotDet as this take some time, in case we did let's make fancy plots
  MACH3LOG_INFO("ND detector parameters: {}", NDParameters);
  if(NDParameters > 0) MakeNDDetPlots();
  
  //KS: Same as above but for FD parameters,
  if(FDParameters > 0) MakeFDDetPlots();
  
  if(OscParameters > 0) MakeOscPlots();

  canv->Print((SaveName+"]").c_str());

  MakeXsecRidgePlots();

  delete canv;
  delete Prefit;
}


inline TGraphAsymmErrors* MakeTGraphAsymmErrors(std::shared_ptr<TFile> File)
{
  double* x =   new double[nBins];
  double* y =   new double[nBins];
  double* exl = new double[nBins];
  double* eyl = new double[nBins];
  double* exh = new double[nBins];
  double* eyh = new double[nBins];

  TH1D* PostHist = (TH1D*)File->Get( ("param_xsec_"+plotType).c_str() )->Clone();

  TVectorD* Errors_HPD_Positive = (TVectorD*)File->Get( "Errors_HPD_Positive" )->Clone();
  TVectorD* Errors_HPD_Negative = (TVectorD*)File->Get( "Errors_HPD_Negative" )->Clone();

  //KS: I am tempted to multithred this...
  for(int i = 0; i < nBins; ++i)
  {
    //KS: We are extrcting value from three object each having different numbering scheme, I have checked carefully so this is correct please don't cahnge all these +1 +0.5 etc. it just work...
    x[i] = i + 0.5;
    y[i] = PostHist->GetBinContent(i+1);

    //KS: We don't want x axis errors as they are confusing in Violin plot
    exh[i] = 0.00001;
    exl[i] = 0.00001;
    eyh[i] = (*Errors_HPD_Positive)(i);
    eyl[i] = (*Errors_HPD_Negative)(i);
  }
  TGraphAsymmErrors* PostGraph = new TGraphAsymmErrors(nBins,x,y,exl,exh,eyl,eyh);
  PostGraph->SetTitle("");

  delete PostHist;
  delete Errors_HPD_Positive;
  delete Errors_HPD_Negative;
  delete[] x;
  delete[] y;
  delete[] exl;
  delete[] eyl;
  delete[] exh;
  delete[] eyh;

  return PostGraph;
}
   
//KS: Make fancy violin plots
void GetViolinPlots(std::string FileName1 = "", std::string FileName2 = "")
{
  //KS: Should be in some config... either way it control whether you plot symetric or assymetric error bars
  bool PlotAssym = true;
    
  //KS: No idea why but ROOT changed treatment of viilin in R6. If you have non uniform binning this will results in very hard to see violin plots.
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
  std::shared_ptr<TFile> File2 = NULL;
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

  Violin = NULL;
  ViolinPre = (TH2D*)File1->Get( "param_violin_prior" );
  Violin = (TH2D*)File1->Get( "param_violin" );
  if(Violin == NULL)
  {
    MACH3LOG_ERROR("Couldn't find violin plot, make sure method from MCMCProcessor is being called");
    return;
  }

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
    
  TH1D* Postfit = (TH1D*)File1->Get( ("param_xsec_"+plotType).c_str() );
  Postfit->SetMarkerColor(kRed);
  Postfit->SetLineColor(kRed);
  Postfit->SetMarkerStyle(7);
  
  TGraphAsymmErrors* PostGraph = MakeTGraphAsymmErrors(File1);
  PostGraph->SetMarkerColor(kBlack);
  PostGraph->SetLineColor(kBlack);
  PostGraph->SetMarkerStyle(7);
  PostGraph->SetLineWidth(2);
  PostGraph->SetLineStyle(kSolid);
 if(File2 != NULL)
  {
    Violin2 = (TH2D*)File2->Get( "param_violin" );
    Violin2->SetMarkerColor(kGreen);
    Violin2->SetLineColor(kGreen);
    Violin2->SetFillColor(kGreen);
    Violin2->SetFillColorAlpha(kGreen, 0.35);
  }
  
  // Make a Legend page
  TLegend *leg = new TLegend(0.0, 0.0, 1.0, 1.0);
  if (ViolinPre != NULL) leg->AddEntry(ViolinPre, "Prior", "lpf");
  if (Violin != NULL) leg->AddEntry(Violin, "Posterior", "lpf");
  if (Violin2 != NULL) leg->AddEntry(Violin2, "Second Violin", "lpf");
  if(PlotAssym) leg->AddEntry(PostGraph, "HPD Assym", "lp");
  else leg->AddEntry(Postfit, "HPD", "lpf");

  canv->cd();
  canv->Clear();
  leg->Draw();
  canv->Print((SaveName+".pdf").c_str());
  delete leg;
  
   // Do some fancy replacements
  PrettifyTitles(ViolinPre);
  //WARNING this is super hardcoded keep this in mind, ranges are diffefernt from one defined for postfit param
  const int XsecPlots = 8;
  int XsecOffset[XsecPlots] = {XsecStartingPos, XsecStartingPos+16, XsecStartingPos+22, XsecStartingPos+22+12, XsecStartingPos+22+12+11, XsecStartingPos+22+12+11+7,  XsecStartingPos+22+12+7+7+22, XsecStartingPos+CrossSectionParameters-FluxParameters};
  std::string titles[XsecPlots-1] = {"CCQE", "Pauli Blocking and Optical Potential", "2p2h", "SPP", "FSI", "CC DIS, CC Multi #pi, CC coh., NC, #nu_{e}", "CCQE Binding Energy"};

  for (int i = 1; i < XsecPlots; ++i) 
  {
    canv->cd();
    ViolinPre->SetTitle(titles[i-1].c_str());
    ViolinPre->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
    ViolinPre->GetYaxis()->SetRangeUser(-1.5, 2.5);

    Violin->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
    Violin->GetYaxis()->SetRangeUser(-1.5, 2.5);
    if(File2 != NULL) Violin2->GetYaxis()->SetRangeUser(-1.5, 2.5);
    
    std::string name = ViolinPre->GetTitle();
    if (name.find("SPP") != std::string::npos)
    { 
      ViolinPre->GetYaxis()->SetRangeUser(-5., 25.); //For RES Eb we need huge range
      Violin->GetYaxis()->SetRangeUser(-5., 25.);
      if(File2 != NULL) Violin2->GetYaxis()->SetRangeUser(-5., 25.);
    }
    else if (name.find("CCQE Binding Energy") != std::string::npos)   
    {
      ViolinPre->GetYaxis()->SetRangeUser(-10., 16.); //For Eb we need bigger range
      Violin->GetYaxis()->SetRangeUser(-10., 16.); //For Eb we need bigger range
      if(File2 != NULL) Violin2->GetYaxis()->SetRangeUser(-10., 16.);
    }
    else if (name.find("FSI") != std::string::npos)   
    {
      ViolinPre->GetYaxis()->SetRangeUser(-0.5, 2.8);
      Violin->GetYaxis()->SetRangeUser(-0.5, 2.8);
      if(File2 != NULL) Violin2->GetYaxis()->SetRangeUser(-0.5, 2.8);
    }
    else if (name.find("CCQE") != std::string::npos)   
    {
      ViolinPre->GetYaxis()->SetRangeUser(-2., 3.);
      Violin->GetYaxis()->SetRangeUser(-2., 3.);
      if(File2 != NULL) Violin2->GetYaxis()->SetRangeUser(-2., 3.);
    }
    //KS: ROOT6 has some additional options, consider updating it. more https://root.cern/doc/master/classTHistPainter.html#HP140b
    ViolinPre->Draw("violinX(03100300)");
    Violin->Draw("violinX(03100300) SAME");
    if(File2 != NULL)
    {
      Violin2->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      Violin2->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      Violin2->Draw("violinX(03100300) SAME");
    }
    if (Postfit != NULL) {
      Postfit->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      if(!PlotAssym) Postfit->Draw("SAME");
    }
    if (PostGraph != NULL) {
      PostGraph->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      if(PlotAssym) PostGraph->Draw("P SAME");
    }
    
    canv->Print((SaveName+".pdf").c_str());
  }
  
  //KS: Now flux
  if(FluxParameters > 0)
  {
    canv->SetBottomMargin(0.1);
    gStyle->SetOptTitle(0);
    const int FluxInterval = 25;
    for (int nFlux = XsecStartingPos+CrossSectionParameters-FluxParameters; nFlux < XsecStartingPos+CrossSectionParameters; nFlux += FluxInterval)
    {
      ViolinPre->GetXaxis()->SetRangeUser(nFlux, nFlux+FluxInterval);
      Violin->GetXaxis()->SetRangeUser(nFlux, nFlux+FluxInterval);
      Postfit->GetXaxis()->SetRangeUser(nFlux, nFlux+FluxInterval);
      PostGraph->GetXaxis()->SetRangeUser(nFlux, nFlux+FluxInterval);

      ViolinPre->GetYaxis()->SetRangeUser(0.7, 1.3);
      Violin->GetYaxis()->SetRangeUser(0.7, 1.3);
      Postfit->GetYaxis()->SetRangeUser(0.7, 1.3);
      PostGraph->GetYaxis()->SetRangeUser(0.7, 1.3);

      //KS: ROOT6 has some additional options, consider updaiting it. more https://root.cern/doc/master/classTHistPainter.html#HP140b
      ViolinPre->Draw("violinX(03100300)");
      Violin->Draw("violinX(03100300) SAME");
      if(File2 != NULL)
      {
        Violin2->GetYaxis()->SetRangeUser(0.7, 1.3);
        Violin2->GetXaxis()->SetRangeUser(nFlux, nFlux+FluxInterval);
        Violin2->Draw("violinX(03100300) SAME");
      }
      if(PlotAssym) PostGraph->Draw("P SAME");
      else Postfit->Draw("SAME");
      canv->Print((SaveName+".pdf").c_str());
    }
  }
  //KS all other parmeters just in case
  ViolinPre->GetXaxis()->SetRangeUser(CrossSectionParameters, nBins);
  Violin->GetXaxis()->SetRangeUser(CrossSectionParameters, nBins);
  Postfit->GetXaxis()->SetRangeUser(CrossSectionParameters, nBins);
  PostGraph->GetXaxis()->SetRangeUser(CrossSectionParameters, nBins);

  ViolinPre->GetYaxis()->SetRangeUser(-3.4, 3.4);
  Violin->GetYaxis()->SetRangeUser(-3.4, 3.4);
  Postfit->GetYaxis()->SetRangeUser(-3.4, 3.4);
  PostGraph->GetYaxis()->SetRangeUser(-3.4, 3.4);
         
  //KS: ROOT6 has some additional options, consider updating it. more https://root.cern/doc/master/classTHistPainter.html#HP140b
  ViolinPre->Draw("violinX(03100300)");
  Violin->Draw("violinX(03100300) SAME");
  if(File2 != NULL)
  {
    Violin2->GetXaxis()->SetRangeUser(CrossSectionParameters, nBins);
    Violin2->GetYaxis()->SetRangeUser(-3.4, 3.4);
    Violin2->Draw("violinX(03100300) SAME");
  }
  if(PlotAssym) PostGraph->Draw("P SAME");
  else Postfit->Draw("SAME");
  canv->Print((SaveName+".pdf").c_str());
     
  canv->Print((SaveName+".pdf]").c_str());
  delete canv;
  delete ViolinPre;
  delete Violin;
  if(Violin2 != NULL) delete Violin2;
  delete Postfit;
  delete PostGraph;
  File1->Close();
  if(File2 != NULL)
  {
    File2->Close();
  }
}

int main(int argc, char *argv[]) 
{
    SetMaCh3LoggerFormat();
    
    man = new MaCh3Plotting::PlottingManager();
    man->parseInputs(argc, argv);
    std::cout << std::endl << std::endl << "====================" << std::endl;
    man->input().GetFile(0).file->ls();
    
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
        //KS: Violin plot currently not supported by three file veriosn although it should be super easy to adapt
    }
    
  return 0;
}

