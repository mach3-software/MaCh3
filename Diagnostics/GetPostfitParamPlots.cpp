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
#include "TGraphAsymmErrors.h"

#include "manager/manager.h"

//How to use ./GetPostfitParamPlots ProcessMCMC_Output1.root <ProcessMCMC_Output2.root> <ProcessMCMC_Output3.root>
//Originally written by Clarence, with some changes by Will, updated by Kamil, made generic plotter by Ewan
//Central postfit value taken is the Highest Posterior Density but can be changed easily to Gaussian etc. Watch out for parameter names and number of parameters per plot being quite hardcoded

//g++ `root-config --cflags` -g -std=c++11 -o GetPostfitParamPlots GetPostfitParamPlots.cpp -I`root-config --incdir` `root-config --glibs --libs`

//Config
YAML::Node config;

std::vector<std::string> FileNames;
std::vector<std::string> FileLabel = {"Postfit", "Second File", "Third File"};
std::vector <TH1D *> PostfitHistVec;

TH1D *Prefit;

TH2D* ViolinPre;
TH2D* Violin;
TH2D* Violin2;

int CrossSectionParameters;
int FluxParameters;
int XsecStartingPos;

int NDParameters;
int NDParametersStartingPos;

int FDParameters;
int FDParametersStartingPos;

int OscParameters;
int OscParametersStartingPos;

std::vector< int > NDSamplesBins;
std::vector< std::string > NDSamplesNames;
int nBins;
TCanvas *canv;
  
std::string SaveName;
TPad *p1;
TPad *p2;

TPad *p3;
TPad *p4;

//KS: Color for 0 - prefit, 1 postfit, 2 another postfit, 3 you know the drill
Color_t PlotColor[] = {kRed, kBlack, kBlue, kGreen};
std::string plotType;


TH1D *getOtherFileHist(TFile *otherFile, TH1D *refHist){
  // EM: get the post fit errors from a comparison file and put them in a TH1D of the same format as the MaCh3 one
  TH1D *retHist = (TH1D *)refHist->Clone();
  retHist->Reset();

  retHist = (TH1D *)(otherFile->Get( ("param_xsec_"+plotType).c_str() ))->Clone();

  return retHist;
}

void copyParToBlockHist(int localBin, std::string paramName, TH1D*blockHist, TH1D*fullHist, bool setLabels = true){
  // EM: copy parameter from the large TH1 to another hist, finding it in the larger one by name
  int globalBin = 0;
  bool Found = false;
  for (; globalBin <= fullHist->GetNbinsX(); globalBin++) {
    if(std::string(fullHist->GetXaxis()->GetBinLabel(globalBin)) == paramName)
    {
      Found = true;
      break;
    }
  }

  if(!Found){
    std::cerr << "ERROR: Could not find parameter " << paramName << ", given as part of the group " << blockHist->GetTitle() << std::endl;
    throw;
  }

  // Set the values in the sub-histograms
  blockHist->SetBinContent(localBin +1, fullHist->GetBinContent(globalBin));
  blockHist->SetBinError(localBin +1, fullHist->GetBinError(globalBin));

  if(setLabels){
    blockHist->GetXaxis()->SetBinLabel(localBin +1, paramName.c_str());
    blockHist->GetXaxis()->LabelsOption("v");
  }
}

void setTH1Style(TH1* hist, std::string styleName){
  (void) styleName;

  hist->SetMarkerColor(632);
  hist->SetMarkerStyle(7);
  hist->SetFillColor(632);
  hist->SetFillStyle(3003);
  hist->SetLineColor(632);
  hist->SetLineStyle(0);

}

inline TH1D* makeRatio(TH1D *PrefitCopy, TH1D *PostfitCopy, bool setAxes){
  // EM: set up the ratio hist
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

    Ratio->SetLineWidth(2);
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

inline void DrawPlots(TCanvas *canvas, TH1D* PrefitCopy, std::vector<TH1D *>PostfitVec, TPad *mainPad, TPad *ratioPad) {
  // EM: Draw!
  canvas->cd();
  mainPad->Draw();
  mainPad->cd();
  PrefitCopy->GetYaxis()->SetTitle("Parameter Value");

  PrefitCopy->GetYaxis()->SetLabelSize(0.);
  PrefitCopy->GetYaxis()->SetTitleSize(0.05);
  PrefitCopy->GetYaxis()->SetTitleOffset(1.3);
  PrefitCopy->Draw("e2");

  for(int fileId=0; fileId < (int)PostfitVec.size(); fileId++){
    TH1D *postFitHist = PostfitVec[fileId];

    postFitHist->SetMarkerColor(PlotColor[fileId+1]);
    postFitHist->SetLineColor(PlotColor[fileId+1]);
    postFitHist->SetMarkerStyle(7);
    postFitHist->SetLineStyle(1+fileId);
    postFitHist->SetLineWidth(2);

    postFitHist->Draw("e1, same");
  }

  canvas->Update();
  TGaxis *axis = new TGaxis(PrefitCopy->GetXaxis()->GetBinLowEdge(PrefitCopy->GetXaxis()->GetFirst()), gPad->GetUymin()+0.01,
                            PrefitCopy->GetXaxis()->GetBinLowEdge(PrefitCopy->GetXaxis()->GetFirst()), gPad->GetUymax(),
                            gPad->GetUymin()+0.01, gPad->GetUymax(), 510, "");
  axis->SetLabelFont(43);
  axis->SetLabelSize(25);
  axis->Draw();

  canvas->cd();
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
    ratioHists[postFitIdx]->SetLineWidth(2);

    ratioHists[postFitIdx]->Draw("p same");
  }

  // draw lines across the plot at +-1 and 0
  TLine line(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), 0.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), 0.0);
  TLine line2(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), 1.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), 1.0);
  TLine line3(ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetFirst()), -1.0, ratioHists[0]->GetXaxis()->GetBinLowEdge(ratioHists[0]->GetXaxis()->GetLast()+1), -1.0);

  line.SetLineColor(kRed);
  line.SetLineStyle(kDashed);
  line.SetLineWidth(2);
  line2.SetLineColor(kRed);
  line2.SetLineStyle(kDashed);
  line2.SetLineWidth(2);
  line3.SetLineColor(kRed);
  line3.SetLineStyle(kDashed);
  line3.SetLineWidth(2);

  line.Draw("same");
  line2.Draw("same");
  line3.Draw("same");

  canvas->Print((SaveName).c_str());

  ratioHists.clear();
  delete axis;
}


inline std::string FancyTitles(std::string origName) {
  std::string prettyName = origName;

  auto const &PrettyNames = config["GetPostfitParamPlots"]["PrettyNames"];

  if(PrettyNames[prettyName]){
    prettyName = PrettyNames[origName].as<std::string>();
  }

  return prettyName;
}


void PrettifyTitles(TH1D *Hist) {
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i)
  {
    std::string title = Hist->GetXaxis()->GetBinLabel(i+1);
    title = FancyTitles(title);
    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
  }
}

void PrettifyTitles(TH2D *Hist) {
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i)
  {
    std::string title = Hist->GetXaxis()->GetBinLabel(i+1);

    title = FancyTitles(title);
    Hist->GetXaxis()->SetBinLabel(i+1, title.c_str());
  }
}

void ReadSettings(TFile *File1)
{
  TTree *Settings = (TTree*)(File1->Get("Settings"));

  Settings->SetBranchAddress("CrossSectionParameters", &CrossSectionParameters);
  Settings->SetBranchAddress("FluxParameters", &FluxParameters);
  Settings->SetBranchAddress("CrossSectionParametersStartingPos", &XsecStartingPos);
  Settings->SetBranchAddress("NDParameters", &NDParameters);
  Settings->SetBranchAddress("NDParametersStartingPos", &NDParametersStartingPos);
  Settings->SetBranchAddress("FDParameters", &FDParameters);
  Settings->SetBranchAddress("FDParametersStartingPos", &FDParametersStartingPos);
  Settings->SetBranchAddress("OscParameters", &OscParameters);
  Settings->SetBranchAddress("OscParametersStartingPos", &OscParametersStartingPos);

  std::vector< int > *NDSamples_Bins =0;
  std::vector< std::string > *NDSamples_Names = 0;
  Settings->SetBranchAddress("NDSamplesNames", &NDSamples_Names);
  Settings->SetBranchAddress("NDSamplesBins", &NDSamples_Bins);

  Settings->GetEntry(0);

  NDSamplesNames = *NDSamples_Names;
  NDSamplesBins = *NDSamples_Bins;
}

void MakeXsecPlots()
{
  // EM: get the names of the blocks of parameters to group together
  std::vector<std::string> const blockNames = config["GetPostfitParamPlots"]["paramGroups"].as<std::vector<std::string>>();
  const int XsecPlots = (int)blockNames.size();

  for (int i = 0; i < XsecPlots; i++)
  {
    // EM: get the configuration for this parameter
    std::string blockName = blockNames[i];
    YAML::Node paramBlock = config["GetPostfitParamPlots"][blockName];

    std::string blockTitle                 = paramBlock[0].as<std::string>();
    std::vector<double>      blockLimits   = paramBlock[1].as<std::vector<double>>();
    std::vector<std::string> blockContents = paramBlock[2].as<std::vector<std::string>>();

    // EM: get num of params in the block
    int nParams = (int)blockContents.size();

    // EM: set some plot things
    TH1D *blockHist_prefit = new TH1D(blockName.c_str(), blockTitle.c_str(), nParams, 0.0, (double)nParams);
    blockHist_prefit->GetYaxis()->SetTitleOffset(Prefit->GetYaxis()->GetTitleOffset()*1.2);
    setTH1Style(blockHist_prefit, "redHatchedError");

    // EM: set the errors for the prefit block hist
    for(int localBin=0; localBin < nParams; localBin ++){
      // the "local" bin is the params index within the group of parameters
      std::string paramName = blockContents[localBin];
      copyParToBlockHist(localBin, paramName, blockHist_prefit, Prefit);
    }

    // EM: now set for the postfit blocks for all files
    std::vector <TH1D *> blockHist_postfit_Vec;
    for(unsigned int fileId = 0; fileId < FileNames.size(); fileId++){
      TH1D *blockHist_postfit = new TH1D((blockName + FileNames[fileId]).c_str(), blockTitle.c_str(), nParams, 0.0, (double)nParams);

      // EM: loop through all the parameters in this block and set the contents in the blocks TH1
      for(int localBin=0; localBin < nParams; localBin ++){
        // the "local" bin is the params index within the group of parameters
        std::string paramName = blockContents[localBin];
        copyParToBlockHist(localBin, paramName, blockHist_postfit, PostfitHistVec[fileId]);
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
  std::vector<std::string> const fluxBlockNames = config["GetPostfitParamPlots"]["FluxGroups"].as<std::vector<std::string>>();

  auto const fluxBinningTable = config["GetPostfitParamPlots"]["FluxBinning"];

  const int FluxPlots = (int)fluxBlockNames.size();

  for (int i = 0; i < FluxPlots; i++)
  {
    std::string fluxBlockName = fluxBlockNames[i];

    // Access the block details
    YAML::Node paramBlock = config["GetPostfitParamPlots"][fluxBlockName];

    std::string blockTitle = paramBlock[0].as<std::string>();
    std::vector<double> blockLimits = paramBlock[1].as<std::vector<double>>();
    std::string blockBinningName = paramBlock[2].as<std::string>();
    std::vector<int> blockContents = paramBlock[3].as<std::vector<int>>();

    // EM: get the binning for this block of flux params
    YAML::Node binningTable = config["YourBinningTableNode"]; // Assuming you have a separate node for binning
    std::vector<double> binning = fluxBinningTable[blockBinningName].as<std::vector<double>>();

    // now make an array cus root hates vectors
    int nBinEdges = (int)binning.size();
    std::vector<double> binArray(nBinEdges);
    for(int edge = 0; edge < nBinEdges; edge ++) binArray[edge] = binning[edge];

    // get num of params in the block
    int nParams = blockContents[1] - blockContents[0] +1;
    // check for sanity
    if(nParams <= 0 || blockContents.size() > 2){
      std::cerr << "ERROR: Invalid flux parameter block endpoints specified for " << fluxBlockName << std::endl;
      std::cerr << "       Should have the form [<low index>, <up index>]" << std::endl;
      throw;
    }
    if(nParams != (int)binning.size() -1){
      std::cerr << "ERROR: Binning provided for flux param block "  << fluxBlockName << " Does not match the number of parameters specified for the block" << std::endl;
      std::cerr << "       Provided " << nParams << " flux parameters, but " << binning.size() -1 << " bins" << std::endl;
      throw;
    }

    TH1D *blockHist_prefit = new TH1D(fluxBlockName.c_str(), blockTitle.c_str(), nBinEdges - 1, binArray.data());
    blockHist_prefit->GetYaxis()->SetTitle("Parameter Variation");
    blockHist_prefit->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    blockHist_prefit->GetXaxis()->SetTitleOffset(blockHist_prefit->GetXaxis()->GetTitleOffset()*1.2);
    setTH1Style(blockHist_prefit, "redHatchedError");
    // set the errors for the prefit block hist
    for(int fluxParId = blockContents[0]; fluxParId <= blockContents[1]; fluxParId++){
      int localBin = fluxParId - blockContents[0];
      std::string paramName = "b_" + std::to_string(fluxParId);
      copyParToBlockHist(localBin, paramName, blockHist_prefit, Prefit, false);
    }

    // now set for the postfit blocks for all files
    std::vector <TH1D *> blockHist_postfit_Vec;
    for(unsigned int fileId = 0; fileId < FileNames.size(); fileId++){
      TH1D *blockHist_postfit = new TH1D(fluxBlockName.c_str(), blockTitle.c_str(), nBinEdges - 1, binArray.data());

      for(int fluxParId = blockContents[0]; fluxParId <= blockContents[1]; fluxParId++){
        int localBin = fluxParId - blockContents[0];
        std::string paramName = "b_" + std::to_string(fluxParId);

        copyParToBlockHist(localBin, paramName, blockHist_postfit, PostfitHistVec[fileId], false);
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
  for (unsigned int i = 0; i < NDSamplesNames.size(); ++i)
  {
    NDbinCounter += NDSamplesBins[i];

    std::string temp = NDSamplesNames[i].c_str();
    while (temp.find("_") != std::string::npos) {
      temp.replace(temp.find("_"), 1, std::string(" "));
    }
    Prefit->SetTitle(temp.c_str());
    Prefit->GetXaxis()->SetRangeUser(Start, NDbinCounter);

    // set the x range for the postfits
    for(unsigned int fileId = 0; fileId < FileNames.size(); fileId++){
      PostfitHistVec[fileId]->GetXaxis()->SetRangeUser(Start, NDbinCounter);
    }

    //KS: We dont' need name for every nd param
    for(int j = 0; j < NDSamplesBins[i]; ++j)
    {
      bool ProductOfTen = false;
      if(j % 10) ProductOfTen = true;
      if(j != 0 && ProductOfTen) Prefit->GetXaxis()->SetBinLabel(Start+j+1, " ");
    }

    Prefit->GetYaxis()->SetRangeUser(0.5, 1.5);

    Start += NDSamplesBins[i];

    DrawPlots(canv, Prefit, PostfitHistVec, p3, p4);
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
    //KS: We don't' need name for every FD param
    for(int j = 0; j < FDSamplesBins[i]; ++j)
    {
      //bool ProductOfTen = false;
      //if(j % 10) ProductOfTen = true;
      //if(j != 0 && ProductOfTen) Prefit->GetXaxis()->SetBinLabel(Start+j+1, " ");
    }
    Prefit->GetYaxis()->SetRangeUser(0.5, 1.5);
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

  Prefit->GetYaxis()->SetRangeUser(-3, +3);
  Prefit->Draw("e2");

  DrawPlots(canv, Prefit, PostfitHistVec, p3, p4);
}


void GetPostfitParamPlots()
{
  std::cout<<"Input Files: "<<std::endl;
  for(unsigned int i = 0; i < FileNames.size(); i++)
  {
    SaveName += FileNames[i] + "_";
    std::cout<< "  " << FileNames[i] << std::endl;
  }
  SaveName += "_PostfitParamPlots.pdf";

  //KS: By default we take HPD values, but by changing "plotType" you can use for example Gauss
  plotType = "HPD";
  //plotType = "gaus";

  std::cout<<"Plotting "<<plotType<<std::endl;

  // open the files
  std::vector <TFile*> fileVec;
  for(unsigned int i = 0; i < FileNames.size(); i++)
  {
    TFile *file = new TFile(FileNames[i].c_str());
    fileVec.push_back(file);
  }

  ReadSettings(fileVec[0]);

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
  p4->SetBottomMargin(0.62);
  p4->SetGrid();

  Prefit = (TH1D*)fileVec[0]->Get("param_xsec_prefit");

  for(unsigned int fileId = 0; fileId < FileNames.size(); fileId ++)
  {
    TH1D *Hist = getOtherFileHist(fileVec[fileId], Prefit);
    Hist->SetMarkerColor(PlotColor[fileId+1]);
    Hist->SetLineColor(PlotColor[fileId+1]);
    Hist->SetMarkerStyle(7);
    Hist->SetLineStyle(1+fileId);

    PostfitHistVec.push_back(Hist);
  }

  nBins = Prefit->GetXaxis()->GetNbins();

  Prefit->SetFillColor(PlotColor[0]);
  Prefit->SetFillStyle(3003);
  Prefit->SetMarkerStyle(7);
  Prefit->SetMarkerColor(PlotColor[0]);
  Prefit->SetLineColor(PlotColor[0]);
  Prefit->SetLineStyle(kSolid);

  // Make a Legend page
  TLegend *leg = new TLegend(0.0, 0.0, 1.0, 1.0);
  if (Prefit != NULL) leg->AddEntry(Prefit, "Prior", "lpf");
  for(unsigned int fileId = 0; fileId < FileNames.size(); fileId++){
    leg->AddEntry(PostfitHistVec[fileId], FileLabel[fileId].c_str(), "lpf");
  }
  canv->cd();
  canv->Clear();
  leg->Draw();
  canv->Print((SaveName).c_str());
  delete leg;

  MakeXsecPlots();

  MakeFluxPlots();

  //KS: By default we don't run ProcessMCMC with PlotDet as this take some time, in case we did let's make fancy plots
  if(NDParameters > 0) MakeNDDetPlots();

  //KS: Same as above but for FD parameters,
  if(FDParameters > 0) MakeFDDetPlots();

  if(OscParameters > 0) MakeOscPlots();

  canv->Print((SaveName+"]").c_str());

  delete canv;
  delete Prefit;
  for(TFile *file: fileVec) file->Close();
  fileVec.clear();
}


inline TGraphAsymmErrors* MakeTGraphAsymmErrors(TFile *File)
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
    
  std::cout<<" Making Violin Plot"<<std::endl;
  if (!FileName1.empty()) std::cout << "File 1 " << FileName1<< std::endl;
  if (!FileName2.empty()) std::cout << "File 2 " << FileName2<< std::endl;
    
  SaveName = FileName1;
  SaveName = SaveName.substr(0, SaveName.find(".root"));
  if(FileName2 != "") SaveName += FileName2;
  if(FileName2 != "") SaveName = SaveName.substr(0, SaveName.find(".root"));
  SaveName += "_Violin";
  if(PlotAssym) SaveName += "_Assym";
   
  TFile *File1 = new TFile(FileName1.c_str());
  TFile *File2 = NULL;
  if(FileName2 != "") File2 = new TFile(FileName2.c_str());
  
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
    std::cout<<"Couldn't find violin plot, make sure method from MCMCProcessor is being called"<<std::endl;
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
    //KS: ROOT6 has some additional options, consider updaiting it. more https://root.cern/doc/master/classTHistPainter.html#HP140b
    ViolinPre->Draw("VIOLIN");
    Violin->Draw("VIOLIN SAME");
    if(File2 != NULL)
    {
      Violin2->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      Violin2->GetXaxis()->SetRangeUser(XsecOffset[i-1], XsecOffset[i]);
      Violin2->Draw("VIOLIN SAME");
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
      ViolinPre->Draw("VIOLIN");
      Violin->Draw("VIOLIN SAME");
      if(File2 != NULL)
      {
        Violin2->GetYaxis()->SetRangeUser(0.7, 1.3);
        Violin2->GetXaxis()->SetRangeUser(nFlux, nFlux+FluxInterval);
        Violin2->Draw("VIOLIN SAME");
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
         
  //KS: ROOT6 has some additional options, consider updaiting it. more https://root.cern/doc/master/classTHistPainter.html#HP140b
  ViolinPre->Draw("VIOLIN");
  Violin->Draw("VIOLIN SAME");
  if(File2 != NULL)
  {
    Violin2->GetXaxis()->SetRangeUser(CrossSectionParameters, nBins);
    Violin2->GetYaxis()->SetRangeUser(-3.4, 3.4);
    Violin2->Draw("VIOLIN SAME");
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
  delete File1;
  if(File2 != NULL)
  {
    File2->Close();
    delete File2;
  }
}

int main(int argc, char *argv[]) 
{
  SetMaCh3LoggerFormat();

  if (argc != 3 && argc != 4 && argc != 5)
  {
    std::cerr << "How to use: "<< argv[0] << "Config.yaml MCMC_Processor_Output.root" << std::endl;
    std::cerr << "You can add up to 3 different files" << std::endl;
    exit(-1);
  }

  config = YAML::LoadFile(std::string(argv[1]));
  if (argc == 3)
  {
    FileNames.push_back(argv[2]);

    GetPostfitParamPlots();
    GetViolinPlots(argv[2]);
  }
  else if (argc == 4)
  {
    FileNames.push_back(std::string(argv[2]));
    FileNames.push_back(std::string(argv[3]));

    GetPostfitParamPlots();
    GetViolinPlots(std::string(argv[2]), std::string(argv[3]));
  }
  else if (argc == 5)
  {
    FileNames.push_back(argv[2]);
    FileNames.push_back(argv[3]);
    FileNames.push_back(argv[4]);

    GetPostfitParamPlots();
    //KS: Violin plot currently not supported by three file version although it should be super easy to adapt
  }

  return 0;
}

