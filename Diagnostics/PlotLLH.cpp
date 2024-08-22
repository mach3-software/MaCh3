// C++
#include <iostream>
#include <algorithm>
#include <iomanip>

// ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TKey.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLine.h"
#include "TROOT.h"

#include "manager/manager.h"

std::vector<std::string> ExtraFileNames;
std::vector<std::string> FileLabels;
std::vector<std::string> FileLabels_default;
std::string PDFName = "LLHScan";
std::string Mach3LLHFileName;
std::string extraDrawOptions = "";

bool splitBySample;
bool plotRatios;
bool drawGrid;

// some options for the plots
double ratioPlotSplit = 0.3;
double yTitleOffset = 1.25;
double sampleLabelThreshold = 0.0275;
int lineWidth = 3;
bool totalOnSplitPlots = false;
bool sameAxis = true;

std::string NOT_FOUND_STR = "__PARAM_NAME_NOT_FOUND__";

// these directories are not considered when looking for LLH broken down by sample
std::vector<std::string> nonSampleDirNames = {
  "Sample_LLH", "Total_LLH", // General
  "XSec_LLH", "ND280_LLH", // ND Specific
  "Xsec_LLH", "SK_LLH" // SK Specific
  };

//EM: the corresponding pre/ post fixes for those directories are:
//    "sam", "full", // General
//    "xs", "nd", // ND Specific
//    "xs", "sk" // SK Specific

//EM: this is the directory that will be used for the base comparisson, maybe make this a command line option in future
std::string scanDirPrefix = "sam";
std::string scanDirPath = "Sample_LLH";

// scale the ratio plot labels by this much to make them same size as the normal plot
double ratioLabelScaling = (1.0/ratioPlotSplit - 1.0);

void parseFileLabels(std::string labelString, std::vector<std::string> &labelVec){
    int end = labelString.find(";"); 
    while (end != -1) { // Loop until no delimiter is left in the string.
        labelVec.push_back(labelString.substr(0, end));
        labelString.erase(labelString.begin(), labelString.begin() + end + 1);
        end = labelString.find(";");
    }
    labelVec.push_back(labelString.substr(0, end));
}


bool getGUNDAMParamName(std::string mach3ParamName, TDirectoryFile *llhDir, std::string &gundamName){
    // try to find the right branch by looking for the mach3 name in the gundam file
    
    gundamName = NOT_FOUND_STR;
    for (TObject *keyObj: *llhDir->GetListOfKeys()) {
      TKey* key = (TKey*)keyObj;
      std::string testString = key->GetName();
           
      if(testString.find(mach3ParamName + "_TGraph") != std::string::npos){
        gundamName = testString;

        break;
      }
      
      // special case for EB_dial
      else if(mach3ParamName.find("EB_dial") != std::string::npos){
        std::string mach3ParamName_reduced = std::string(mach3ParamName);
        mach3ParamName_reduced.replace(mach3ParamName_reduced.find("EB_dial"), std::string("EB_dial").size(), std::string("EB_bin")); 
        
        if(testString.find(mach3ParamName_reduced + "_TGraph") != std::string::npos){
          gundamName = testString;

          break;
        }
      }

      else {
        // if we didnt find it, try again with potential extra bits of string removed
        std::vector<std::string> mach3OnlyLabels = {"NDS_"};
        for(std::string extraBit: mach3OnlyLabels){
          std::string mach3ParamName_reduced = std::string(mach3ParamName);
          if(mach3ParamName_reduced.find(extraBit) != std::string::npos)
            mach3ParamName_reduced.erase(mach3ParamName_reduced.find(extraBit), extraBit.length());
          else 
            continue;

          if(testString.find(mach3ParamName_reduced + "_TGraph") != std::string::npos){
            gundamName = testString;

            break;
          }
        }
      }
    }

    return !(gundamName == NOT_FOUND_STR); // true if found, false if not
}

bool getParameterHist(TFile *file, TH1D *retHist, std::string mach3ParamName, std::string directory = "sam", std::string directoryPath = "Sample_LLH", bool printNames=false){ // default to sample directory
  // check whats in the file to determine which fitter it cam from 
  TDirectoryFile *sampleLLHDir = (TDirectoryFile*)file->Get("Sample_LLH");
  TDirectoryFile *FitterEngine = (TDirectoryFile*)file->Get("FitterEngine");

  if(sampleLLHDir != NULL){ // this is a mach3 LLH scan
    TH1D *Hist;

    // see if we recognise the directory as a non sample one
    bool recognisedDir = false;
    for(std::string compStr: nonSampleDirNames) if(compStr == directoryPath) recognisedDir=true;
    
    if(recognisedDir){
      TDirectoryFile *LLHDir = (TDirectoryFile*)file->Get(directoryPath.c_str());
      Hist = (TH1D*)LLHDir->Get(Form("%s_%s", mach3ParamName.c_str(), directory.c_str()));

      if(!Hist){ // might be an SK LLH scan
        std::cout<<Form("%s_%s", directory.c_str(), mach3ParamName.c_str())<<std::endl;
        Hist = (TH1D*)LLHDir->Get(Form("%s_%s", directory.c_str(), mach3ParamName.c_str()));
      }
    }

    else{ // assume directory for a specific sample
      TDirectoryFile *dir = (TDirectoryFile*)file->Get(Form("%s_LLH", directory.c_str())); 
      Hist = (TH1D*)dir->Get(Form("%s%s", mach3ParamName.c_str(), directory.c_str()));
    }

    // if we havent found it by now, give up
    if(!Hist) return false;

    *retHist = TH1D(*Hist);
  }

  else if(FitterEngine != NULL){ // this is a gundam LLH scan
    TDirectoryFile *postFit = (TDirectoryFile*)FitterEngine->Get("preFit");
    TDirectoryFile *scan    = (TDirectoryFile*)postFit->Get("scan");

    // get the llh folder from the correct directory
    TDirectoryFile *llh;
    if(directory == "sam") llh = (TDirectoryFile*)scan->Get("llhStat");
    else if (directory == "full") llh = (TDirectoryFile*)scan->Get("llh");
    else if (directory == "xs") llh = (TDirectoryFile*)scan->Get("llhPenalty");
    else{
      std::cerr<<"ERROR: unknown directory option: "<<directory<<std::endl;
      throw;
    }

    std::string gundamName;
    bool found = getGUNDAMParamName(mach3ParamName, llh, gundamName);
    if(printNames)std::cout << std::setw(40) << mach3ParamName << " | "<< std::setw(40) << gundamName << std::endl;

    if(!found) return false; // stop here if we couldnt find an equivalent parameter

    TKey *paramKey = llh->FindKey(gundamName.c_str());
    std::cout<<gundamName<<std::endl;
    
    TGraph *paramGraph = (TGraph *)paramKey->ReadObj();

    // fill the hist from the graph that we loaded
    for(int idx =0; idx <= retHist->GetNbinsX(); idx++){
      retHist->SetBinContent(idx, paramGraph->Eval(retHist->GetBinCenter(idx)));
    }
  }

  else { 
    std::cerr<<"ERROR: Hey, really sorry but I'm not smart enough just now to read from the file you provided: "<<file->GetName()<<std::endl;
    std::cerr<<"       I omly know about MaCh3 and gundam files just now but You can teach me about yours if you like by adding the code to parse it here: "<<std::endl;
    std::cerr<<"       "<<__FILE__<<":"<<__LINE__<<std::endl;
    throw;
  }

  return true;
}


void set_plot_style()
{
    //EM: define the colour palette taken from https://www.color-hex.com/color-palette/49436
    const Int_t NCont = 255;
    
    const Int_t NRGBs = 4;
    Double_t stops[NRGBs] = { 0.0        , 0.25       , 0.5        , 0.75       }; //, 1.00        };
    Double_t red[NRGBs]   = { 213.0/256.0, 204.0/256.0, 0.0        , 240.0/256.0}; //, 0.0         };
    Double_t green[NRGBs] = { 94.0/256.0 , 121.0/256.0, 114.0/256.0, 228.0/256.0}; //, 158.0/256.0 };
    Double_t blue[NRGBs]  = { 0.0        , 167.0/256.0, 178.0/256.0, 66.0/256.0 }; //, 115.0/256.0 };
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}


std::string prettify_name(std::string fullName){
    std::string newName = fullName;
    
    std::string numu = "#nu_{#mu}";
    std::string numuBar = "#bar{#nu}_{#mu}";
    std::string pi = "#pi";
    std::string gamma = "#gamma";
    std::string nuBar = "#bar{#nu}";
    std::string boldP = "#bf{P}";
    
    if(newName.find("anti-numu") <= newName.length())newName.replace(newName.find("anti-numu"), 9, numuBar);
    if(newName.find("numu")      <= newName.length())newName.replace(newName.find("numu"), 4, numu);
    if(newName.find("NuMu")      <= newName.length())newName.replace(newName.find("NuMu"), 4, numu);
    if(newName.find("photon")    <= newName.length())newName.replace(newName.find("photon"), 6, gamma);
    if(newName.find("AntiNu")    <= newName.length())newName.replace(newName.find("AntiNu"), 6, nuBar);
    if(newName.find("protons")   <= newName.length())newName.replace(newName.find("protons"), 7, boldP);
    
    if(newName.find("Mode")      <= newName.length())newName.erase(newName.find("Mode"), 4);
    
    return newName;
}

void getSplitSampleStack(TFile *file, std::vector<std::string> sampleVector, std::string parameterName, TH1D LLH_allSams, std::vector<Double_t> &cumSums, std::vector<bool> &drawLabel, THStack *sampleStack, TLegend *splitSamplesLegend, float baselineLLH_sam = 0.00001){
  float LLH_sam_integ = LLH_allSams.Integral();
  float cumSum = 0.0;
  int nBins = LLH_allSams.GetNbinsX();

  for(uint i=0; i < sampleVector.size(); i++){
    std::string sampName = sampleVector[i];   
    
    TH1D *LLH_indivSam = new TH1D(Form("%s_%s_%s", parameterName.c_str(), sampName.c_str(), file->GetName()), "", 
                                  nBins, LLH_allSams.GetBinLowEdge(1), LLH_allSams.GetBinLowEdge(nBins+1));

    getParameterHist(file, LLH_indivSam, parameterName, sampName, sampName); // NB: This will only work if the file is a mach3 file

    LLH_indivSam->SetStats(0);
    LLH_indivSam->SetLineColor(TColor::GetColorPalette(floor((float)i * TColor::GetNumberOfColors()/ (float)sampleVector.size())));
    LLH_indivSam->SetFillColor(TColor::GetColorPalette(floor((float)i * TColor::GetNumberOfColors()/ (float)sampleVector.size())));
    sampleStack->Add(LLH_indivSam); 
    splitSamplesLegend->AddEntry(LLH_indivSam, prettify_name(sampName).c_str(), "lf");

    cumSum += LLH_indivSam->GetBinContent(nBins);

    cumSums.push_back(cumSum);
    
    if((LLH_indivSam->Integral()/ LLH_sam_integ > sampleLabelThreshold) && (LLH_indivSam->Integral() / baselineLLH_sam > sampleLabelThreshold) ){ //dont draw a label if the likelihood contribution is less than threshold%
      drawLabel.push_back(true);
    }
    else drawLabel.push_back(false);
  }

  return;
}


int PlotLLH(){
    // open the additional files
    std::vector<TFile *> ExtraCompFiles(ExtraFileNames.size());
    for(uint extraFileId = 0; extraFileId < ExtraFileNames.size(); extraFileId ++){
      ExtraCompFiles[extraFileId] = new TFile(ExtraFileNames[extraFileId].c_str());
    }

    //gStyle->SetPalette(1);
    set_plot_style();
    
    TFile* file1 = new TFile(Mach3LLHFileName.c_str(),"READ");
    
    // get the names of the directories
    TList *dirListFull = file1->GetListOfKeys();

    //EM: get the names of the samples from the file
    std::vector<std::string> sampleDirList; // the names of the directories for each sample
    std::vector<std::string> sampleList; // the names of the samples

    // get the names of the parameters from the Sample_LLH directory
    TDirectory *sampleDir = (TDirectory*)file1->Get(scanDirPath.c_str());
    std::vector<std::string> paramList;
    
    std::cout<<"found parameters:"<<std::endl;
    TList *paramKeys = sampleDir->GetListOfKeys();
    for(int i=0; i < paramKeys->GetSize(); i++){
        std::string paramName = ((TKey*)paramKeys->At(i))->GetName();
        
        // pop out the _sam at the end of the name
        std::string removeString = std::string(Form("_%s", scanDirPrefix.c_str()));
        if(paramName.find(removeString) != std::string::npos){
          paramName.erase(paramName.find(removeString), removeString.length());
        }

        // and just in case this is an SK scan
        removeString = std::string(Form("%s_", scanDirPrefix.c_str()));
        if(paramName.find(removeString) != std::string::npos){
          paramName.erase(paramName.find(removeString), removeString.length());
        }

        std::cout<<"\t"<<paramName<<std::endl;
        paramList.push_back(paramName);
    }
    
    if(splitBySample){
        std::cout<<"Found samples: "<<std::endl;
        for(int i=0; i < dirListFull->GetSize(); i++)
        {
            std::string dirName = ((TKey*)dirListFull->At(i))->GetName();
            if (!(std::find(nonSampleDirNames.begin(), nonSampleDirNames.end(), dirName) != nonSampleDirNames.end()))
            {
                std::cout<<"\t"<<dirName<<" : ";
                sampleDirList.push_back(dirName);
                
                // pop out the _LLH at the end of the directory name
                std::string removeString = "_LLH";
                dirName.erase(dirName.find(removeString), dirName.length());
                std::cout<<"\t"<<dirName<<std::endl;
                sampleList.push_back(dirName);
            }        
        }
        
        if(sampleDirList.size()==0){
            std::cout<<"found no sample directories, if you want to split the LLH by sample you neet to run LLH scan with LLH_SCAN_BY_SAMPLE = true option in your config file"<<std::endl;
            return -1;
        }
    }

    TCanvas *blankCanv = new TCanvas("blankCanv","",1024,1024);
    TCanvas *splitSampBlankCanv = new TCanvas("splitSampBlankCanv","",1024 * (1 + ExtraFileNames.size()),1024);
    blankCanv->SaveAs(Form("%s.pdf[", PDFName.c_str()));
    if(splitBySample)splitSampBlankCanv->SaveAs(Form("%s_bySample.pdf[", PDFName.c_str()));

    // loop over the spline parameters
    for(std::string paramName : paramList){
      std::cout<<"working on parameter "<<paramName<<std::endl;
      // ###############################################################
      // First lets do just the straight up likelihoods from all samples
      // ###############################################################
      
      // make the canvas and other plotting stuff
      TCanvas* allSamplesCanv = new TCanvas("AllSampCanv","",1024,1024);
      TLegend* legend = new TLegend(0.3,0.6,0.7,0.8);
      THStack *compStack   = new THStack(Form("%s_%s", paramName.c_str(), scanDirPrefix.c_str()), Form("%s_%s", paramName.c_str(), scanDirPrefix.c_str()));
      THStack *ratioCompStack   = new THStack(Form("%s_sam_ratio", paramName.c_str()), "");

      // split the canvas if plotting ratios
      TPad *AllSamplesPad, *AllSamplesRatioPad;
      if(plotRatios){
        AllSamplesPad = new TPad("AllSampPad", "AllSampPad", 0.0, ratioPlotSplit, 1.0, 1.0);
        AllSamplesPad->SetBottomMargin(0.0);
        AllSamplesRatioPad = new TPad("AllSampRatioPad", "AllSampRatioPad", 0.0, 0.0 , 1.0, ratioPlotSplit);
        AllSamplesRatioPad->SetTopMargin(0.0);
        AllSamplesRatioPad->SetBottomMargin(0.3);
      }
      else{
        AllSamplesPad = new TPad("AllSampPad", "AllSampPad", 0.0, 0.0 ,1.0, 1.0);
        AllSamplesRatioPad = new TPad("AllSampRatioPad", "AllSampRatioPad",0.0 ,0.0 , 0.0, 0.0);
      }

      // get the sample reweight hist from the main file 
      TH1D *LLH_sam;
      LLH_sam = (TH1D*)sampleDir->Get(Form("%s_%s", paramName.c_str(), scanDirPrefix.c_str()));
      
      // maybe SK
      if(!LLH_sam)LLH_sam = (TH1D*)sampleDir->Get(Form("%s_%s", scanDirPrefix.c_str(), paramName.c_str()));
      
      if(!LLH_sam){
        std::cerr<<"ERROR: failed to find the LLH scan hist in "<<scanDirPath<<std::endl;
        throw;
      }


      LLH_sam->SetStats(0);
      LLH_sam->SetLineColor(kBlack);
      compStack->Add(LLH_sam);
      legend->AddEntry(LLH_sam,FileLabels[0].c_str(),"l");
      
      int nBins = LLH_sam->GetNbinsX();
      
      // go through the other files 
      for(uint extraFileIdx = 0; extraFileIdx < ExtraFileNames.size(); extraFileIdx ++ ){
        // first get the corresponding histogram from the comparisson file
        TFile *compFile = ExtraCompFiles[extraFileIdx];
        
        TH1D *compHist = new TH1D(Form("%s_%s", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()), Form("%s_%s", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()), 
                                  nBins, LLH_sam->GetBinLowEdge(1), LLH_sam->GetBinLowEdge(nBins+1));
        TH1D *divHist = new TH1D(Form("%s_%s_div", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()), Form("%s_%s_div", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()), 
                                  nBins, LLH_sam->GetBinLowEdge(1), LLH_sam->GetBinLowEdge(nBins+1));
                           
        if(!getParameterHist(compFile, compHist, paramName, scanDirPrefix, scanDirPath)){
          std::cout<<"WARNING: didnt find equivalent of parameter "<< paramName<<" in file "<<compFile->GetName()<<" wont be included in the comparisson plot"<<std::endl;
          continue;
        }
        
        // make them look different to each other
        compHist->SetLineColor(TColor::GetColorPalette(floor((float)extraFileIdx * TColor::GetNumberOfColors()/ (float)(ExtraFileNames.size() +1))));
        compHist->SetLineStyle(2+ extraFileIdx %9);
        compHist->SetLineWidth(lineWidth);
 
        if(plotRatios){
          // get the ratio hist
          divHist->Divide(compHist, LLH_sam);
          divHist->SetLineColor(TColor::GetColorPalette(floor((float)extraFileIdx * TColor::GetNumberOfColors()/ (float)(ExtraFileNames.size() +1))));
          divHist->SetLineStyle(2+ extraFileIdx %9);
          divHist->SetLineWidth(lineWidth);
          ratioCompStack->Add(divHist);
        }

        // add it to the comparisson hstack and legend
        compStack->Add(compHist); 
        legend->AddEntry(compHist, FileLabels[1 + extraFileIdx].c_str(), "l");
      }

      // draw the log likelihoods
      allSamplesCanv->cd();
      allSamplesCanv->Draw();
      AllSamplesPad->Draw();
      AllSamplesPad->cd();
      if(drawGrid)AllSamplesPad->SetGrid();
      compStack->Draw(Form("NOSTACK%s", extraDrawOptions.c_str()));
      if(!plotRatios)compStack->GetXaxis()->SetTitle("Parameter Variation");
      compStack->GetYaxis()->SetTitle(Form("-2LLH_{%s}", scanDirPrefix.c_str()));
      compStack->GetYaxis()->SetTitleOffset(yTitleOffset);
      legend->Draw();
      
      // add the ratio plot if specified
      if(plotRatios){
        allSamplesCanv->cd();
        AllSamplesRatioPad->Draw();
        AllSamplesRatioPad->cd();
        if(drawGrid)AllSamplesRatioPad->SetGrid();
        
        // do this so 1.0 is in the middle of the plot vertically
        double stackMax, stackMin;
        stackMax = ratioCompStack->GetMaximum("NOSTACK");
        stackMin = ratioCompStack->GetMinimum("NOSTACK");

        double stackLim = std::max(std::abs(1.0 - stackMax), std::abs(1.0 - stackMin));

        ratioCompStack->SetMinimum(1.0 - 1.05 * stackLim);
        ratioCompStack->SetMaximum(1.0 + 1.05 * stackLim);

        // draw it
        ratioCompStack->Draw(Form("NOSTACK%s", extraDrawOptions.c_str()));

        // make it look a bit nicer
        ratioCompStack->GetXaxis()->SetLabelSize(ratioLabelScaling * compStack->GetXaxis()->GetLabelSize());
        ratioCompStack->GetXaxis()->SetTitleSize(ratioLabelScaling * compStack->GetXaxis()->GetLabelSize());
        ratioCompStack->GetXaxis()->SetTitle("Parameter Variation");
        
        ratioCompStack->GetYaxis()->SetLabelSize(ratioLabelScaling * compStack->GetYaxis()->GetLabelSize());
        ratioCompStack->GetYaxis()->SetTitleSize(ratioLabelScaling * compStack->GetYaxis()->GetLabelSize());
        ratioCompStack->GetYaxis()->SetTitleOffset(yTitleOffset);
        ratioCompStack->GetYaxis()->SetTitle("Parameter Variation");
        ratioCompStack->GetYaxis()->SetNdivisions(5,2,0);

        // add horizontal line at 1 for reference
        TLine line = TLine();
        line.SetLineColor(kBlack);
        line.SetLineWidth(lineWidth);
        line.DrawLine(LLH_sam->GetBinLowEdge(1), 1.0, LLH_sam->GetBinLowEdge(nBins+1), 1.0);

      }

      // save to the output file
      allSamplesCanv->SaveAs(Form("%s.pdf", PDFName.c_str()));


      if(splitBySample){
        // #########################################
        // ## now lets make plots split by sample ##
        // #########################################

        // split up the canvas so can have side by side plots, one for each file
        TCanvas* splitSamplesCanv = new TCanvas("splitSampCanv","",4096 * (1 + ExtraFileNames.size()), 4096);
        splitSamplesCanv->Divide(1 + ExtraFileNames.size());
        splitSamplesCanv->Draw();

        THStack *baseSplitSamplesStack = new THStack(paramName.c_str(), Form("%s - %s", paramName.c_str(), FileLabels[0].c_str()));
        TLegend* baseSplitSamplesLegend = new TLegend(0.37,0.475,0.63,0.9);

        splitSamplesCanv->cd(1)->SetGrid();

        std::vector<Double_t> cumSums;
        std::vector<bool> drawLabel;

        getSplitSampleStack(file1, sampleList, paramName, *LLH_sam, cumSums, drawLabel, baseSplitSamplesStack, baseSplitSamplesLegend);

        baseSplitSamplesStack->Draw(extraDrawOptions.c_str());

        if(totalOnSplitPlots){
            LLH_sam->SetLineWidth(1); // undo SetLineWidth that was done above
            LLH_sam->Draw(Form("same%s", extraDrawOptions.c_str()));
            baseSplitSamplesLegend->AddEntry(LLH_sam, "All Samples", "l");
          }
        
        baseSplitSamplesLegend->Draw();

        TLatex *label = new TLatex;
        // format the label
        label->SetTextAlign(11);
        label->SetTextAngle(-55);
        label->SetTextSize(0.012);

        // need to draw the labels after other stuff or they dont show up
        for(uint i=0; i < sampleDirList.size(); i++){
          std::string sampName = sampleList[i];
          if(!drawLabel[i])continue;
          label->DrawLatex(LLH_sam->GetBinLowEdge(nBins+1), cumSums[i], Form("#leftarrow%s", prettify_name(sampName).c_str()) );
        }
        
        
        // now we plot the comparisson file plots
        for(uint extraFileIdx = 0; extraFileIdx < ExtraFileNames.size(); extraFileIdx ++ ){
          splitSamplesCanv->cd(2 + extraFileIdx);

          // split the canvas if plotting ratios
          TPad *splitSamplesPad, *splitSamplesRatioPad;
          if(plotRatios){
            splitSamplesPad = new TPad("splitSampPad", "splitSampPad", 0.0, ratioPlotSplit, 1.0, 1.0);
            splitSamplesPad->SetBottomMargin(0.0);
            splitSamplesRatioPad = new TPad("splitSampRatioPad", "splitSampRatioPad", 0.0, 0.0 , 1.0, ratioPlotSplit);
            splitSamplesRatioPad->SetTopMargin(0.0);
            splitSamplesRatioPad->SetBottomMargin(0.3);
          }
          else{
            splitSamplesPad = new TPad("splitSampPad", "splitSampPad", 0.0, 0.0 ,1.0, 1.0);
            splitSamplesRatioPad = new TPad("splitSampRatioPad", "splitSampRatioPad",0.0 ,0.0 , 0.0, 0.0);
          }

          std::vector<Double_t> cum_Sums;
          std::vector<bool> draw_Label;

          THStack *splitSamplesStack = new THStack(paramName.c_str(), Form("%s - %s", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()));
          TLegend *splitSamplesLegend = new TLegend(0.37,0.475,0.63,0.9);
          
          TFile *compFile = ExtraCompFiles[extraFileIdx];
          TH1D *compLLH_sam = new TH1D(Form("%s_%s", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()), Form("%s_%s", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()), 
                                  nBins, LLH_sam->GetBinLowEdge(1), LLH_sam->GetBinLowEdge(nBins+1));

          if(!getParameterHist(compFile, compLLH_sam, paramName, scanDirPrefix, scanDirPath)){
            std::cout<<"WARNING: didnt find equivalent of parameter "<< paramName<<" in file "<<compFile->GetName()<<" wont be included in the split sample comparisson plot"<<std::endl;
            continue;
          }

          // if on the same y axis, also check that the contribution of the sample compared to the baseline LLH integral is above the threshold
          // otherwise the labels might get very crowded if the comparisson LLH is much smaller than the baseline one
          if(sameAxis)getSplitSampleStack(compFile, sampleList, paramName, *compLLH_sam, cum_Sums, draw_Label, splitSamplesStack, splitSamplesLegend, LLH_sam->Integral());
          else getSplitSampleStack(compFile, sampleList, paramName, *compLLH_sam, cum_Sums, draw_Label, splitSamplesStack, splitSamplesLegend);

          // go to the pad for the histograms
          if(drawGrid)splitSamplesPad->SetGrid();
          splitSamplesPad->Draw();
          splitSamplesPad->cd();

          splitSamplesStack->Draw(extraDrawOptions.c_str());
          splitSamplesStack->GetXaxis()->SetTitle("Parameter Variation");
          splitSamplesStack->GetYaxis()->SetTitle("-2LLH_{sam}");
          if(sameAxis)splitSamplesStack->SetMaximum(baseSplitSamplesStack->GetMaximum());

          if(totalOnSplitPlots){
            compLLH_sam->Draw(Form("same%s", extraDrawOptions.c_str()));
            compLLH_sam->SetLineColor(kBlack);
            splitSamplesLegend->AddEntry(compLLH_sam, "All Samples", "l");
          }
          splitSamplesLegend->Draw();
          
          // need to draw the labels after other stuff or they dont show up
          for(uint i=0; i < sampleDirList.size(); i++){
            std::string sampName = sampleList[i];
            if(!draw_Label[i])continue;

            label->DrawLatex(compLLH_sam->GetBinLowEdge(compLLH_sam->GetNbinsX()+1), cum_Sums[i], Form("#leftarrow%s", prettify_name(sampName).c_str()) );
          }

          if(plotRatios){
            splitSamplesCanv->cd(2 + extraFileIdx);
            if(drawGrid)splitSamplesRatioPad->SetGrid();
            splitSamplesRatioPad->Draw();
            splitSamplesRatioPad->cd();

            THStack *splitSamplesStackRatios = new THStack(paramName.c_str(), "");

            TList *baselineHistList = baseSplitSamplesStack->GetHists();
            TList *compHistList     = splitSamplesStack->GetHists();

            if(baselineHistList->GetEntries() != compHistList->GetEntries()){
              std::cout<<"WARNING: Number of samples in file "<<compFile->GetName()<<" is not the same as in the baseline file, "<<file1->GetName()<<std::endl;
              std::cout<<"         Cannot do direct comparisson of the contributions to make the ratio plot, will skip."<<std::endl;
              continue;
            }

            for(int sampleIdx = 0; sampleIdx < baselineHistList->GetEntries(); sampleIdx++){
              TH1D *divHist = new TH1D(Form("%s_%s_splitDiv", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()), Form("%s_%s_splitDiv", paramName.c_str(), FileLabels[extraFileIdx+1].c_str()), 
                                      compLLH_sam->GetNbinsX(), compLLH_sam->GetBinLowEdge(1), compLLH_sam->GetBinLowEdge(nBins+1));
              divHist->Divide((TH1D*)compHistList->At(sampleIdx), (TH1D*)baselineHistList->At(sampleIdx));
              splitSamplesStackRatios->Add(divHist);
              divHist->SetLineColor(((TH1D*)compHistList->At(sampleIdx))->GetLineColor());
            }

            splitSamplesStackRatios->Draw(Form("NOSTACK%s", extraDrawOptions.c_str()));

            splitSamplesStackRatios->GetXaxis()->SetLabelSize(ratioLabelScaling * splitSamplesStack->GetXaxis()->GetLabelSize());
            splitSamplesStackRatios->GetXaxis()->SetTitleSize(ratioLabelScaling * splitSamplesStack->GetXaxis()->GetLabelSize());
            splitSamplesStackRatios->GetXaxis()->SetTitle("Parameter Variation");
            
            splitSamplesStackRatios->GetYaxis()->SetLabelSize(ratioLabelScaling * splitSamplesStack->GetYaxis()->GetLabelSize());
            splitSamplesStackRatios->GetYaxis()->SetTitleSize(ratioLabelScaling * splitSamplesStack->GetYaxis()->GetLabelSize());
            splitSamplesStackRatios->GetYaxis()->SetNdivisions(5,2,0);

          }
        }

        splitSamplesCanv->SaveAs(Form("%s_bySample.pdf", PDFName.c_str()));

        delete baseSplitSamplesStack;
        delete baseSplitSamplesLegend;
      }
      
      //delete sampleCanv;   
      delete compStack;    
      delete allSamplesCanv; 
    }
    
    blankCanv->SaveAs(Form("%s.pdf]", PDFName.c_str()));
    if(splitBySample)splitSampBlankCanv->SaveAs(Form("%s_bySample.pdf]", PDFName.c_str()));
                             
    return 0;
}

int main(int argc, char **argv)
{
  SetMaCh3LoggerFormat();

  // parse the inputs
  int c;
  while ((c = getopt(argc, argv, "o:i:l:d:srg"))!= -1) {
    if (c < 0)
      break;
    switch (c) {
      case 'o': {
        PDFName = optarg;
        break;
      }
      case 'i': {
        Mach3LLHFileName = optarg;
        FileLabels_default.push_back(optarg);
        break;
      }
      case 's': {
        splitBySample=true;
        break;
      }
      case 'r': {
        plotRatios = true;
        break;
      }
      case 'g': {
        drawGrid = true;
        break;
      }
      case 'l': {
        parseFileLabels(optarg, FileLabels);
        std::cout<<"INFO: Specified file labels {";
        for(std::string label: FileLabels)std::cout<<label<<", ";
        std::cout<<"}"<<std::endl;
        break;
      }
      case 'd': {
        extraDrawOptions = optarg;
        break;
      }
    }
  }

  std::cout<<std::endl<<"LLH files to compare to mach3 scan: \n{";
  // optind is for the extra arguments that are not parsed by the program
  for(; optind < argc; optind++){
      if(Mach3LLHFileName == ""){ //didnt provide an explicit input file, just take the first on the list
        Mach3LLHFileName = argv[optind];
        FileLabels_default.push_back(argv[optind]);
        continue;
      }

      ExtraFileNames.push_back(argv[optind]);
      std::cout<<argv[optind]<<", ";
      FileLabels_default.push_back(argv[optind]);
  }
  std::cout<<"}"<<std::endl<<std::endl;

  std::cout << "Now plotting LLH scan "<< std::endl;
  if(splitBySample) std::cout<<"Splitting by sample"<<std::endl;

  if(plotRatios && ExtraFileNames.size() == 0){
    std::cerr<<"ERROR: you specified -r <plotRatios> = true but didnt specify any files to compare against, was this a mistake?"<<std::endl;
    throw;
  }

  if(FileLabels.size() == 0){
    std::cout<<"INFO: No file labels specified, will just use the file names"<<std::endl;
    FileLabels = FileLabels_default;
  }

  if((FileLabels.size() != 0) &&(FileLabels.size() != 1+ ExtraFileNames.size())){
    std::cerr<<"ERROR: hmmm, you gave me "<< FileLabels.size() << " labels but "<< ExtraFileNames.size() +1<<" files"<<std::endl;
    std::cerr<<"       that doesn\'t seem right to me, did you forget a file? or a label maybe?"<<std::endl;
    throw;
  }

  return PlotLLH();
}

