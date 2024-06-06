// C++ includes
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

// ROOT includes
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TMatrixT.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TObject.h"
#include "TChain.h"
#include "TH1.h"
#include "TColor.h"
#include "TObjString.h"

#ifdef MULTITHREAD
#include "omp.h"
#endif

#include "manager/manager.h"

//KS: Simple script to get whatever penalty term you like, since flux and xsec are on systematic we cannot just directly take it from the chain
//g++ `root-config --cflags` -std=c++11 -g -o GetPenaltyTerm GetPenaltyTerm.cpp -I`root-config --incdir` `root-config --glibs --libs`

void GetPenaltyTerm(std::string inputFile, std::string configFile);
void ReadXSecFile(std::string inputFile);

std::vector <int> nominal;
std::vector <bool> isFlat;
std::vector<TString> BranchNames;
std::vector<std::string> ParamNames;
int size;

double** invCovMatrix;

int main(int argc, char *argv[]) 
{
  SetMaCh3LoggerFormat();
  if (argc != 3 )
  {
    std::cerr<< " Something went wrong " << std::endl;
    std::cerr << "./GetPenaltyTerm root_file_to_analyse.root " << std::endl;
    exit(-1);
  }
  std::string filename = argv[1];
  std::string config = argv[2];
  GetPenaltyTerm(filename, config);

  return 0;
}

//KS: This should really be moved to MCMC Processor
void GetPenaltyTerm(std::string inputFile, std::string configFile)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 0, 0, 1024, 1024);
  canvas->SetGrid();
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetBottomMargin(0.1);
  canvas->SetTopMargin(0.02);
  canvas->SetRightMargin(0.08);
  canvas->SetLeftMargin(0.15);

  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  gStyle->SetPalette(51);

  ReadXSecFile(inputFile);
    
  // Open the Chain
  TChain* Chain = new TChain("posteriors","");
  Chain->Add(inputFile.c_str()); 
    
  // Get the list of branches
  TObjArray* brlis = (TObjArray*)(Chain->GetListOfBranches());

  // Get the number of branches
  int nBranches = brlis->GetEntries();
  int RelevantBranches = 0;
  for (int i = 0; i < nBranches; i++) 
  {
    // Get the TBranch and its name
    TBranch* br = (TBranch*)brlis->At(i);
    TString bname = br->GetName();

    // If we're on beam systematics
    if(bname.BeginsWith("xsec_"))
    {
      BranchNames.push_back(bname);
      RelevantBranches++;
    }
  }
  
  // Set all the branches to off
  Chain->SetBranchStatus("*", false);
  
  double* fParProp = new double[RelevantBranches];
  // Turn on the branches which we want for parameters
  for (int i = 0; i < RelevantBranches; ++i) 
  {
    Chain->SetBranchStatus(BranchNames[i].Data(), true);
    Chain->SetBranchAddress(BranchNames[i].Data(), &fParProp[i]);
  }

  YAML::Node Settings = YAML::LoadFile(configFile.c_str());

  std::vector<std::string> SetsNames;
  std::vector<std::vector<std::string>> RemoveNames;
  std::vector<bool> Exclude;
  std::vector<std::string> FancyTittle;

  std::vector<std::vector<bool>> isRelevantParam;
  std::vector<std::string> node = Settings["GetPenaltyTerm"]["PenaltySets"].as<std::vector<std::string>>();
  for (unsigned int i = 0; i < node.size(); i++)
  {
    std::string ParName = node[i];
    SetsNames.push_back(ParName);

    const auto& Set = Settings["GetPenaltyTerm"][ParName];

    RemoveNames.push_back(Set[0].as<std::vector<std::string>>());
    Exclude.push_back(Set[1].as<bool>());
    FancyTittle.push_back(Set[2].as<std::string>());
  }

  const int NSets = SetsNames.size();

  isRelevantParam.resize(NSets);
  //Loop over sets in the config
  for(int i = 0; i < NSets; i++)
  {
    isRelevantParam[i].resize(size);
    int counter = 0;
    //Loop over paramters in the Covariance object
    for (int j = 0; j < size; j++)
    {
      isRelevantParam[i][j] = false;

      //KS: Here we loop over all names and if parameters wasn't matched then we set it is relevant. For xsec it is easier to do it liek this
      if(Exclude[i])
      {
        bool found = false;
        for (unsigned int k = 0; k < RemoveNames[i].size(); k++)
        {
          if (ParamNames[j].rfind(RemoveNames[i][k], 0) == 0)
          {
            found = true;
          }
        }
        if(!found)
        {
          isRelevantParam[i][j] = true;
          counter++;
        }
      }
      //KS: Here is much simpler, if parameter matched then it is relevant
      else
      {
        for (unsigned int k = 0; k < RemoveNames[i].size(); k++)
        {
          if (ParamNames[j].rfind(RemoveNames[i][k], 0) == 0)
          {
            isRelevantParam[i][j] = true;
            counter++;
            break;
          }
        }
      }
    }
    std::cout<<" Found "<<counter<<" params for set "<<SetsNames[i]<<std::endl;
  }

  int AllEvents = Chain->GetEntries();
  TH1D **hLogL = new TH1D *[NSets];
  for(int i = 0; i < NSets; i++)
  {
    std::string NameTemp = "LogL_" + SetsNames[i];
    hLogL[i] = new TH1D(NameTemp.c_str(), NameTemp.c_str(), AllEvents, 0 , AllEvents);
    hLogL[i]->SetLineColor(kBlue);
  }
  double* logL = new double[NSets]();
  for(int n = 0;  n < AllEvents; n++)
  {
    if(n%10000 == 0) std::cout<<n<<std::endl;
      
    Chain->GetEntry(n);

    for(int k = 0; k < NSets; ++k) logL[k] = 0.;
#ifdef MULTITHREAD
    // The per-thread array
    double *logL_private = nullptr;

    // Declare the omp parallel region
    // The parallel region needs to stretch beyond the for loop!
    #pragma omp parallel private(logL_private)
    {
      logL_private = new double[NSets]();
      for(int k = 0; k < NSets; ++k) logL_private[k] = 0.;

      #pragma omp for
      for (int i = 0; i < size; i++)
      {
        for (int j = 0; j <= i; ++j)
        {
          //check if flat prior
          if (!isFlat[i] && !isFlat[j])
          {
            for(int k = 0; k < NSets; ++k)
            {
              //Check if parameter is relevant for this set
              if (isRelevantParam[k][i] && isRelevantParam[k][j])
              {
                //KS: Since matrix is symmetric we can calcaute non diagonal elements only once and multiply by 2, can bring up to factor speed decrease.
                int scale = 1;
                if(i != j) scale = 2;
                logL_private[k] += scale * 0.5*(fParProp[i] - nominal[i])*(fParProp[j] - nominal[j])*invCovMatrix[i][j];
              }
            }
          }
        }
      }
      // Now we can write the individual arrays from each thread to the main array
      for(int k = 0; k < NSets; ++k)
      {
        #pragma omp atomic
        logL[k] += logL_private[k];
      }
      //Delete private arrayss
      delete[] logL_private;
    }//End omp range

#else
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j <= i; ++j)
      {
        //check if flat prior
        if (!isFlat[i] && !isFlat[j])
        {
          for(int k = 0; k < NSets; ++k)
          {
            //Check if parameter is relevant for this set
            if (isRelevantParam[k][i] && isRelevantParam[k][j])
            {
              //KS: Since matrix is symetric we can calcaute non daigonal elements only once and multiply by 2, can bring up to factor speed decrease.
              int scale = 1;
              if(i != j) scale = 2;
              logL[k] += scale * 0.5*(fParProp[i] - nominal[i])*(fParProp[j] - nominal[j])*invCovMatrix[i][j];
            }
          }
        }
      }
    }
#endif
    for(int k = 0; k < NSets; ++k)
    {
      hLogL[k]->SetBinContent(n, logL[k]);
    }
  }//End loop over steps

  delete[] logL;
  // Directory for posteriors
  std::string OutputName = inputFile + "_PenaltyTerm" +".root";
  TFile* OutputFile = new TFile(OutputName.c_str(), "recreate");
  TDirectory *PenaltyTermDir = OutputFile->mkdir("PenaltyTerm");

  canvas->Print(Form("%s_PenaltyTerm.pdf[",inputFile.c_str()), "pdf");
  for(int i = 0; i < NSets; i++)
  {
    double Maximum = hLogL[i]->GetMaximum();
    hLogL[i]->GetYaxis()->SetRangeUser(0., Maximum*1.2);
    hLogL[i]->SetTitle(FancyTittle[i].c_str());
    hLogL[i]->GetXaxis()->SetTitle("Step");
    hLogL[i]->GetYaxis()->SetTitle(FancyTittle[i].c_str());
    hLogL[i]->GetYaxis()->SetTitleOffset(1.4);

    hLogL[i]->Draw("");

    PenaltyTermDir->cd();
    hLogL[i]->Write();

    canvas->Print(Form("%s_PenaltyTerm.pdf",inputFile.c_str()), "pdf");
    delete hLogL[i];
  }
  canvas->Print(Form("%s_PenaltyTerm.pdf]",inputFile.c_str()), "pdf");
  delete[] hLogL;
  delete[] fParProp;
  delete Chain;
  delete canvas;

  for (int i = 0; i < size; i++)
  {
    delete[] invCovMatrix[i];
  }
  delete[] invCovMatrix;
}

void ReadXSecFile(std::string inputFile)
{
  // Now read the MCMC file
  TFile *TempFile = new TFile(inputFile.c_str(), "open");

  // Get the matrix
  TMatrixDSym *XSecMatrix = (TMatrixDSym*)(TempFile->Get("CovarianceFolder/xsec_cov"));

  // Get the settings for the MCMC
  TMacro *Config = (TMacro*)(TempFile->Get("MaCh3_Config"));
  if (Config == nullptr) {
    MACH3LOG_ERROR("Didn't find MaCh3_Config tree in MCMC file! {}", inputFile);
    TempFile->ls();
    throw;
  }

  YAML::Node Settings = TMacroToYAML(*Config);

  //CW: Get the xsec Covariance matrix
  std::vector<std::string> XsecCovPos = GetFromManager<std::vector<std::string>>(Settings["General"]["Systematics"]["XsecCovFile"], {"none"});
  if(XsecCovPos.back() == "none")
  {
    MACH3LOG_WARN("Couldn't find XsecCov branch in output");
    std::cout<<Settings<<std::endl;
    throw;
  }

  //KS:Most inputs are in ${MACH3}/inputs/blarb.root
  if (std::getenv("MACH3") != nullptr) {
    MACH3LOG_INFO("Found MACH3 environment variable: {}", std::getenv("MACH3"));
    for(unsigned int i = 0; i < XsecCovPos.size(); i++)
      XsecCovPos[i].insert(0, std::string(std::getenv("MACH3"))+"/");
  }

  YAML::Node XSecFile;
  XSecFile["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for(unsigned int i = 0; i < XsecCovPos.size(); i++)
  {
    YAML::Node YAMLDocTemp = YAML::LoadFile(XsecCovPos[i]);
    for (const auto& item : YAMLDocTemp["Systematics"]) {
      XSecFile["Systematics"].push_back(item);
    }
  }


  size = XSecMatrix->GetNrows();

  auto systematics = XSecFile["Systematics"];
  for (auto it = systematics.begin(); it != systematics.end(); ++it)
  {
    auto const &param = *it;

    ParamNames.push_back(param["Systematic"]["Names"]["FancyName"].as<std::string>());
    nominal.push_back( param["Systematic"]["ParameterValues"]["PreFitValue"].as<double>() );

    bool flat = false;
    if (param["Systematic"]["FlatPrior"]) { flat = param["Systematic"]["FlatPrior"].as<bool>(); }
    isFlat.push_back( flat );
  }
    
  XSecMatrix->Invert();
  //KS: Let's use double as it is faster than TMatrix
  invCovMatrix = new double*[size];
  for (int i = 0; i < size; i++) 
  {
    invCovMatrix[i] = new double[size];
    for (int j = 0; j < size; ++j) 
    {
      invCovMatrix[i][j] = -999;
    }
  }
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; ++j) 
    {
      invCovMatrix[i][j] = (*XSecMatrix)(i,j);
    }
  }  
  
  TempFile->Close();
  delete TempFile;
}
