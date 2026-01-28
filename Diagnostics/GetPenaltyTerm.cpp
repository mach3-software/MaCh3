// MaCh3 includes
#include "Manager/Manager.h"
#include "Samples/SampleStructs.h"
#include "Samples/HistogramUtils.h"
#include "Parameters/ParameterHandlerUtils.h"

_MaCh3_Safe_Include_Start_ //{
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
#include "TROOT.h"
_MaCh3_Safe_Include_End_ //}

/// @file GetPenaltyTerm.cpp
/// @brief KS: This file contains the implementation of the function to extract specific penalty terms from systematic chains.
///
/// This script is designed to retrieve penalty terms from various sources, such as flux and cross-section systematic chains.
/// Since flux and cross-section uncertainties are handled systematically, the penalty term cannot be taken directly from the chain.
///
/// @todo KS: This should really be moved to MCMC Processor

void ReadCovFile(const std::string& inputFile,
                 std::vector <double>& Prior,
                 std::vector <bool>& isFlat,
                 std::vector<std::string>& ParamNames,
                 std::vector<std::vector<double>>& invCovMatrix,
                 int& nParams)
{
  // Now read the MCMC file
  TFile *TempFile = M3::Open(inputFile, "open", __FILE__, __LINE__);

  // Get the matrix
  TDirectory* CovarianceFolder = TempFile->Get<TDirectory>("CovarianceFolder");
  TMatrixDSym *CovMatrix = M3::GetCovMatrixFromChain(CovarianceFolder);

  // Get the settings for the MCMC
  TMacro *Config = TempFile->Get<TMacro>("MaCh3_Config");
  if (Config == nullptr) {
    MACH3LOG_ERROR("Didn't find MaCh3_Config tree in MCMC file! {}", inputFile);
    TempFile->ls();
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  YAML::Node Settings = TMacroToYAML(*Config);

  //CW: Get the Covariance matrix
  std::vector<std::string> CovPos = GetFromManager<std::vector<std::string>>(Settings["General"]["Systematics"]["XsecCovFile"], {"none"});
  if(CovPos.back() == "none")
  {
    MACH3LOG_WARN("Couldn't find Cov branch in output");
    MaCh3Utils::PrintConfig(Settings);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  //KS:Most inputs are in ${MACH3}/inputs/blarb.root
  if (std::getenv("MACH3") != nullptr) {
    MACH3LOG_INFO("Found MACH3 environment variable: {}", std::getenv("MACH3"));
    for(unsigned int i = 0; i < CovPos.size(); i++)
      CovPos[i].insert(0, std::string(std::getenv("MACH3"))+"/");
  }

  YAML::Node CovFile;
  CovFile["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for(unsigned int i = 0; i < CovPos.size(); i++)
  {
    YAML::Node YAMLDocTemp = M3OpenConfig(CovPos[i]);
    for (const auto& item : YAMLDocTemp["Systematics"]) {
      CovFile["Systematics"].push_back(item);
    }
  }

  nParams = CovMatrix->GetNrows();

  auto systematics = CovFile["Systematics"];
  for (auto it = systematics.begin(); it != systematics.end(); ++it)
  {
    auto const &param = *it;

    ParamNames.push_back(param["Systematic"]["Names"]["FancyName"].as<std::string>());
    Prior.push_back( param["Systematic"]["ParameterValues"]["PreFitValue"].as<double>() );

    bool flat = false;
    if (param["Systematic"]["FlatPrior"]) { flat = param["Systematic"]["FlatPrior"].as<bool>(); }
    isFlat.push_back( flat );
  }

  CovMatrix->Invert();
  //KS: Let's use double as it is faster than TMatrix
  invCovMatrix.resize(nParams, std::vector<double>(nParams, -999));

  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2)
  #endif
  for (int i = 0; i < nParams; i++)
  {
    for (int j = 0; j < nParams; ++j)
    {
      invCovMatrix[i][j] = (*CovMatrix)(i,j);
    }
  }

  TempFile->Close();
  delete TempFile;
}

void LoadSettings(YAML::Node& Settings,
                  std::vector<std::string>& SetsNames,
                  std::vector<std::string>& FancyTitle,
                  std::vector<std::vector<bool>>& isRelevantParam,
                  const std::vector<std::string>& ParamNames,
                  const int nParams)
{
  std::vector<std::string> node = Settings["GetPenaltyTerm"]["PenaltySets"].as<std::vector<std::string>>();
  std::vector<std::vector<std::string>> RemoveNames;
  std::vector<bool> Exclude;

  for (unsigned int i = 0; i < node.size(); i++)
  {
    std::string ParName = node[i];
    SetsNames.push_back(ParName);

    const auto& Set = Settings["GetPenaltyTerm"][ParName];

    RemoveNames.push_back(Set[0].as<std::vector<std::string>>());
    Exclude.push_back(Set[1].as<bool>());
    FancyTitle.push_back(Set[2].as<std::string>());
  }

  const int NSets = int(SetsNames.size());

  isRelevantParam.resize(NSets);
  //Loop over sets in the config
  for(int i = 0; i < NSets; i++)
  {
    isRelevantParam[i].resize(nParams);
    int counter = 0;
    //Loop over parameters in the Covariance object
    for (int j = 0; j < nParams; j++)
    {
      isRelevantParam[i][j] = false;

      //KS: Here we loop over all names and if parameters wasn't matched then we set it is relevant.
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
    MACH3LOG_INFO(" Found {} params for set {}", counter, SetsNames[i]);
  }
}

void GetPenaltyTerm(const std::string& inputFile, const std::string& configFile)
{
  auto canvas = std::make_unique<TCanvas>("canvas", "canvas", 0, 0, 1024, 1024);
  canvas->SetGrid();
  canvas->SetTickx();
  canvas->SetTicky();

  canvas->SetBottomMargin(0.1f);
  canvas->SetTopMargin(0.02f);
  canvas->SetRightMargin(0.08f);
  canvas->SetLeftMargin(0.15f);

  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  gStyle->SetPalette(51);

  std::vector <double> Prior;
  std::vector <bool> isFlat;
  std::vector<std::string> ParamNames;
  std::vector<std::vector<double>> invCovMatrix;
  int nParams;
  ReadCovFile(inputFile, Prior, isFlat, ParamNames, invCovMatrix, nParams);

  std::vector<TString> BranchNames;

  // Open the Chain
  TChain* Chain = new TChain("posteriors","");
  Chain->Add(inputFile.c_str()); 
    
  // Get the list of branches
  TObjArray* brlis = Chain->GetListOfBranches();

  // Get the number of branches
  int nBranches = brlis->GetEntries();
  int RelevantBranches = 0;
  for (int i = 0; i < nBranches; i++) 
  {
    // Get the TBranch and its name
    TBranch* br = static_cast<TBranch*>(brlis->At(i));
    if(!br){
      MACH3LOG_ERROR("Invalid branch at position {}", i);
      throw MaCh3Exception(__FILE__,__LINE__);
    }
    TString bname = br->GetName();

    // If we're on beam systematics
    if(bname.BeginsWith("param_"))
    {
      BranchNames.push_back(bname);
      RelevantBranches++;
    }
  }
  
  // Set all the branches to off
  Chain->SetBranchStatus("*", false);
  
  std::vector<double> fParProp(RelevantBranches);
  // Turn on the branches which we want for parameters
  for (int i = 0; i < RelevantBranches; ++i) 
  {
    Chain->SetBranchStatus(BranchNames[i].Data(), true);
    Chain->SetBranchAddress(BranchNames[i].Data(), &fParProp[i]);
  }

  YAML::Node Settings = M3OpenConfig(configFile);
  std::vector<std::string> SetsNames;
  std::vector<std::string> FancyTitle;
  std::vector<std::vector<bool>> isRelevantParam;

  LoadSettings(Settings, SetsNames, FancyTitle, isRelevantParam, ParamNames, nParams);

  const int NSets = int(SetsNames.size());
  int AllEvents = int(Chain->GetEntries());
  std::vector<std::unique_ptr<TH1D>> hLogL(NSets);
  for (int i = 0; i < NSets; i++) {
    std::string NameTemp = "LogL_" + SetsNames[i];
    hLogL[i] = std::make_unique<TH1D>(NameTemp.c_str(), NameTemp.c_str(), AllEvents, 0, AllEvents);
    hLogL[i]->SetLineColor(kBlue);
  }
  std::vector<double> logL(NSets, 0.0);
  for(int n = 0; n < AllEvents; ++n)
  {
    if(n%10000 == 0) MaCh3Utils::PrintProgressBar(n, AllEvents);
      
    Chain->GetEntry(n);

    for(int k = 0; k < NSets; ++k) logL[k] = 0.;
#ifdef MULTITHREAD
    // The per-thread array
    double *logL_private = nullptr;

    // Declare the omp parallel region
    // The parallel region needs to stretch beyond the for loop!
    #pragma omp parallel private(logL_private)
    {
      logL_private = new double[NSets];
      for(int k = 0; k < NSets; ++k) logL_private[k] = 0.;

      #pragma omp for
      for (int i = 0; i < nParams; i++)
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
                //KS: Since matrix is symmetric we can calculate non diagonal elements only once and multiply by 2, can bring up to factor speed decrease.
                int scale = 1;
                if(i != j) scale = 2;
                logL_private[k] += scale * 0.5*(fParProp[i] - Prior[i])*(fParProp[j] - Prior[j])*invCovMatrix[i][j];
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
      //Delete private arrays
      delete[] logL_private;
    }//End omp range

#else
    for (int i = 0; i < nParams; i++)
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
              //KS: Since matrix is symmetric we can calculate non diagonal elements only once and multiply by 2, can bring up to factor speed decrease.
              int scale = 1;
              if(i != j) scale = 2;
              logL[k] += scale * 0.5*(fParProp[i] - Prior[i])*(fParProp[j] - Prior[j])*invCovMatrix[i][j];
            }
          }
        }
      }
    }
#endif // end MULTITHREAD
    for(int k = 0; k < NSets; ++k)
    {
      hLogL[k]->SetBinContent(n, logL[k]);
    }
  }//End loop over steps

  // Directory for posteriors
  std::string OutputName = inputFile + "_PenaltyTerm" +".root";
  TFile *OutputFile = M3::Open(OutputName, "recreate", __FILE__, __LINE__);
  TDirectory *PenaltyTermDir = OutputFile->mkdir("PenaltyTerm");

  canvas->Print(Form("%s_PenaltyTerm.pdf[",inputFile.c_str()), "pdf");
  for(int i = 0; i < NSets; i++)
  {
    const double Maximum = hLogL[i]->GetMaximum();
    hLogL[i]->GetYaxis()->SetRangeUser(0., Maximum*1.2);
    hLogL[i]->SetTitle(FancyTitle[i].c_str());
    hLogL[i]->GetXaxis()->SetTitle("Step");
    hLogL[i]->GetYaxis()->SetTitle(FancyTitle[i].c_str());
    hLogL[i]->GetYaxis()->SetTitleOffset(1.4f);

    hLogL[i]->Draw("");

    PenaltyTermDir->cd();
    hLogL[i]->Write();

    canvas->Print(Form("%s_PenaltyTerm.pdf",inputFile.c_str()), "pdf");
  }
  canvas->Print(Form("%s_PenaltyTerm.pdf]",inputFile.c_str()), "pdf");
  delete Chain;

  OutputFile->Close();
  delete OutputFile;
}

int main(int argc, char *argv[])
{
  SetMaCh3LoggerFormat();
  MaCh3Utils::MaCh3Welcome();
  if (argc != 3 )
  {
    MACH3LOG_WARN("Something went wrong ");
    MACH3LOG_WARN("{} root_file_to_analyse.root", argv[1]);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  std::string filename = argv[1];
  std::string config = argv[2];
  GetPenaltyTerm(filename, config);

  return 0;
}
