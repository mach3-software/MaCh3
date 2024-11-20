// ROOT includes
#include "TObjArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TColor.h"
#include "TStyle.h"
#include "TROOT.h"

#ifdef MULTITHREAD
#include "omp.h"
#endif

// MaCh3 includes
#include "manager/manager.h"

/// @file RHat.cpp
/// @brief This executable calculates the \f$ \hat{R} \f$ estimator for Markov Chain Monte Carlo (MCMC) convergence.
///
/// KS: This exe is meant to calculate the \f$ \hat{R} \f$ estimator. For a well-converged chain, this distribution
/// should be centered at one. The \f$ \hat{R} \f$ statistic is used to assess the convergence of MCMC simulations
/// and helps determine whether the chains have reached a stable distribution.
///
/// @cite gelman2019.

// *******************
int Ntoys;
int Nchains;

int nDraw;

std::vector<TString> BranchNames;
std::vector<std::string> MCMCFile;
std::vector<bool> ValidPar;

double ***Draws;

double** Mean;
double** StandardDeviation;

double* MeanGlobal;
double* StandardDeviationGlobal;

double* BetweenChainVariance;
double* MarginalPosteriorVariance;
double* RHat;
double* EffectiveSampleSize;

double ***DrawsFolded;
double* MedianArr;

double** MeanFolded;
double** StandardDeviationFolded;

double* MeanGlobalFolded;
double* StandardDeviationGlobalFolded;

double* BetweenChainVarianceFolded;
double* MarginalPosteriorVarianceFolded;
double* RHatFolded;
double* EffectiveSampleSizeFolded;
// *******************
void PrepareChains();
void InitialiseArrays();

void RunDiagnostic();
void CalcRhat();

void SaveResults();
void DestroyArrays();
double CalcMedian(double arr[], int size);

void CapVariable(double var, double cap);

// *******************
int main(int argc, char *argv[]) {
// *******************

  SetMaCh3LoggerFormat();
  MaCh3Utils::MaCh3Welcome();

  Draws = nullptr;
  Mean = nullptr;
  StandardDeviation = nullptr;

  MeanGlobal = nullptr;
  StandardDeviationGlobal = nullptr;

  BetweenChainVariance = nullptr;
  MarginalPosteriorVariance = nullptr;
  RHat = nullptr;
  EffectiveSampleSize = nullptr;

  DrawsFolded = nullptr;
  MedianArr = nullptr;
  MeanFolded = nullptr;
  StandardDeviationFolded = nullptr;

  MeanGlobalFolded = nullptr;
  StandardDeviationGlobalFolded = nullptr;

  BetweenChainVarianceFolded = nullptr;
  MarginalPosteriorVarianceFolded = nullptr;
  RHatFolded = nullptr;
  EffectiveSampleSizeFolded = nullptr;

  Nchains = 0;

  if (argc == 1 || argc == 2)
  {
    MACH3LOG_ERROR("Wrong arguments");
    MACH3LOG_ERROR("./RHat Ntoys MCMCchain_1.root MCMCchain_2.root MCMCchain_3.root ... [how many you like]");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  Ntoys = atoi(argv[1]);

  //KS Gelman suggests to diagnose on more than one chain
  for (int i = 2; i < argc; i++)
  {
    MCMCFile.push_back(std::string(argv[i]));
    MACH3LOG_INFO("Adding file: {}", MCMCFile.back());
    Nchains++;
  }

  if(Ntoys < 1)
  {
    MACH3LOG_ERROR("You specified {} specify larger greater than 0", Ntoys);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if(Nchains == 1)
  {
    MACH3LOG_WARN("Gelman is going to be sad :(. He suggested you should use more than one chain (at least 4). Code works fine for one chain, however, estimator might be biased.");
    MACH3LOG_WARN("Multiple chains are more likely to reveal multimodality and poor adaptation or mixing:");
  }
  MACH3LOG_INFO("Diagnosing {} chains, with {} toys", Nchains, Ntoys);

  PrepareChains();

  InitialiseArrays();

  //KS: Main function
  RunDiagnostic();

  SaveResults();

  DestroyArrays();

  return 0;
}

// *******************
// Load chain and prepare toys
void PrepareChains() {
// *******************
  auto rnd = std::make_unique<TRandom3>(0);

  MACH3LOG_INFO("Generating {}", Ntoys);

  TStopwatch clock;
  clock.Start();

  std::vector<int> BurnIn(Nchains);
  std::vector<int> nEntries(Nchains);
  std::vector<int> nBranches(Nchains);
  std::vector<int> step(Nchains);

  Draws = new double**[Nchains]();
  DrawsFolded = new double**[Nchains]();

  // KS: This can reduce time necessary for caching even by half
  #ifdef MULTITHREAD
  //ROOT::EnableImplicitMT();
  #endif

  // Open the Chain
  //It is tempting to multithread here but unfortunately, ROOT files are not thread safe :(
  for (int m = 0; m < Nchains; m++)
  {
    TChain* Chain = new TChain("posteriors");
    Chain->Add(MCMCFile[m].c_str());
    MACH3LOG_INFO("On file: {}", MCMCFile[m].c_str());
    nEntries[m] = int(Chain->GetEntries());

    // Set the step cut to be 20%
    BurnIn[m] = nEntries[m]/5;

    // Get the list of branches
    TObjArray* brlis = Chain->GetListOfBranches();

    // Get the number of branches
    nBranches[m] = brlis->GetEntries();

    if(m == 0) BranchNames.reserve(nBranches[m]);

    // Set all the branches to off
    Chain->SetBranchStatus("*", false);

    // Loop over the number of branches
    // Find the name and how many of each systematic we have
    for (int i = 0; i < nBranches[m]; i++)
    {
      // Get the TBranch and its name
      TBranch* br = static_cast<TBranch *>(brlis->At(i));
      if(!br){
        MACH3LOG_ERROR("Invalid branch at position {}", i);
        throw MaCh3Exception(__FILE__,__LINE__);
      }
      TString bname = br->GetName();

      // Read in the step
      if (bname == "step") {
        Chain->SetBranchStatus(bname, true);
        Chain->SetBranchAddress(bname, &step[m]);
      }
      //Count all branches
      else if (bname.BeginsWith("PCA_") || bname.BeginsWith("accProb") || bname.BeginsWith("stepTime") )
      {
        continue;
      }
      else
      {
        //KS: Save branch name only for one chain, we assume all chains have the same branches, otherwise this doesn't make sense either way
        if(m == 0)
        {
          BranchNames.push_back(bname);
          //KS: We calculate R Hat also for LogL, just in case, however we plot them separately
          if(bname.BeginsWith("LogL"))
          {
            ValidPar.push_back(false);
          }
          else
          {
            ValidPar.push_back(true);
          }
        }
        Chain->SetBranchStatus(bname, true);
        MACH3LOG_DEBUG("{}", bname);
      }
    }

    if(m == 0) nDraw = int(BranchNames.size());

    //TN: Qualitatively faster sanity check, with the very same outcome (all chains have the same #branches)
    if(m > 0)
    {
      if(nBranches[m] != nBranches[0])
      {
        MACH3LOG_ERROR("Ups, something went wrong, chain {} called {} has {} branches, while 0 called {} has {} branches", m, MCMCFile[m], nBranches[m], MCMCFile[0], nBranches[0]);
        MACH3LOG_ERROR("All chains should have the same number of branches");
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
    }

    //TN: move the Draws here, so we need to iterate over every chain only once
    Draws[m] = new double*[Ntoys]();
    DrawsFolded[m] = new double*[Ntoys]();
    for(int i = 0; i < Ntoys; i++)
    {
      Draws[m][i] = new double[nDraw]();
      DrawsFolded[m][i] = new double[nDraw]();
      for(int j = 0; j < nDraw; j++)
      {
        Draws[m][i][j] = 0.;
        DrawsFolded[m][i][j] = 0.;
      }
    }

    //TN: move looping over toys here, so we don't need to loop over chains more than once
    if(BurnIn[m] >= nEntries[m])
    {
      MACH3LOG_ERROR("You are running on a chain shorter than BurnIn cut");
      MACH3LOG_ERROR("Number of entries {} BurnIn cut {}", nEntries[m], BurnIn[m]);
      MACH3LOG_ERROR("You will run into the infinite loop");
      MACH3LOG_ERROR("You can make a new chain or modify BurnIn cut");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    for (int i = 0; i < Ntoys; i++)
    {
      // Get a random entry after burn in
      int entry = int(nEntries[m]*rnd->Rndm());

      Chain->GetEntry(entry);

      // If we have combined chains by hadd need to check the step in the chain
      // Note, entry is not necessarily the same as the step due to merged ROOT files, so can't choose an entry in the range BurnIn - nEntries :(
      if (step[m] < BurnIn[m])
      {
        i--;
        continue;
      }

      // Output some info for the user
      if (Ntoys > 10 && i % (Ntoys/10) == 0) {
        MaCh3Utils::PrintProgressBar(i+m*Ntoys, static_cast<Long64_t>(Ntoys)*Nchains);
        MACH3LOG_DEBUG("Getting random entry {}", entry);
      }

      // Set the branch addresses for params
      for (int j = 0; j < nDraw; ++j)
      {
        Chain->SetBranchAddress(BranchNames[j].Data(), &Draws[m][i][j]);
      }
      Chain->GetEntry(entry);

    }//end loop over toys

    //TN: There, we now don't need to keep the chain in memory anymore
    delete Chain;
  }

  //KS: Now prepare folded draws, quoting Gelman
  //"We propose to report the maximum of rank normalized split-Rb and rank normalized folded-split-Rb for each parameter"
  MedianArr = new double[nDraw]();
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int j = 0; j < nDraw; j++)
  {
    MedianArr[j] = 0.;
    std::vector<double> TempDraws(static_cast<size_t>(Ntoys) * Nchains);
    for(int m = 0; m < Nchains; m++)
    {
      for(int i = 0; i < Ntoys; i++)
      {
        const int im = i+m;
        TempDraws[im] = Draws[m][i][j];
      }
    }
    MedianArr[j] = CalcMedian(TempDraws.data(), Ntoys*Nchains);
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(3)
  #endif
  for(int m = 0; m < Nchains; m++)
  {
    for(int i = 0; i < Ntoys; i++)
    {
      for(int j = 0; j < nDraw; j++)
      {
        DrawsFolded[m][i][j] = std::fabs(Draws[m][i][j] - MedianArr[j]);
      }
    }
  }
  clock.Stop();
  MACH3LOG_INFO("Finished calculating Toys, it took {:.2f}s to finish", clock.RealTime());
}

// *******************
// Create all arrays we are going to use later
void InitialiseArrays() {
// *******************

  MACH3LOG_INFO("Initialising arrays");
  Mean = new double*[Nchains]();
  StandardDeviation = new double*[Nchains]();

  MeanGlobal = new double[nDraw]();
  StandardDeviationGlobal = new double[nDraw]();
  BetweenChainVariance = new double[nDraw]();

  MarginalPosteriorVariance = new double[nDraw]();
  RHat = new double[nDraw]();
  EffectiveSampleSize = new double[nDraw]();

  MeanFolded = new double*[Nchains]();
  StandardDeviationFolded = new double*[Nchains]();

  MeanGlobalFolded = new double[nDraw]();
  StandardDeviationGlobalFolded = new double[nDraw]();
  BetweenChainVarianceFolded = new double[nDraw]();

  MarginalPosteriorVarianceFolded = new double[nDraw]();
  RHatFolded = new double[nDraw]();
  EffectiveSampleSizeFolded = new double[nDraw]();

  for (int m = 0; m < Nchains; ++m)
  {
    Mean[m] = new double[nDraw]();
    StandardDeviation[m] = new double[nDraw]();

    MeanFolded[m] = new double[nDraw]();
    StandardDeviationFolded[m] = new double[nDraw]();
    for (int j = 0; j < nDraw; ++j)
    {
      Mean[m][j] = 0.;
      StandardDeviation[m][j] = 0.;

      MeanFolded[m][j] = 0.;
      StandardDeviationFolded[m][j] = 0.;
      if(m == 0)
      {
        MeanGlobal[j] = 0.;
        StandardDeviationGlobal[j] = 0.;
        BetweenChainVariance[j] = 0.;
        MarginalPosteriorVariance[j] = 0.;
        RHat[j] = 0.;
        EffectiveSampleSize[j] = 0.;

        MeanGlobalFolded[j] = 0.;
        StandardDeviationGlobalFolded[j] = 0.;
        BetweenChainVarianceFolded[j] = 0.;
        MarginalPosteriorVarianceFolded[j] = 0.;
        RHatFolded[j] = 0.;
        EffectiveSampleSizeFolded[j] = 0.;
      }
    }
  }
}

// *******************
void RunDiagnostic() {
// *******************
  CalcRhat();
  //In case in future we expand this
}

// *******************
//KS: Based on Gelman et. al. arXiv:1903.08008v5
// Probably most of it could be moved cleverly to MCMC Processor, keep it separate for now
void CalcRhat() {
// *******************

  TStopwatch clock;
  clock.Start();

  //KS: Start parallel region
  // If we would like to do this for thousands of chains we might consider using GPU for this
  #ifdef MULTITHREAD
  #pragma omp parallel
  {
  #endif

    #ifdef MULTITHREAD
    #pragma omp for collapse(2)
    #endif
    //KS: loop over chains and draws are independent so might as well collapse for sweet cache hits
    //Calculate the mean for each parameter within each considered chain
    for (int m = 0; m < Nchains; ++m)
    {
      for (int j = 0; j < nDraw; ++j)
      {
        for(int i = 0; i < Ntoys; i++)
        {
          Mean[m][j] += Draws[m][i][j];
          MeanFolded[m][j] += DrawsFolded[m][i][j];
        }
        Mean[m][j] = Mean[m][j]/Ntoys;
        MeanFolded[m][j] = MeanFolded[m][j]/Ntoys;
      }
    }

    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    //Calculate the mean for each parameter global means we include information from several chains
    for (int j = 0; j < nDraw; ++j)
    {
      for (int m = 0; m < Nchains; ++m)
      {
        MeanGlobal[j] += Mean[m][j];
        MeanGlobalFolded[j] += MeanFolded[m][j];
      }
      MeanGlobal[j] = MeanGlobal[j]/Nchains;
      MeanGlobalFolded[j] = MeanGlobalFolded[j]/Nchains;
    }


    #ifdef MULTITHREAD
    #pragma omp for collapse(2)
    #endif
    //Calculate the standard deviation for each parameter within each considered chain
    for (int m = 0; m < Nchains; ++m)
    {
      for (int j = 0; j < nDraw; ++j)
      {
        for(int i = 0; i < Ntoys; i++)
        {
          StandardDeviation[m][j] += (Draws[m][i][j] - Mean[m][j])*(Draws[m][i][j] - Mean[m][j]);
          StandardDeviationFolded[m][j] += (DrawsFolded[m][i][j] - MeanFolded[m][j])*(DrawsFolded[m][i][j] - MeanFolded[m][j]);
        }
        StandardDeviation[m][j] = StandardDeviation[m][j]/(Ntoys-1);
        StandardDeviationFolded[m][j] = StandardDeviationFolded[m][j]/(Ntoys-1);
      }
    }

    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    //Calculate the standard deviation for each parameter combining information from all chains
    for (int j = 0; j < nDraw; ++j)
    {
      for (int m = 0; m < Nchains; ++m)
      {
        StandardDeviationGlobal[j] += StandardDeviation[m][j];
        StandardDeviationGlobalFolded[j] += StandardDeviationFolded[m][j];
      }
      StandardDeviationGlobal[j] = StandardDeviationGlobal[j]/Nchains;
      StandardDeviationGlobalFolded[j] = StandardDeviationGlobalFolded[j]/Nchains;
    }

    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    for (int j = 0; j < nDraw; ++j)
    {
      //KS: This term only makes sense if we have at least 2 chains
      if(Nchains == 1)
      {
        BetweenChainVariance[j] = 0.;
        BetweenChainVarianceFolded[j] = 0.;
      }
      else
      {
        for (int m = 0; m < Nchains; ++m)
        {
          BetweenChainVariance[j] += ( Mean[m][j] - MeanGlobal[j])*( Mean[m][j] - MeanGlobal[j]);
          BetweenChainVarianceFolded[j] += ( MeanFolded[m][j] - MeanGlobalFolded[j])*( MeanFolded[m][j] - MeanGlobalFolded[j]);
        }
        BetweenChainVariance[j] = BetweenChainVariance[j]*Ntoys/(Nchains-1);
        BetweenChainVarianceFolded[j] = BetweenChainVarianceFolded[j]*Ntoys/(Nchains-1);
      }
    }

    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    for (int j = 0; j < nDraw; ++j)
    {
      MarginalPosteriorVariance[j] = (Ntoys-1) * StandardDeviationGlobal[j] /(Ntoys) + BetweenChainVariance[j]/Ntoys;
      MarginalPosteriorVarianceFolded[j] = (Ntoys-1) * StandardDeviationGlobalFolded[j] /(Ntoys) + BetweenChainVarianceFolded[j]/Ntoys;
    }

    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    //Finally calculate our estimator
    for (int j = 0; j < nDraw; ++j)
    {
      RHat[j] = sqrt(MarginalPosteriorVariance[j]/StandardDeviationGlobal[j]);
      RHatFolded[j] = sqrt(MarginalPosteriorVarianceFolded[j]/StandardDeviationGlobalFolded[j]);

      //KS: For flat params values can be crazy so cap at 0
      CapVariable(RHat[j], 0);
      CapVariable(RHatFolded[j], 0);
    }

    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    //KS: Additionally calculates effective step size which is an estimate of the sample size required to achieve the same level of precision if that sample was a simple random sample.
    for (int j = 0; j < nDraw; ++j)
    {
      if(Nchains > 1) EffectiveSampleSize[j] = Nchains * Ntoys * MarginalPosteriorVariance[j] / BetweenChainVariance[j];
      if(Nchains > 1) EffectiveSampleSizeFolded[j] = Nchains * Ntoys * MarginalPosteriorVarianceFolded[j] / BetweenChainVarianceFolded[j];

      //KS: For flat params values can be crazy so cap at 0
      CapVariable(EffectiveSampleSize[j], 0);
      CapVariable(EffectiveSampleSizeFolded[j], 0);
    }
  #ifdef MULTITHREAD
  } //End parallel region
  #endif

  clock.Stop();
  MACH3LOG_INFO("Finished calculating RHat, it took {:.2f}s to finish", clock.RealTime());
}


// *******************
void SaveResults() {
// *******************
  #pragma GCC diagnostic ignored "-Wfloat-conversion"

  std::string NameTemp = "";
  //KS: If we run over many many chains there is danger that name will be so absurdly long we run over system limit and job will be killed :(
  if(Nchains < 5)
  {
    for (int i = 0; i < Nchains; i++)
    {
      std::string temp = MCMCFile[i];

      while (temp.find(".root") != std::string::npos) {
        temp = temp.substr(0, temp.find(".root"));
      }

      NameTemp = NameTemp + temp + "_";
    }
  }
  else {
    NameTemp = std::to_string(Nchains) + "Chains" + "_";
  }
  NameTemp += "diag.root";

  TFile* DiagFile = new TFile(NameTemp.c_str(), "recreate");

  DiagFile->cd();

  TH1D *StandardDeviationGlobalPlot = new TH1D("StandardDeviationGlobalPlot", "StandardDeviationGlobalPlot", 200, 0, 2);
  TH1D *BetweenChainVariancePlot = new TH1D("BetweenChainVariancePlot", "BetweenChainVariancePlot", 200, 0, 2);
  TH1D *MarginalPosteriorVariancePlot = new TH1D("MarginalPosteriorVariancePlot", "MarginalPosteriorVariancePlot", 200, 0, 2);
  TH1D *RhatPlot = new TH1D("RhatPlot", "RhatPlot", 200, 0, 2);
  TH1D *EffectiveSampleSizePlot = new TH1D("EffectiveSampleSizePlot", "EffectiveSampleSizePlot", 400, 0, 10000);

  TH1D *StandardDeviationGlobalFoldedPlot = new TH1D("StandardDeviationGlobalFoldedPlot", "StandardDeviationGlobalFoldedPlot", 200, 0, 2);
  TH1D *BetweenChainVarianceFoldedPlot = new TH1D("BetweenChainVarianceFoldedPlot", "BetweenChainVarianceFoldedPlot", 200, 0, 2);
  TH1D *MarginalPosteriorVarianceFoldedPlot = new TH1D("MarginalPosteriorVarianceFoldedPlot", "MarginalPosteriorVarianceFoldedPlot", 200, 0, 2);
  TH1D *RhatFoldedPlot = new TH1D("RhatFoldedPlot", "RhatFoldedPlot", 200, 0, 2);
  TH1D *EffectiveSampleSizeFoldedPlot = new TH1D("EffectiveSampleSizeFoldedPlot", "EffectiveSampleSizeFoldedPlot", 400, 0, 10000);

  TH1D *RhatLogPlot = new TH1D("RhatLogPlot", "RhatLogPlot", 200, 0, 2);
  TH1D *RhatFoldedLogPlot = new TH1D("RhatFoldedLogPlot", "RhatFoldedLogPlot", 200, 0, 2);

  int Criterium = 0;
  int CiteriumFolded = 0;
  for(int j = 0; j < nDraw; j++)
  {
    //KS: Fill only valid parameters
    if(ValidPar[j])
    {
      StandardDeviationGlobalPlot->Fill(StandardDeviationGlobal[j]);
      BetweenChainVariancePlot->Fill(BetweenChainVariance[j]);
      MarginalPosteriorVariancePlot->Fill(MarginalPosteriorVariance[j]);
      RhatPlot->Fill(RHat[j]);
      EffectiveSampleSizePlot->Fill(EffectiveSampleSize[j]);
      if(RHat[j] > 1.1) Criterium++;


      StandardDeviationGlobalFoldedPlot->Fill(StandardDeviationGlobalFolded[j]);
      BetweenChainVarianceFoldedPlot->Fill(BetweenChainVarianceFolded[j]);
      MarginalPosteriorVarianceFoldedPlot->Fill(MarginalPosteriorVarianceFolded[j]);
      RhatFoldedPlot->Fill(RHatFolded[j]);
      EffectiveSampleSizeFoldedPlot->Fill(EffectiveSampleSizeFolded[j]);
      if(RHatFolded[j] > 1.1) CiteriumFolded++;
    }
    else
    {
      RhatLogPlot->Fill(RHat[j]);
      RhatFoldedLogPlot->Fill(RHatFolded[j]);
    }
  }
  //KS: We set criterium of 1.1 based on Gelman et al. (2003) Bayesian Data Analysis
  MACH3LOG_WARN("Number of parameters which has R hat greater than 1.1 is {}({:.2f}%) while for R hat folded {}({:.2f}%)", Criterium, 100*double(Criterium)/double(nDraw), CiteriumFolded, 100*double(CiteriumFolded)/double(nDraw));
  for(int j = 0; j < nDraw; j++)
  {
    if( (RHat[j] > 1.1 || RHatFolded[j] > 1.1) && ValidPar[j])
    {
      MACH3LOG_CRITICAL("Parameter {} has R hat higher than 1.1", BranchNames[j]);
    }
  }
  StandardDeviationGlobalPlot->Write();
  BetweenChainVariancePlot->Write();
  MarginalPosteriorVariancePlot->Write();
  RhatPlot->Write();
  EffectiveSampleSizePlot->Write();

  StandardDeviationGlobalFoldedPlot->Write();
  BetweenChainVarianceFoldedPlot->Write();
  MarginalPosteriorVarianceFoldedPlot->Write();
  RhatFoldedPlot->Write();
  EffectiveSampleSizeFoldedPlot->Write();

  RhatLogPlot->Write();
  RhatFoldedLogPlot->Write();

  //KS: Now we make fancy canvases, consider some function to have less copy pasting
  auto TempCanvas = std::make_unique<TCanvas>("Canvas", "Canvas", 1024, 1024);
  gStyle->SetOptStat(0);
  TempCanvas->SetGridx();
  TempCanvas->SetGridy();

  // Random line to write useful information to TLegend
  auto TempLine = std::make_unique<TLine>(0, 0, 0, 0);
  TempLine->SetLineColor(kBlack);

  RhatPlot->GetXaxis()->SetTitle("R hat");
  RhatPlot->SetLineColor(kRed);
  RhatPlot->SetFillColor(kRed);
  RhatFoldedPlot->SetLineColor(kBlue);
  RhatFoldedPlot->SetFillColor(kBlue);

  TLegend *Legend = new TLegend(0.55, 0.6, 0.9, 0.9);
  Legend->SetTextSize(0.04);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);

  Legend->AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", Ntoys, Nchains), "");
  Legend->AddEntry(RhatPlot, "Rhat Gelman 2013", "l");
  Legend->AddEntry(RhatFoldedPlot, "Rhat-Folded Gelman 2021", "l");

  RhatPlot->Draw();
  RhatFoldedPlot->Draw("same");
  Legend->Draw("same");
  TempCanvas->Write("Rhat");
  delete Legend;
  Legend = nullptr;

  //Now R hat for log L
  RhatLogPlot->GetXaxis()->SetTitle("R hat for LogL");
  RhatLogPlot->SetLineColor(kRed);
  RhatLogPlot->SetFillColor(kRed);
  RhatFoldedLogPlot->SetLineColor(kBlue);
  RhatFoldedLogPlot->SetFillColor(kBlue);

  Legend = new TLegend(0.55, 0.6, 0.9, 0.9);
  Legend->SetTextSize(0.04);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);

  Legend->AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", Ntoys, Nchains), "");
  Legend->AddEntry(RhatLogPlot, "Rhat Gelman 2013", "l");
  Legend->AddEntry(RhatFoldedLogPlot, "Rhat-Folded Gelman 2021", "l");

  RhatLogPlot->Draw();
  RhatFoldedLogPlot->Draw("same");
  Legend->Draw("same");
  TempCanvas->Write("RhatLog");
  delete Legend;
  Legend = nullptr;

  //Now canvas for effective sample size
  EffectiveSampleSizePlot->GetXaxis()->SetTitle("S_{eff, BDA2}");
  EffectiveSampleSizePlot->SetLineColor(kRed);
  EffectiveSampleSizeFoldedPlot->SetLineColor(kBlue);

  Legend = new TLegend(0.45, 0.6, 0.9, 0.9);
  Legend->SetTextSize(0.03);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);

  const double Mean1 = EffectiveSampleSizePlot->GetMean();
  const double RMS1 = EffectiveSampleSizePlot->GetRMS();
  const double Mean2 = EffectiveSampleSizeFoldedPlot->GetMean();
  const double RMS2 = EffectiveSampleSizeFoldedPlot->GetRMS();

  Legend->AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", Ntoys, Nchains), "");
  Legend->AddEntry(EffectiveSampleSizePlot, Form("S_{eff, BDA2} #mu = %.2f, #sigma = %.2f",Mean1 ,RMS1), "l");
  Legend->AddEntry(EffectiveSampleSizeFoldedPlot, Form("S_{eff, BDA2} Folded, #mu = %.2f, #sigma = %.2f",Mean2 ,RMS2), "l");

  EffectiveSampleSizePlot->Draw();
  EffectiveSampleSizeFoldedPlot->Draw("same");
  Legend->Draw("same");
  TempCanvas->Write("EffectiveSampleSize");

  //Fancy memory cleaning
  delete StandardDeviationGlobalPlot;
  delete BetweenChainVariancePlot;
  delete MarginalPosteriorVariancePlot;
  delete RhatPlot;
  delete EffectiveSampleSizePlot;

  delete StandardDeviationGlobalFoldedPlot;
  delete BetweenChainVarianceFoldedPlot;
  delete MarginalPosteriorVarianceFoldedPlot;
  delete RhatFoldedPlot;
  delete EffectiveSampleSizeFoldedPlot;

  delete Legend;

  delete RhatLogPlot;
  delete RhatFoldedLogPlot;

  DiagFile->Close();
  delete DiagFile;

  MACH3LOG_INFO("Finished and wrote results to {}", NameTemp);
}

// *******************
//KS: Pseudo destructor
void DestroyArrays() {
// *******************

  MACH3LOG_INFO("Killing all arrays");
  delete[] MeanGlobal;
  delete[] StandardDeviationGlobal;
  delete[] BetweenChainVariance;
  delete[] MarginalPosteriorVariance;
  delete[] RHat;
  delete[] EffectiveSampleSize;

  delete[] MeanGlobalFolded;
  delete[] StandardDeviationGlobalFolded;
  delete[] BetweenChainVarianceFolded;
  delete[] MarginalPosteriorVarianceFolded;
  delete[] RHatFolded;
  delete[] EffectiveSampleSizeFolded;

  for(int m = 0; m < Nchains; m++)
  {
    for(int i = 0; i < Ntoys; i++)
    {
      delete[] Draws[m][i];
      delete[] DrawsFolded[m][i];
    }
    delete[] Draws[m];
    delete[] Mean[m];
    delete[] StandardDeviation[m];

    delete[] DrawsFolded[m];
    delete[] MeanFolded[m];
    delete[] StandardDeviationFolded[m];
  }
  delete[] Draws;
  delete[] Mean;
  delete[] StandardDeviation;

  delete[] DrawsFolded;
  delete[] MedianArr;
  delete[] MeanFolded;
  delete[] StandardDeviationFolded;
}

// *******************
//calculate median
double CalcMedian(double arr[], const int size) {
// *******************
  std::sort(arr, arr+size);
  if (size % 2 != 0)
    return arr[size/2];
  return (arr[(size-1)/2] + arr[size/2])/2.0;
}

// *******************
//calculate median
void CapVariable(double var, const double cap) {
// *******************
  if(std::isnan(var) || !std::isfinite(var)) var = cap;
}
