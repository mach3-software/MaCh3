// MaCh3 includes
#include "Manager/Manager.h"
#include "Samples/SampleStructs.h"
#include "Samples/HistogramUtils.h"

_MaCh3_Safe_Include_Start_ //{
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
_MaCh3_Safe_Include_End_ //}

/// @file RHat.cpp
/// @brief This executable calculates the \f$ \hat{R} \f$ estimator for Markov Chain Monte Carlo (MCMC) convergence.
///
/// KS: This exe is meant to calculate the \f$ \hat{R} \f$ estimator. For a well-converged chain, this distribution
/// should be centered at one. The \f$ \hat{R} \f$ statistic is used to assess the convergence of MCMC simulations
/// and helps determine whether the chains have reached a stable distribution.
///
/// @cite gelman2019.
///
/// MJR: Update -- Improved memory usage so that whole chains can be quickly loaded without requiring copious amounts
///      of RAM. This comes at the cost of not being able to calculate the Folded RHat since finding the median
///      requires the loading of full chains at a time. The method has been validated to give identical results to the
///      "High Memory" (original) version at a fraction of the runtime and resources.
///
///      The input format is also slightly altered; since we can now load entire chains, there's less need to
///      specify how many toys are desired for a sub-sample, so the Ntoys input has been removed.
///
/// @author Kamil Skwarczynski
/// @author Michael Reh

// *******************
int* Ntoys_requested;
int* Ntoys_filled;
int TotToys;
unsigned int NThin;
int Nchains;
int nDraw;

std::vector<TString> BranchNames;
std::vector<std::string> MCMCFile;
std::vector<bool> ValidPar;

double* S1_global; // Sum_i^N x_i   | total
double* S2_global; // Sum_i^N x_i^2 | total
double** S1_chain; // Sum_i^N x_i   | for each chain
double** S2_chain; // Sum_i^N x_i^2 | for each chain

double** Mean;
double** StandardDeviation;

double* MeanGlobal;
double* StandardDeviationGlobal;

double* BetweenChainVariance;
double* MarginalPosteriorVariance;
double* RHat;
double* EffectiveSampleSize;

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

  Mean = nullptr;
  StandardDeviation = nullptr;

  MeanGlobal = nullptr;
  StandardDeviationGlobal = nullptr;

  BetweenChainVariance = nullptr;
  MarginalPosteriorVariance = nullptr;
  RHat = nullptr;
  EffectiveSampleSize = nullptr;

  Nchains = 0;

  if (argc < 2)
  {
    MACH3LOG_ERROR("Wrong arguments");
    MACH3LOG_ERROR("./RHat NThin MCMCchain_1.root MCMCchain_2.root MCMCchain_3.root ... [how many you like]");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  NThin = atoi(argv[1]);

  //KS Gelman suggests to diagnose on more than one chain
  for (int i = 2; i < argc; i++)
  {
    MCMCFile.push_back(std::string(argv[i]));
    MACH3LOG_INFO("Adding file: {}", MCMCFile.back());
    Nchains++;
  }

  if(Nchains == 1)
  {
    MACH3LOG_WARN("Gelman is going to be sad :(. He suggested you should use more than one chain (at least 4). Code works fine for one chain, however, estimator might be biased.");
    MACH3LOG_WARN("Multiple chains are more likely to reveal multimodality and poor adaptation or mixing:");
  }
  MACH3LOG_INFO("Diagnosing {} chains", Nchains);

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

  TStopwatch clock;
  clock.Start();

  Ntoys_requested = new int[Nchains]();
  Ntoys_filled = new int[Nchains]();
  TotToys = 0;
  std::vector<unsigned int> BurnIn(Nchains);
  std::vector<unsigned int> nEntries(Nchains);
  std::vector<int> nBranches(Nchains);
  std::vector<unsigned int> step(Nchains);

  S1_chain = new double*[Nchains]();
  S2_chain = new double*[Nchains]();

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

    nEntries[m] = static_cast<unsigned int>(Chain->GetEntries());
    Ntoys_requested[m] = nEntries[m]/NThin;
    Ntoys_filled[m] = 0;

    MACH3LOG_INFO("On file: {}", MCMCFile[m].c_str());
    MACH3LOG_INFO("Generating {} Toys", Ntoys_requested[m]);

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

    // MJR: Initialize quantities needed for calculating RHat
    S1_chain[m] = new double[nDraw]();
    S2_chain[m] = new double[nDraw]();
    if (m == 0)
    {
      S1_global = new double[nDraw]();
      S2_global = new double[nDraw]();
    }
    for (int id = 0; id < nDraw; ++id)
    {
      S1_chain[m][id] = 0.0;
      S2_chain[m][id] = 0.0;
      if (m == 0)
      {
        S1_global[id] = 0.0;
        S2_global[id] = 0.0;
      }
    }

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

    // MJR: Create an array to hold branch values. Resetting branch addresses
    //      for every step is very expensive.
    double* branch_values = new double[nDraw]();
    for (int id = 0; id < nDraw; ++id)
    {
      Chain->SetBranchAddress(BranchNames[id].Data(), &branch_values[id]);
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

    MACH3LOG_INFO("Loading chain {} / {}...", m, Nchains);
    for (int i = 0; i < Ntoys_requested[m]; i++)
    {
      // This is here as a placeholder in case we want to do some thinning later
      int entry = i*NThin;

      Chain->GetEntry(entry);

      // If we have combined chains by hadd need to check the step in the chain
      // Note, entry is not necessarily the same as the step due to merged ROOT files, so can't choose an entry in the range BurnIn - nEntries :(
      if (step[m] < BurnIn[m])
      {
        continue;
      }

      // Output some info for the user
      if (Ntoys_requested[m] > 10 && i % (Ntoys_requested[m]/10) == 0) {
        MaCh3Utils::PrintProgressBar(i+m*Ntoys_requested[m], static_cast<Long64_t>(Ntoys_requested[m])*Nchains);
        MACH3LOG_DEBUG("Getting random entry {}", entry);
      }

      // MJR: Fill running quantities instead of loading everything into RAM.
      //      This is where we save big on both memory and time (resetting
      //      branch addresses and calling GetEntry() again here is very slow).
      for (int j = 0; j < nDraw; ++j)
      {
        S1_global[j] += branch_values[j];
        S2_global[j] += branch_values[j]*branch_values[j];
        S1_chain[m][j] += branch_values[j];
        S2_chain[m][j] += branch_values[j]*branch_values[j];
      }

      // Increment counters
      Ntoys_filled[m]++;
      TotToys++;

    }//end loop over toys

    //TN: There, we now don't need to keep the chain in memory anymore
    delete Chain;
    delete[] branch_values;
    MACH3LOG_INFO("Finished loading chain {}!", m);
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

  for (int m = 0; m < Nchains; ++m)
  {
    Mean[m] = new double[nDraw]();
    StandardDeviation[m] = new double[nDraw]();

    for (int j = 0; j < nDraw; ++j)
    {
      Mean[m][j] = 0.;
      StandardDeviation[m][j] = 0.;

      if(m == 0)
      {
        MeanGlobal[j] = 0.;
        StandardDeviationGlobal[j] = 0.;
        BetweenChainVariance[j] = 0.;
        MarginalPosteriorVariance[j] = 0.;
        RHat[j] = 0.;
        EffectiveSampleSize[j] = 0.;
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
    #pragma omp for
    #endif
    //KS: loop over chains and draws are independent so might as well collapse for sweet cache hits
    //Calculate the mean for each parameter within each considered chain
    // MJR: Calculate using running totals to massively save on time and memory
    for (int m = 0; m < Nchains; ++m)
    {
      for (int j = 0; j < nDraw; ++j)
      {
        Mean[m][j] = S1_chain[m][j] / static_cast<double>(Ntoys_filled[m]);
        StandardDeviation[m][j] = S2_chain[m][j]/static_cast<double>(Ntoys_filled[m]) - Mean[m][j]*Mean[m][j];
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
        StandardDeviationGlobal[j] += StandardDeviation[m][j];
      }
      MeanGlobal[j] = S1_global[j] / static_cast<double>(TotToys);
      StandardDeviationGlobal[j] = StandardDeviationGlobal[j] / static_cast<double>(Nchains);
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
      }
      else
      {
        for (int m = 0; m < Nchains; ++m)
        {
          BetweenChainVariance[j] += ( Mean[m][j] - MeanGlobal[j])*( Mean[m][j] - MeanGlobal[j]) * Ntoys_filled[m];
        }
        BetweenChainVariance[j] = BetweenChainVariance[j]/(Nchains-1);
      }
    }

    int avgNtoys = TotToys/Nchains;
    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    for (int j = 0; j < nDraw; ++j)
    {
      MarginalPosteriorVariance[j] = (avgNtoys-1) * StandardDeviationGlobal[j] / (avgNtoys) + BetweenChainVariance[j]/avgNtoys;
    }

    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    //Finally calculate our estimator
    for (int j = 0; j < nDraw; ++j)
    {
      RHat[j] = sqrt(MarginalPosteriorVariance[j]/StandardDeviationGlobal[j]);

      //KS: For flat params values can be crazy so cap at 0
      CapVariable(RHat[j], 0);
    }

    #ifdef MULTITHREAD
    #pragma omp for
    #endif
    //KS: Additionally calculates effective step size which is an estimate of the sample size required to achieve the same level of precision if that sample was a simple random sample.
    for (int j = 0; j < nDraw; ++j)
    {
      if(Nchains > 1) EffectiveSampleSize[j] = TotToys * MarginalPosteriorVariance[j] / BetweenChainVariance[j];

      //KS: For flat params values can be crazy so cap at 0
      CapVariable(EffectiveSampleSize[j], 0);
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
      // Strip directory path
      const auto slash = temp.find_last_of("/\\");
      if (slash != std::string::npos) {
        temp = temp.substr(slash + 1);
      }

      NameTemp = NameTemp + temp + "_";
    }
  }
  else {
    NameTemp = std::to_string(Nchains) + "Chains" + "_";
  }
  NameTemp += "diag.root";

  TFile *DiagFile = M3::Open(NameTemp, "recreate", __FILE__, __LINE__);
  DiagFile->cd();

  TH1D *StandardDeviationGlobalPlot = new TH1D("StandardDeviationGlobalPlot", "StandardDeviationGlobalPlot", nDraw, 0, nDraw);
  TH1D *BetweenChainVariancePlot = new TH1D("BetweenChainVariancePlot", "BetweenChainVariancePlot", nDraw, 0, nDraw);
  TH1D *MarginalPosteriorVariancePlot = new TH1D("MarginalPosteriorVariancePlot", "MarginalPosteriorVariancePlot", nDraw, 0, nDraw);
  TH1D *RhatPlot = new TH1D("RhatPlot", "RhatPlot", 200, 0, 2);
  TH1D *EffectiveSampleSizePlot = new TH1D("EffectiveSampleSizePlot", "EffectiveSampleSizePlot", 400, 0, 10000);

  TH1D *RhatLogPlot = new TH1D("RhatLogPlot", "RhatLogPlot", 200, 0, 2);

  int Criterium = 0;
  for(int j = 0; j < nDraw; j++)
  {
    //KS: Fill only valid parameters
    if(ValidPar[j])
    {
      StandardDeviationGlobalPlot->Fill(j,StandardDeviationGlobal[j]);
      BetweenChainVariancePlot->Fill(j,BetweenChainVariance[j]);
      MarginalPosteriorVariancePlot->Fill(j,MarginalPosteriorVariance[j]);
      RhatPlot->Fill(RHat[j]);
      EffectiveSampleSizePlot->Fill(EffectiveSampleSize[j]);
      if(RHat[j] > 1.1) Criterium++;
    }
    else
    {
      RhatLogPlot->Fill(RHat[j]);
    }
  }
  //KS: We set criterium of 1.1 based on Gelman et al. (2003) Bayesian Data Analysis
  MACH3LOG_WARN("Number of parameters which has R hat greater than 1.1 is {}({:.2f}%)", Criterium, 100*double(Criterium)/double(nDraw));
  for(int j = 0; j < nDraw; j++)
  {
    if( (RHat[j] > 1.1) && ValidPar[j])
    {
      MACH3LOG_CRITICAL("Parameter {} has R hat higher than 1.1", BranchNames[j]);
    }
  }
  StandardDeviationGlobalPlot->Write();
  BetweenChainVariancePlot->Write();
  MarginalPosteriorVariancePlot->Write();
  RhatPlot->Write();
  EffectiveSampleSizePlot->Write();

  RhatLogPlot->Write();

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

  TLegend *Legend = new TLegend(0.55, 0.6, 0.9, 0.9);
  Legend->SetTextSize(0.04);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);

  Legend->AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", TotToys, Nchains), "");
  Legend->AddEntry(RhatPlot, "Rhat Gelman 2013", "l");

  RhatPlot->Draw();
  Legend->Draw("same");
  TempCanvas->Write("Rhat");
  delete Legend;
  Legend = nullptr;

  //Now R hat for log L
  RhatLogPlot->GetXaxis()->SetTitle("R hat for LogL");
  RhatLogPlot->SetLineColor(kRed);
  RhatLogPlot->SetFillColor(kRed);

  Legend = new TLegend(0.55, 0.6, 0.9, 0.9);
  Legend->SetTextSize(0.04);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);

  Legend->AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", TotToys, Nchains), "");
  Legend->AddEntry(RhatLogPlot, "Rhat Gelman 2013", "l");

  RhatLogPlot->Draw();
  Legend->Draw("same");
  TempCanvas->Write("RhatLog");
  delete Legend;
  Legend = nullptr;

  //Now canvas for effective sample size
  EffectiveSampleSizePlot->GetXaxis()->SetTitle("S_{eff, BDA2}");
  EffectiveSampleSizePlot->SetLineColor(kRed);

  Legend = new TLegend(0.45, 0.6, 0.9, 0.9);
  Legend->SetTextSize(0.03);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);

  const double Mean1 = EffectiveSampleSizePlot->GetMean();
  const double RMS1 = EffectiveSampleSizePlot->GetRMS();

  Legend->AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", TotToys, Nchains), "");
  Legend->AddEntry(EffectiveSampleSizePlot, Form("S_{eff, BDA2} #mu = %.2f, #sigma = %.2f",Mean1 ,RMS1), "l");

  EffectiveSampleSizePlot->Draw();
  Legend->Draw("same");
  TempCanvas->Write("EffectiveSampleSize");

  //Fancy memory cleaning
  delete StandardDeviationGlobalPlot;
  delete BetweenChainVariancePlot;
  delete MarginalPosteriorVariancePlot;
  delete RhatPlot;
  delete EffectiveSampleSizePlot;

  delete Legend;

  delete RhatLogPlot;

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

  for(int m = 0; m < Nchains; m++)
  {
    delete[] Mean[m];
    delete[] StandardDeviation[m];
    delete[] S1_chain[m];
    delete[] S2_chain[m];
  }
  delete[] Mean;
  delete[] StandardDeviation;
  delete[] S1_chain;
  delete[] S2_chain;
  delete[] S1_global;
  delete[] S2_global;

  delete[] Ntoys_requested;
  delete[] Ntoys_filled;
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
