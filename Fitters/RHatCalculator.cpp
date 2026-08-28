#include "RHatCalculator.h"
#include "MCMCProcessor.h"
#include <filesystem>
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

// ****************************
RHatCalculator::RHatCalculator(bool HighMemory, std::vector<std::string>& Inputs, int entries) {
// ****************************
  SetMaCh3LoggerFormat();

  HighMemoryMode = HighMemory;

  MCMCFile = Inputs;
  Nchains = static_cast<int>(MCMCFile.size());
  if(HighMemoryMode){
    Ntoys = entries;
  } else{
    NThin = entries;
  }
}

// ****************************
// The destructor
RHatCalculator::~RHatCalculator() {
// ****************************
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

// *******************
void RHatCalculator::RunDiagnostic() {
// *******************
  if(HighMemoryMode) {
    PrepareChains_HighMem();
  } else {
    PrepareChains();
  }

  InitialiseArrays();

  //KS: Main function
  if(HighMemoryMode) {
    CalcRhat_HighMem();
  } else {
    CalcRhat();
  }
  SaveResults();
}

// *******************
// Load chain and prepare toys
void RHatCalculator::PrepareChains_HighMem() {
// *******************
  auto rnd = std::make_unique<TRandom3>(0);

  MACH3LOG_INFO("Generating {}", Ntoys);

  TStopwatch clock;
  clock.Start();

  std::vector<unsigned int> BurnIn(Nchains);
  std::vector<unsigned int> nEntries(Nchains);
  std::vector<int> nBranches(Nchains);
  std::vector<unsigned int> step(Nchains);

  Draws.resize(Nchains);
  DrawsFolded.resize(Nchains);

  std::unique_ptr<MCMCProcessor> Processor;

  // Open the Chain
  //It is tempting to multithread here but unfortunately, ROOT files are not thread safe :(
  for (int m = 0; m < Nchains; m++)
  {
    MACH3LOG_INFO("On file: {}", MCMCFile[m].c_str());
    if (!std::filesystem::exists(MCMCFile[m])) {
      MACH3LOG_ERROR("File: {}, doesn't exist", MCMCFile[m]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    TChain* Chain = new TChain("posteriors");
    Chain->Add(MCMCFile[m].c_str());
    nEntries[m] = static_cast<unsigned int>(Chain->GetEntries());

    // Set the step cut to be 20%
    BurnIn[m] = nEntries[m]/5;

    // Get the list of branches
    TObjArray* brlis = Chain->GetListOfBranches();

    // Get the number of branches
    nBranches[m] = brlis->GetEntries();

    // Set all the branches to off
    Chain->SetBranchStatus("*", false);

    std::vector<TString> SampleLLH;
    std::vector<TString> SystLLH;

    if(m == 0) {
      Processor = std::make_unique<MCMCProcessor>(MCMCFile[m]);
      Processor->Initialise();
      nDraw = Processor->GetNParams();
      SampleLLH = Processor->GetSampleBranchNames();
      SystLLH = Processor->GetSystBranchNames();
      nDraw = nDraw + static_cast<int>(SampleLLH.size()) + static_cast<int>(SystLLH.size());
      BranchNames.reserve(nDraw);
    }

    // Set all the branches to off
    Chain->SetBranchStatus("*", false);

    // Loop over the number of branches
    // Find the name and how many of each systematic we have
    for (int i = 0; i < Processor->GetNParams(); i++)
    {
      TString bname = Processor->GetBranchNames()[i];
      //KS: Save branch name only for one chain, we assume all chains have the same branches, otherwise this doesn't make sense either way
      if(m == 0) {
        BranchNames.push_back(bname);
        //KS: We calculate R Hat also for LogL, just in case, however we plot them separately
        ValidPar.push_back(true);
      }
      Chain->SetBranchStatus(bname, true);
      // Get the TBranch and its name
      if (!Chain->GetBranch(bname)) {
        MACH3LOG_ERROR("Branch '{}' does not exist in the TChain", bname.Data());
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      MACH3LOG_DEBUG("{}", bname);
    }

    for (size_t i = 0; i < SampleLLH.size(); ++i) {
      TString bname = SampleLLH[i];
      Chain->SetBranchStatus(bname, true);
      if(m == 0) {
        BranchNames.push_back(bname);
        //KS: We calculate R Hat also for LogL, just in case, however we plot them separately
        ValidPar.push_back(false);
      }
    }
    for (size_t i = 0; i < SystLLH.size(); ++i) {
      TString bname = SystLLH[i];
      Chain->SetBranchStatus(bname, true);
      if(m == 0) {
        BranchNames.push_back(bname);
        //KS: We calculate R Hat also for LogL, just in case, however we plot them separately
        ValidPar.push_back(false);
      }
    }
    // Read in the step
    Chain->SetBranchStatus("step", true);
    Chain->SetBranchAddress("step", &step[m]);
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
    Draws[m].resize(Ntoys);
    DrawsFolded[m].resize(Ntoys);
    for(int i = 0; i < Ntoys; i++)
    {
      Draws[m][i].resize(nDraw, 0.0);
      DrawsFolded[m][i].resize(nDraw, 0.0);
    }

    // MJR: array to hold branch values; SetBranchAddress in every step is very
    //      expensive, so doing it once only here saves time
    std::vector<double> branch_values(nDraw, 0.0);
    for (int j = 0; j < nDraw; ++j) {
      Chain->SetBranchAddress(BranchNames[j].Data(), &branch_values[j]);
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
      if (step[m] < BurnIn[m]) {
        i--;
        continue;
      }

      // Output some info for the user
      if (Ntoys > 10 && i % (Ntoys/10) == 0) {
        M3::Utils::PrintProgressBar(i+m*Ntoys, static_cast<Long64_t>(Ntoys)*Nchains);
        MACH3LOG_DEBUG("Getting random entry {}", entry);
      }

      // Set the branch addresses for params
      for (int j = 0; j < nDraw; ++j) {
        Draws[m][i][j] = branch_values[j];
      }
    }//end loop over toys

    //TN: There, we now don't need to keep the chain in memory anymore
    delete Chain;
  }

  //KS: Now prepare folded draws, quoting Gelman
  //"We propose to report the maximum of rank normalized split-Rb and rank normalized folded-split-Rb for each parameter"
  MedianArr.resize(nDraw);
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
// Load chain and prepare toys
void RHatCalculator::PrepareChains() {
// *******************
  TStopwatch clock;
  clock.Start();

  Ntoys_requested.resize(Nchains);
  Ntoys_filled.resize(Nchains);
  TotToys = 0;
  std::vector<unsigned int> BurnIn(Nchains);
  std::vector<unsigned int> nEntries(Nchains);
  std::vector<int> nBranches(Nchains);
  std::vector<unsigned int> step(Nchains);

  S1_chain.resize(Nchains);
  S2_chain.resize(Nchains);

  std::unique_ptr<MCMCProcessor> Processor;

  // Open the Chain
  //It is tempting to multithread here but unfortunately, ROOT files are not thread safe :(
  for (int m = 0; m < Nchains; m++)
  {
    MACH3LOG_INFO("On file: {}", MCMCFile[m].c_str());
    if (!std::filesystem::exists(MCMCFile[m])) {
      MACH3LOG_ERROR("File: {}, doesn't exist", MCMCFile[m]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
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

    // Set all the branches to off
    Chain->SetBranchStatus("*", false);

    std::vector<TString> SampleLLH;
    std::vector<TString> SystLLH;

    if(m == 0) {
      Processor = std::make_unique<MCMCProcessor>(MCMCFile[m]);
      Processor->Initialise();
      nDraw = Processor->GetNParams();
      SampleLLH = Processor->GetSampleBranchNames();
      SystLLH = Processor->GetSystBranchNames();
      nDraw = nDraw + static_cast<int>(SampleLLH.size()) + static_cast<int>(SystLLH.size());
      BranchNames.reserve(nDraw);
    }

    // Set all the branches to off
    Chain->SetBranchStatus("*", false);

    // Loop over the number of branches
    // Find the name and how many of each systematic we have
    for (int i = 0; i < Processor->GetNParams(); i++)
    {
      TString bname = Processor->GetBranchNames()[i];
      //KS: Save branch name only for one chain, we assume all chains have the same branches, otherwise this doesn't make sense either way
      if(m == 0) {
        BranchNames.push_back(bname);
        //KS: We calculate R Hat also for LogL, just in case, however we plot them separately
        ValidPar.push_back(true);
      }
      Chain->SetBranchStatus(bname, true);

      // Get the TBranch and its name
      if (!Chain->GetBranch(bname)) {
        MACH3LOG_ERROR("Branch '{}' does not exist in the TChain", bname.Data());
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      MACH3LOG_DEBUG("{}", bname);
    }

    for (size_t i = 0; i < SampleLLH.size(); ++i) {
      TString bname = SampleLLH[i];
      Chain->SetBranchStatus(bname, true);
      if(m == 0) {
        BranchNames.push_back(bname);
        //KS: We calculate R Hat also for LogL, just in case, however we plot them separately
        ValidPar.push_back(false);
      }
    }
    for (size_t i = 0; i < SystLLH.size(); ++i) {
      TString bname = SystLLH[i];
      Chain->SetBranchStatus(bname, true);
      if(m == 0) {
        BranchNames.push_back(bname);
        //KS: We calculate R Hat also for LogL, just in case, however we plot them separately
        ValidPar.push_back(false);
      }
    }
    // Read in the step
    Chain->SetBranchStatus("step", true);
    Chain->SetBranchAddress("step", &step[m]);
    // MJR: Initialize quantities needed for calculating RHat
    S1_chain[m].resize(nDraw, 0);
    S2_chain[m].resize(nDraw, 0);
    if (m == 0)
    {
      S1_global.resize(nDraw, 0);
      S2_global.resize(nDraw, 0);
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
    std::vector<double> branch_values(nDraw, 0.0);
    for (int id = 0; id < nDraw; ++id) {
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
        M3::Utils::PrintProgressBar(i+m*Ntoys_requested[m], static_cast<Long64_t>(Ntoys_requested[m])*Nchains);
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
    MACH3LOG_INFO("Finished loading chain {}!", m);
  }

  clock.Stop();
  MACH3LOG_INFO("Finished calculating Toys, it took {:.2f}s to finish", clock.RealTime());
}

// *******************
// Create all arrays we are going to use later
void RHatCalculator::InitialiseArrays() {
// *******************
  MACH3LOG_INFO("Starting {}", __func__);
  Mean.resize(Nchains);
  StandardDeviation.resize(Nchains);

  for (int m = 0; m < Nchains; ++m) {
    Mean[m].resize(nDraw, 0);
    StandardDeviation[m].resize(nDraw, 0);
  }

  MeanGlobal.resize(nDraw, 0);
  StandardDeviationGlobal.resize(nDraw, 0);
  BetweenChainVariance.resize(nDraw, 0);

  MarginalPosteriorVariance.resize(nDraw, 0);
  RHat.resize(nDraw, 0);
  EffectiveSampleSize.resize(nDraw, 0);

  if(HighMemoryMode){
    MeanFolded.resize(Nchains);
    StandardDeviationFolded.resize(Nchains);

    for (int m = 0; m < Nchains; ++m) {
      MeanFolded[m].resize(nDraw, 0);
      StandardDeviationFolded[m].resize(nDraw, 0);
    }

    MeanGlobalFolded.resize(nDraw, 0);
    StandardDeviationGlobalFolded.resize(nDraw, 0);
    BetweenChainVarianceFolded.resize(nDraw, 0);

    MarginalPosteriorVarianceFolded.resize(nDraw, 0);
    RHatFolded.resize(nDraw, 0);
    EffectiveSampleSizeFolded.resize(nDraw, 0);
  }
}

// *******************
//KS: Based on Gelman et. al. arXiv:1903.08008v5
void RHatCalculator::CalcRhat() {
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
//KS: Based on Gelman et. al. arXiv:1903.08008v5
// Probably most of it could be moved cleverly to MCMC Processor, keep it separate for now
void RHatCalculator::CalcRhat_HighMem() {
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
void RHatCalculator::SaveResults() {
// *******************
  #pragma GCC diagnostic push
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

  TH1D *StandardDeviationGlobalFoldedPlot = nullptr;
  TH1D *BetweenChainVarianceFoldedPlot = nullptr;
  TH1D *MarginalPosteriorVarianceFoldedPlot = nullptr;
  TH1D *RhatFoldedPlot = nullptr;
  TH1D *EffectiveSampleSizeFoldedPlot = nullptr;
  TH1D *RhatFoldedLogPlot = nullptr;

  if (HighMemoryMode)
  {
    StandardDeviationGlobalFoldedPlot = new TH1D("StandardDeviationGlobalFoldedPlot", "StandardDeviationGlobalFoldedPlot", nDraw, 0, nDraw);
    BetweenChainVarianceFoldedPlot = new TH1D("BetweenChainVarianceFoldedPlot", "BetweenChainVarianceFoldedPlot", nDraw, 0, nDraw);
    MarginalPosteriorVarianceFoldedPlot = new TH1D("MarginalPosteriorVarianceFoldedPlot", "MarginalPosteriorVarianceFoldedPlot", nDraw, 0, nDraw);
    RhatFoldedPlot = new TH1D("RhatFoldedPlot", "RhatFoldedPlot", 200, 0, 2);
    EffectiveSampleSizeFoldedPlot = new TH1D("EffectiveSampleSizeFoldedPlot", "EffectiveSampleSizeFoldedPlot", 400, 0, 10000);
    RhatFoldedLogPlot = new TH1D("RhatFoldedLogPlot", "RhatFoldedLogPlot", 200, 0, 2);
  }

  int Criterium = 0;
  int CiteriumFolded = 0;
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
      if(HighMemoryMode) {
        StandardDeviationGlobalFoldedPlot->Fill(j,StandardDeviationGlobalFolded[j]);
        BetweenChainVarianceFoldedPlot->Fill(j,BetweenChainVarianceFolded[j]);
        MarginalPosteriorVarianceFoldedPlot->Fill(j,MarginalPosteriorVarianceFolded[j]);
        RhatFoldedPlot->Fill(RHatFolded[j]);
        EffectiveSampleSizeFoldedPlot->Fill(EffectiveSampleSizeFolded[j]);
        if(RHatFolded[j] > 1.1) CiteriumFolded++;
      }
    }
    else
    {
      RhatLogPlot->Fill(RHat[j]);
      if(HighMemoryMode) RhatFoldedLogPlot->Fill(RHatFolded[j]);
    }
  }

  if (HighMemoryMode)
  {
    //KS: We set criterium of 1.1 based on Gelman et al. (2003) Bayesian Data Analysis
    MACH3LOG_WARN("Number of parameters which has R hat greater than 1.1 is {}({:.2f}%) while for R hat folded {}({:.2f}%)",
                   Criterium, 100*double(Criterium)/double(nDraw), CiteriumFolded, 100*double(CiteriumFolded)/double(nDraw));
    for(int j = 0; j < nDraw; j++)
    {
      if( (RHat[j] > 1.1 || RHatFolded[j] > 1.1) && ValidPar[j])
      {
        MACH3LOG_CRITICAL("Parameter {} has R hat higher than 1.1", BranchNames[j]);
      }
    }
  }
  else
  {
    //KS: We set criterium of 1.1 based on Gelman et al. (2003) Bayesian Data Analysis
    MACH3LOG_WARN("Number of parameters which has R hat greater than 1.1 is {}({:.2f}%)", Criterium, 100*double(Criterium)/double(nDraw));
    for(int j = 0; j < nDraw; j++)
    {
      if( (RHat[j] > 1.1) && ValidPar[j])
      {
        MACH3LOG_CRITICAL("Parameter {} has R hat higher than 1.1", BranchNames[j]);
      }
    }
  }

  StandardDeviationGlobalPlot->Write();
  BetweenChainVariancePlot->Write();
  MarginalPosteriorVariancePlot->Write();
  RhatPlot->Write();
  EffectiveSampleSizePlot->Write();

  RhatLogPlot->Write();

  if (HighMemoryMode)
  {
    StandardDeviationGlobalFoldedPlot->Write();
    BetweenChainVarianceFoldedPlot->Write();
    MarginalPosteriorVarianceFoldedPlot->Write();
    RhatFoldedPlot->Write();
    EffectiveSampleSizeFoldedPlot->Write();

    RhatFoldedLogPlot->Write();
  }

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

  if(HighMemoryMode){
    RhatFoldedPlot->SetLineColor(kBlue);
    RhatFoldedPlot->SetFillColor(kBlue);
  }

  TLegend Legend(0.55, 0.6, 0.9, 0.9);
  Legend.SetTextSize(0.04);
  Legend.SetFillColor(0);
  Legend.SetFillStyle(0);
  Legend.SetLineWidth(0);
  Legend.SetLineColor(0);

  Legend.AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", TotToys, Nchains), "");
  Legend.AddEntry(RhatPlot, "Rhat Gelman 2013", "l");
  if(HighMemoryMode) Legend.AddEntry(RhatFoldedPlot, "Rhat-Folded Gelman 2021", "l");

  RhatPlot->Draw();
  Legend.Draw("same");
  if(HighMemoryMode) RhatFoldedPlot->Draw("same");
  TempCanvas->Write("Rhat");

  //Now R hat for log L
  RhatLogPlot->GetXaxis()->SetTitle("R hat for LogL");
  RhatLogPlot->SetLineColor(kRed);
  RhatLogPlot->SetFillColor(kRed);
  if(HighMemoryMode) RhatFoldedLogPlot->SetLineColor(kBlue);
  if(HighMemoryMode) RhatFoldedLogPlot->SetFillColor(kBlue);

  TLegend LegendFolded(0.55, 0.6, 0.9, 0.9);

  LegendFolded.SetTextSize(0.04);
  LegendFolded.SetFillColor(0);
  LegendFolded.SetFillStyle(0);
  LegendFolded.SetLineWidth(0);
  LegendFolded.SetLineColor(0);

  LegendFolded.AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", TotToys, Nchains), "");
  LegendFolded.AddEntry(RhatLogPlot, "Rhat Gelman 2013", "l");
  if(HighMemoryMode) LegendFolded.AddEntry(RhatFoldedLogPlot, "Rhat-Folded Gelman 2021", "l");

  RhatLogPlot->Draw();
  LegendFolded.Draw("same");
  TempCanvas->Write("RhatLog");

  //Now canvas for effective sample size
  EffectiveSampleSizePlot->GetXaxis()->SetTitle("S_{eff, BDA2}");
  EffectiveSampleSizePlot->SetLineColor(kRed);
  if(HighMemoryMode) EffectiveSampleSizeFoldedPlot->SetLineColor(kBlue);

  TLegend LegendESS(0.45, 0.6, 0.9, 0.9);
  LegendESS.SetTextSize(0.03);
  LegendESS.SetFillColor(0);
  LegendESS.SetFillStyle(0);
  LegendESS.SetLineWidth(0);
  LegendESS.SetLineColor(0);

  const double Mean1 = EffectiveSampleSizePlot->GetMean();
  const double RMS1 = EffectiveSampleSizePlot->GetRMS();

  LegendESS.AddEntry(TempLine.get(), Form("Number of throws=%.0i, Number of chains=%.1i", TotToys, Nchains), "");
  LegendESS.AddEntry(EffectiveSampleSizePlot, Form("S_{eff, BDA2} #mu = %.2f, #sigma = %.2f",Mean1 ,RMS1), "l");
  if (HighMemoryMode)
  {
    const double Mean2 = EffectiveSampleSizeFoldedPlot->GetMean();
    const double RMS2  = EffectiveSampleSizeFoldedPlot->GetRMS();
    LegendESS.AddEntry(EffectiveSampleSizeFoldedPlot, Form("S_{eff, BDA2} Folded, #mu = %.2f, #sigma = %.2f", Mean2, RMS2), "l");
  }
  EffectiveSampleSizePlot->Draw();
  LegendESS.Draw("same");
  if(HighMemoryMode) EffectiveSampleSizeFoldedPlot->Draw("same");
  TempCanvas->Write("EffectiveSampleSize");


  //Fancy memory cleaning
  delete StandardDeviationGlobalPlot;
  delete BetweenChainVariancePlot;
  delete MarginalPosteriorVariancePlot;
  delete RhatPlot;
  delete EffectiveSampleSizePlot;

  delete RhatLogPlot;

  if (HighMemoryMode)
  {
    delete StandardDeviationGlobalFoldedPlot;
    delete BetweenChainVarianceFoldedPlot;
    delete MarginalPosteriorVarianceFoldedPlot;
    delete RhatFoldedPlot;
    delete EffectiveSampleSizeFoldedPlot;

    delete RhatFoldedLogPlot;
  }

  DiagFile->Close();
  delete DiagFile;

  MACH3LOG_INFO("Finished and wrote results to {}", NameTemp);
  #pragma GCC diagnostic pop
}
