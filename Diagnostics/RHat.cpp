// C++ includes
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>

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

#ifdef MULTITHREAD
#include "omp.h"
#endif

//KS: This exe is meant to calculate R hat estimator. For a well converged this distribution should be centred at one. 
//Based on Gelman et. al. arXiv:1903.08008v5

// *******************
int Ntoys;
int Nchains;

TChain *Chain;
int nDraw;

bool VERBOSE;
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
    Draws = NULL;
    Mean = NULL;
    StandardDeviation = NULL;
    
    MeanGlobal = NULL;
    StandardDeviationGlobal = NULL;
     
    BetweenChainVariance = NULL;
    MarginalPosteriorVariance = NULL;
    RHat = NULL;
    EffectiveSampleSize = NULL;
    
    DrawsFolded = NULL;
    MedianArr = NULL;
    MeanFolded = NULL;
    StandardDeviationFolded = NULL;
    
    MeanGlobalFolded = NULL;
    StandardDeviationGlobalFolded = NULL;
     
    BetweenChainVarianceFolded = NULL;
    MarginalPosteriorVarianceFolded = NULL;
    RHatFolded = NULL;
    EffectiveSampleSizeFolded = NULL;
    
    VERBOSE = false;
    Nchains = 0;
    
    if (argc == 1 && argc == 2) 
    {
        std::cerr << "Wrong arguments" << std::endl;
        std::cerr << "./RHat Ntoys MCMCchain_1.root MCMCchain_2.root MCMCchain_3.root ... [how many you like]" << std::endl;
        exit(-1);
    }
    
    Ntoys = atoi(argv[1]);
    
    //KS Gelman suggests to diagnose on more than one chain
    for (int i = 2; i < argc; i++)
    {
        MCMCFile.push_back(std::string(argv[i]));
        std::cout<<"Adding file: "<<MCMCFile.back()<<std::endl;
        Nchains++;
    }
    
    if(Ntoys < 1) 
    {
        std::cerr << "You specified " << Ntoys << " specify larger greater than 0" << std::endl;
        exit(-1);
    }
        
    if(Nchains == 1)
    {
        std::cout<<" Gelman is going to be sad :(. He suggested you should use more than one chain (at least 4). Code works fine for one chain, however, estimator might be biased."<<std::endl;
        std::cout<<" Multiple chains are more likely to reveal multimodality and poor adaptation or mixing:"<<std::endl;
    }
    std::cout<<"Diagnosing "<<Nchains<<" chians, with "<<Ntoys<<" toys"<<std::endl;
    
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
    
    TRandom3 *rnd = new TRandom3(0);

    std::cout << "Generating " << Ntoys << std::endl;
    TStopwatch clock;
    clock.Start();
    
    TChain **Chain = new TChain*[Nchains];
    
    int *BurnIn = new int[Nchains]();
    int *nEntries = new int[Nchains]();
    int *nBranches = new int[Nchains]();
    
    int step[Nchains];
    // Open the Chain
    //It is tempting to multithread here but unfortunately, ROOT files are not thread safe :(
    for (int m = 0; m < Nchains; m++) 
    {
        Chain[m] = new TChain("posteriors");
        Chain[m]->Add(MCMCFile[m].c_str());

        nEntries[m] = Chain[m]->GetEntries();

        // Set the step cut to be 20%
        BurnIn[m] = nEntries[m]/5;
    
        // Get the list of branches
        TObjArray* brlis = (TObjArray*)(Chain[m]->GetListOfBranches());

        // Get the number of branches
        nBranches[m] = brlis->GetEntries();

        if(m == 0) BranchNames.reserve(nBranches[m]);

        // Set all the branches to off
        Chain[m]->SetBranchStatus("*", false);
        
        // Loop over the number of branches
        // Find the name and how many of each systematic we have
        for (int i = 0; i < nBranches[m]; i++) 
        {
            // Get the TBranch and its name
            TBranch* br = (TBranch*)brlis->At(i);
            TString bname = br->GetName();

            // Read in the step
            if (bname == "step") {
                Chain[m]->SetBranchStatus(bname, true);
                Chain[m]->SetBranchAddress(bname, &step[m]);
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
                Chain[m]->SetBranchStatus(bname, true);
                if (VERBOSE) std::cout<<bname<<std::endl;
            }
        } 
    }
    nDraw = BranchNames.size();

    //KS: Sanity check that there is the same number of branches in all considered chains
    for (int m = 0; m < Nchains; m++) 
    {
      for (int mi = 0; mi < Nchains; mi++) 
      {
        if(nBranches[m] != nBranches[mi])
        {
            std::cerr<<"Ups, something went wrong, chain "<<m<<" called "<< MCMCFile[m] <<" has "<<nBranches[m]<<" branches, while "<<mi<<" called "<<MCMCFile[mi]<<" has "<<nBranches[mi]<<" branches"<<std::endl;
            std::cerr<<"All chains should have the same number of branches"<<std::endl;
            throw;
        }
      }
    }
    
    Draws = new double**[Nchains]();
    DrawsFolded = new double**[Nchains]();
    for(int m = 0; m < Nchains; m++)
    {
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
    }
  
  
    //KS: Now we will sample through MCMC posteriors distributions
    //It is tempting to multithread here but unfortunately, ROOT files are not thread safe :(
    for (int m = 0; m < Nchains; m++) 
    {
        if(BurnIn[m] >= nEntries[m])
        {
            std::cerr<<"You are running on a chain shorter than BurnIn cut"<<std::endl; 
            std::cerr<<"Number of entries "<<nEntries[m]<<" BurnIn cut "<<BurnIn<<std::endl;
            std::cerr<<"You will run into the infinite loop"<<std::endl;
            std::cerr<<"You can make a new chain or modify BurnIn cut"<<std::endl;
            std::cerr<<"Find me here: "<<__FILE__ << ":" << __LINE__ << std::endl;
            throw;
        }
            
        for (int i = 0; i < Ntoys; i++) 
        {
            // Get a random entry after burn in
            int entry = (int)(nEntries[m]*rnd->Rndm());
            
            Chain[m]->GetEntry(entry);
    
            // If we have combined chains by hadd need to check the step in the chain
            // Note, entry is not necessarily the same as the step due to merged ROOT files, so can't choose an entry in the range BurnIn - nEntries :(
            if (step[m] < BurnIn[m]) 
            {
                i--;
                continue;
            }
            
            // Output some info for the user
            if (Ntoys > 10 && i % (Ntoys/10) == 0) {
                std::cout << "On toy " << i << "/" << Ntoys << " (" << int(double(i)/double(Ntoys)*100.0) << "%)" << std::endl;
                std::cout << "   Getting random entry " << entry << std::endl;
            }
        
            // Set the branch addresses for params
            for (int j = 0; j < nDraw; ++j) 
            {
                Chain[m]->SetBranchAddress(BranchNames[j].Data(), &Draws[m][i][j]);
            }
            
            Chain[m]->GetEntry(entry);

        }//end loop over toys
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
        double* TempDraws = new double[Ntoys*Nchains]();
        for(int m = 0; m < Nchains; m++)
        {
            for(int i = 0; i < Ntoys; i++)
            {
                const int im = i+m;
                TempDraws[im] = Draws[m][i][j];
            }
        }
        MedianArr[j] = CalcMedian(TempDraws, Ntoys*Nchains);
        delete[] TempDraws;
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
    delete rnd;
    
    delete[] BurnIn;
    delete[] nEntries;
    delete[] nBranches;
    
    for (int m = 0; m < Nchains; m++) 
    {
        delete Chain[m];
    }
    delete[] Chain;
    
    clock.Stop();
    std::cout << "Finsihed producing Toys, it took " << clock.RealTime() << "s to finish" << std::endl;
}

// *******************
// Create all arrays we are going to use later
void InitialiseArrays() {
// *******************
  
    std::cout << "Initialising arrays " << std::endl;

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
    //KS: Additionaly calculates effective step size which is an estimate of the sample size required to achieve the same level of precision if that sample was a simple random sample.
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
    std::cout << "Finsihed calcaulting RHat, it took " << clock.RealTime() << "s to finish" << std::endl;
}


// *******************
void SaveResults() {
// *******************    
    std::string NameTemp = "";
    //KS: If we run over many many chains there is danger that name will be so absurdly long we run over systmat limit and job will be killed :(
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
    else
    {
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
    std::cout<<" Number of parameters which has R hat greater than 1.1 is "<<Criterium<<"("<<100*double(Criterium)/double(nDraw)<<"%) while for R hat folded "<<CiteriumFolded<<"("<<100*double(CiteriumFolded)/double(nDraw)<<"%)"<<std::endl;

    for(int j = 0; j < nDraw; j++)
    {
      if( (RHat[j] > 1.1 || RHatFolded[j] > 1.1) && ValidPar[j])
      {
        std::cout<<"Parameter "<<BranchNames[j]<<" has R hat higher than 1.1"<<std::endl;
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
    TCanvas *TempCanvas = new TCanvas("Canvas", "Canvas", 1024, 1024);
    gStyle->SetOptStat(0);
    TempCanvas->SetGridx();
    TempCanvas->SetGridy();
    
    // Random line to write useful information to TLegend
    TLine *TempLine = new TLine(0 , 0, 0, 0);
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
    
    Legend->AddEntry(TempLine, Form("Number of throws=%.0i, Number of chains=%.1i", Ntoys, Nchains), "");
    Legend->AddEntry(RhatPlot, "Rhat Gelman 2013", "l");
    Legend->AddEntry(RhatFoldedPlot, "Rhat-Folded Gelman 2021", "l");
      
    RhatPlot->Draw();
    RhatFoldedPlot->Draw("same");
    Legend->Draw("same");
    TempCanvas->Write("Rhat");
    delete Legend;
    Legend = NULL;
    
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

    Legend->AddEntry(TempLine, Form("Number of throws=%.0i, Number of chains=%.1i", Ntoys, Nchains), "");
    Legend->AddEntry(RhatLogPlot, "Rhat Gelman 2013", "l");
    Legend->AddEntry(RhatFoldedLogPlot, "Rhat-Folded Gelman 2021", "l");

    RhatLogPlot->Draw();
    RhatFoldedLogPlot->Draw("same");
    Legend->Draw("same");
    TempCanvas->Write("RhatLog");
    delete Legend;
    Legend = NULL;

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
    
    Legend->AddEntry(TempLine, Form("Number of throws=%.0i, Number of chains=%.1i", Ntoys, Nchains), "");
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
    delete TempCanvas;
    delete TempLine;
    
    delete RhatLogPlot;
    delete RhatFoldedLogPlot;

    DiagFile->Close();
    delete DiagFile;
    
    std::cout<<"Finished and wrote results to "<<NameTemp<<std::endl;
}

// *******************
//KS: Pseudo destructor
void DestroyArrays() {
// *******************
    
    std::cout<<"Killing all arrays"<<std::endl;
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
double CalcMedian(double arr[], int size) {
// *******************   
   std::sort(arr, arr+size);
   if (size % 2 != 0)
      return (double)arr[size/2];
   return (double)(arr[(size-1)/2] + arr[size/2])/2.0;
}


// *******************
//calculate median
void CapVariable(double var, double cap) {
// *******************   

    if(std::isnan(var) || !std::isfinite(var)) var = cap;
}
