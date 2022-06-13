// Compile with something like
//
// g++ `root-config --cflags` -g -o DiagMCMC DiagMCMC.cpp -I`root-config --incdir` `root-config --glibs --libs`
//
// Run on ND280 only fits to get xsec and beam posteriors prettified. Essentially a compiled makeTransferMatrixAll
//
// ./DiagMCMC will tell you how
//
// Mostly self-contains class and executable
//
// MCMC stuff to implement:
// Trace plots            -- DONE
// LogL vs step plots     -- DONE
// Acceptance probability -- DONE
// Autocorrelation        -- DONE
// _Batched Means_        -- DONE

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <math.h>

#include "TChain.h"
#include "TFile.h"

#include "TObjArray.h"
#include "TBranch.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"
#include "TStopwatch.h"

#include "TF1.h"
#include "TH1.h"

#include "TVectorD.h"
#include "TCanvas.h"
#include "TColor.h"

#ifdef MULTITHREAD
#include "omp.h"
#endif

class MCMCDiag_Handler {

  public:
    MCMCDiag_Handler(std::string FileName);
    ~MCMCDiag_Handler() { std::cout << "I wrote to " << OutputFileName << std::endl; }

    // Diagnose the MC
    void DiagMCMC();

    // Set the cross-section covariance matrix from external
    void SetXsecCov(std::string FileName) {XsecCov = FileName;};

    // Set the number of entries externally
    void SetNEntries(const int UserEntries) {
      std::cout << "User specified " << UserEntries << " events, over-riding" << std::endl;
      nEntries = UserEntries;
    };

  private:
    // Checks the XsecCov string is set
    void CheckXsecCov();

    // Compare the beginning of one string with another
    bool BeginString(std::string NameString, std::string CompareString);

    // The input file from which we read everything
    TFile *InputFile;

    // The TChain for the diagnoser
    TChain *Chain;
    std::string XsecCov;

    // The output file which we write everything to
    TFile *OutputFile;
    std::string OutputFileName;

    // Vector of title for the cross-section parameters
    std::vector<std::string> XsecNames;

    // Some helpers
    void ReadXsecCov();
    void CountBranches();
    void ReadInputTree();
    void ParamTraces();
    void AutoCorrelation();
    void BatchedMeans();
    void AcceptanceProbabilities();
    
    // Number of parameters
    int nParams;
    int nXsec;
    int nFlux;
    int nSamples;
    int nSysts;
    int nBranches;
    int nEntries;

    int nBatches;

    // Limits for cross-section
    std::vector<double> XsecLimLow;
    std::vector<double> XsecLimHigh;
    std::vector<double> XsecCentral;

    // Holds all the parameter variations
    double **ParamValues;
    double *ParamSums;
    double **BatchedAverages;

    // Holds the sample values
    double **SampleValues;
    // Holds the systs values
    double **SystValues;

    // Holds all accProb
    double *AccProbValues;
    double *AccProbBatchedAverages;
    
    // Trace plots
    TH1D **TraceParamPlots;
    TH1D **TraceSamplePlots;
    TH1D **TraceSystsPlots;
    TH1D **BatchedParamPlots;

    // LagK autocorrelation plots
    TH1D **LagKPlots;

    // Acceptance Prob Plots
    TH1D *AcceptanceProbPlot;
    TH1D *BatchedAcceptanceProblot;
    
    // Vector of the total logL for each step
    std::vector<double> LogL_v;
    std::vector<std::string> ParamName_v;
    std::vector<std::string> SampleName_v;
    std::vector<std::string> SystName_v;

};
