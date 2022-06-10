#ifndef __MCMCPROCESSOR_H_
#define __MCMCPROCESSOR_H_

#ifndef __UNDEF__
#define __UNDEF__ 1234567890
#endif

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TObjArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TColor.h"
#include "TStyle.h"

// Class to process MCMC output produced by mcmc::runMCMC
// Useful for when we want to extract values from a previous MCMC 
//  and want to send off another MCMC with those values, or perhaps to LLH scans at central values
// Mostly taken from nd280_utils/DrawComp.cpp
//
// Return: Postfit parameters, Postfit covariance matrix
//
// Make postfit pmu cosmu dists
// Make LLH scans around output

enum ParameterEnum {
  kFluxPar = 0,
  kXSecPar = 1,
  kNearDetPar = 2
};

class MCMCProcessor {
  public:
    MCMCProcessor(const std::string &InputFile, bool MakePostfitCorr);
    ~MCMCProcessor();

    // Get the post-fit results (arithmetic and Gaussian)
    void GetPostfit(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Central_Gauss, TVectorD *&Errors_Gauss, TVectorD *&Peaks);
    // Get the post-fit covariances and correlations
    void GetCovariance(TMatrixDSym *&Cov, TMatrixDSym *&Corr);
    // Or the individual post-fits
    void GetPostfit_Ind(TVectorD *&Central, TVectorD *&Errors, TVectorD *&Peaks, ParameterEnum kParam);

    void MakePostfit();

    // Get the number of parameters
    int GetNParams() { return nDraw; };
    int GetNFlux() { return nFlux; };
    int GetNXSec() { return nXSec; };
    int GetNND() { return nNear; };

    std::string const & GetFluxCov() const { return FluxCov; };
    std::string const & GetXSecCov() const { return XSecCov; };
    std::string const & GetNearCov() const { return NearCov; };
    std::string const & GetNDruns() const { return NDruns; };
    std::vector<std::string> const & GetNDsel() const {return NDsel;};

    // Draw the post-fit comparisons
    void DrawPostfit();
    // Draw the post-fit covariances
    void DrawCovariance();

    // Get the vector of branch names
    const std::vector<TString>& GetBranchNames() const { return BranchNames;};

    // Set the step cutting
    // Either by string
    void SetStepCut(std::string Cuts);
    // Or by int
    void SetStepCut(int Cuts);

  private:
    inline TH1D* MakePrefit();
    inline void MakeCovariance();
    inline void MakeOutputFile();

    inline void ReadInputCov();
    inline void FindInputFiles();
    inline void ReadFluxFile();
    inline void ReadXSecFile();
    inline void ReadNearFile();

    inline void ScanInput();
    inline void SetupOutput();

    std::string MCMCFile;

    // Cross-section covariance matrix name position
    std::string XSecCov;
    // Flux covariance matrix name position
    std::string FluxCov;
    // Near cov
    std::string NearCov;
    // ND runs
    std::string NDruns;
    // ND selections
    std::vector<std::string> NDsel;

    TChain *Chain;
    std::string StepCut;
    int nBranches;

    std::vector<TString> BranchNames;
    std::vector<TString> FluxNames;
    std::vector<TString> XSecNames;
    std::vector<TString> NearNames;

    // Is the ith parameter varied
    std::vector<bool> IamVaried;

    std::vector<double> FluxCentral;
    std::vector<double> FluxErrors;
    std::vector<double> XSecCentral;
    std::vector<double> XSecErrors;
    std::vector<double> NearCentral;
    std::vector<double> NearErrors;

    std::string OutputName;
    TString CanvasName;

    bool PlotFlux;
    bool PlotXSec;
    bool PlotDet;

    bool MakeCorr;
    bool MadePostfit;

    // The output file
    TFile *OutputFile;
    // Directory for posteriors
    TDirectory *PostDir;

    int nDraw;
    int nFlux;
    int nXSec;
    int nNear;
    int nEntries;

    // Gaussian fitter
    TF1 *Gauss;

    // Make an enum for which class this parameter belongs to so we don't have to keep string comparing
    std::vector<ParameterEnum> ParamType;

    TCanvas *Posterior;

    TVectorD *Means;
    TVectorD *Errors;
    TVectorD *Means_Gauss;
    TVectorD *Errors_Gauss;
    TVectorD *Peaks;

    TMatrixDSym *Covariance;
    TMatrixDSym *Correlation;

    // Number of bins
    int nBins;
    // Drawrange for SetMaximum
    double DrawRange;
};

#endif
