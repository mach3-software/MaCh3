#pragma once

// MaCh3 includes
#include "Fitters/StatisticalUtils.h"
#include "Samples/HistogramUtils.h"
#include "Parameters/ParameterHandlerUtils.h"

/// @brief Class responsible for calculating \f$ \hat{R} \f$ estimator for Markov Chain Monte Carlo (MCMC) convergence.
///
/// @cite gelman2019.
///
/// @author Kamil Skwarczynski
/// @author Michael Reh
/// @author Tomas Nosek
class RHatCalculator {
  public:
    /// @brief Constructor
    /// @param HighMemory High mem is generally slower but has calculates additional info
    /// @param Inputs vector of paths to MCMC PrepareChains
    /// @param entries if in high mem mode this number of toys per chain, for low mem thinning setting
    RHatCalculator(bool HighMemory, std::vector<std::string>& Inputs, int entries);
    /// @brief Destroys the RHatCalculator object.
    virtual ~RHatCalculator();

    void RunDiagnostic();
  private:
    /// @brief Load chain and prepare toys
    void PrepareChains_HighMem();
    /// @brief Load chain and prepare toys
    void PrepareChains();

    /// @brief Create all arrays we are going to use later
    void InitialiseArrays();

    /// @brief KS: Based on Gelman et. al. arXiv:1903.08008v5
    void CalcRhat_HighMem();
    /// @brief KS: Based on Gelman et. al. arXiv:1903.08008v5
    void CalcRhat();

    void SaveResults();

    ///// Common /////
    bool HighMemoryMode;
    int nDraw;
    int Nchains;

    std::vector<TString> BranchNames;
    std::vector<std::string> MCMCFile;
    std::vector<bool> ValidPar;

    double** Mean;
    double** StandardDeviation;

    double* MeanGlobal;
    double* StandardDeviationGlobal;

    double* BetweenChainVariance;
    double* MarginalPosteriorVariance;
    double* RHat;
    double* EffectiveSampleSize;

    ///// High mem only /////
    int Ntoys;

    double ***Draws;

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

    ///// Low mem Only /////
    int* Ntoys_requested;
    int* Ntoys_filled;
    int TotToys;
    unsigned int NThin;

    /// Sum_i^N x_i   | total
    double* S1_global;
    /// Sum_i^N x_i^2 | total
    double* S2_global;
    /// Sum_i^N x_i   | for each chain
    double** S1_chain;
    /// Sum_i^N x_i^2 | for each chain
    double** S2_chain;
};
