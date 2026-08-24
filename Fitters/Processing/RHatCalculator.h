#pragma once

// MaCh3 includes
#include "Fitters/Processing/StatisticalUtils.h"
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

    std::vector<std::vector<double>> Mean;
    std::vector<std::vector<double>> StandardDeviation;

    std::vector<double> MeanGlobal;
    std::vector<double> StandardDeviationGlobal;

    std::vector<double> BetweenChainVariance;
    std::vector<double> MarginalPosteriorVariance;
    std::vector<double> RHat;
    std::vector<double> EffectiveSampleSize;

    ///// High mem only /////
    int Ntoys;
    std::vector<std::vector<std::vector<double>>> Draws;
    std::vector<std::vector<std::vector<double>>> DrawsFolded;

    std::vector<double> MedianArr;

    std::vector<std::vector<double>> MeanFolded;
    std::vector<std::vector<double>> StandardDeviationFolded;

    std::vector<double> MeanGlobalFolded;
    std::vector<double> StandardDeviationGlobalFolded;

    std::vector<double> BetweenChainVarianceFolded;
    std::vector<double> MarginalPosteriorVarianceFolded;
    std::vector<double> RHatFolded;
    std::vector<double> EffectiveSampleSizeFolded;

    ///// Low mem Only /////
    std::vector<int> Ntoys_requested;
    std::vector<int> Ntoys_filled;
    int TotToys;
    unsigned int NThin;

    /// Sum_i^N x_i | total
    std::vector<double> S1_global;
    /// Sum_i^N x_i^2 | total
    std::vector<double> S2_global;
    /// Sum_i^N x_i | for each chain
    std::vector<std::vector<double>> S1_chain;
    /// Sum_i^N x_i^2 | for each chain
    std::vector<std::vector<double>> S2_chain;
};
