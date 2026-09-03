/// @file ProcessMCMCModule.hpp
/// @brief Module for processing MCMC chains and producing diagnostic plots

//MaCh3 includes
#pragma once
// yaml Includes
#include "yaml-cpp/yaml.h"
#include "Fitters/Processing/MCMCProcessor.h"
#include "CLI/API/plugin.hpp"

namespace M3{

/// @file ProcessMCMCModule.hpp
/// @brief Main executable responsible for different types of MCMC processing like drawing posteriors, triangle plots etc. Actual implementation of methods is in MCMCProcessor
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author Kamil Skwarczynski

/// @class ProcessMCMCModule
/// @brief Module for processing MCMC chains and producing various diagnostic plots
///
/// This module provides functionality for:
/// - Drawing posterior distributions
/// - Creating triangle plots
/// - Computing Bayes factors
/// - Performing Savage-Dickey tests
/// - Calculating parameter evolution
/// - Comparing multiple MCMC chains
  class ProcessMCMCModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~ProcessMCMCModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC processing
      /// @return Exit code (0 on success)
      int Run() override;


    private:
      /// @brief Parse custom binning edges from YAML configuration
      /// @param Settings YAML node containing CustomBinEdges section
      /// @return Map of parameter names to their (min, max) bin edges
      std::map<std::string, std::pair<double, double>> GetCustomBinning(const YAML::Node& Settings);

      /// @brief Process a single MCMC chain
      /// @param inputFile Path to the MCMC chain ROOT file
      void ProcessMCMC(const std::string& inputFile);

      /// @brief Compare and process multiple MCMC chains
      void MultipleProcessMCMC();

      /// @brief Calculate Bayes factors for hypothesis testing
      /// @param Processor Pointer to the MCMCProcessor instance
      void CalcBayesFactor(MCMCProcessor* Processor);

      /// @brief Calculate Savage-Dickey ratios for Bayes factor estimation
      /// @param Processor Pointer to the MCMCProcessor instance
      void CalcSavageDickey(MCMCProcessor* Processor);

      /// @brief Calculate parameter evolution over MCMC steps
      /// @param Processor Pointer to the MCMCProcessor instance
      void CalcParameterEvolution(MCMCProcessor* Processor);

      /// @brief Create bipolar plots for parameter visualization
      /// @param Processor Pointer to the MCMCProcessor instance
      void CalcBipolarPlot(MCMCProcessor* Processor);

      /// @brief Generate triangle plots showing parameter correlations
      /// @param Processor Pointer to the MCMCProcessor instance
      void GetTrianglePlot(MCMCProcessor* Processor);

      /// @brief Diagnose covariance matrix stability across burn-in cuts
      /// @param Processor Pointer to the MCMCProcessor instance
      /// @param inputFile Path to the input file for naming output
      void DiagnoseCovarianceMatrix(MCMCProcessor* Processor, const std::string& inputFile);

      /// @brief Convert TMatrixDSym to TH2D histogram for plotting
      /// @param Matrix Pointer to the symmetric matrix
      /// @param title Title for the resulting histogram
      /// @return Pointer to the created TH2D histogram
      TH2D* TMatrixIntoTH2D(TMatrixDSym* Matrix, const std::string& title);

      /// @brief Perform Kolmogorov-Smirnov test between posterior distributions
      /// @param Processor Vector of MCMCProcessor instances
      /// @param Posterior Canvas for drawing results
      /// @param canvasname Name for the output PDF
      void KolmogorovSmirnovTest(const std::vector<std::unique_ptr<MCMCProcessor>>& Processor,
                                 const std::unique_ptr<TCanvas>& Posterior,
                                 const TString& canvasname);

      int nFiles;                                ///< Number of MCMC files being processed
      std::vector <std::string> FileNames;       ///< List of MCMC chain file paths
      std::vector <std::string> TitleNames;      ///< List of titles for each chain
      std::string config;                        ///< Path to the configuration file
  };
}
