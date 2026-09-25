/// @file ReweightModule.hpp
/// @brief Module for MCMC reweighting analysis

#pragma once
#include "Fitters/Processing/MCMCProcessor.h"
#include "CLI/API/plugin.hpp"
_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TGraph.h"
_MaCh3_Safe_Include_End_ //}

// C++ includes
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <map>

namespace M3{

    /// @brief Types of chain reweighting available
    enum kReweightType {
        kGaussian,     //!< Assumes gaussian prior
        kTGraph,       //!< Calculates Likelihood based on TGraph
        kTGraph2D,     //!< Calculates Likelihood based on TGraph2D
        kReweightTypes //!< This only enumerates
    };

    /// Structure to hold reweight configuration
    struct ReweightConfig {
        std::string key;       ///< The YAML key for this reweight
        std::string name;
        M3::kReweightType type;  ///< "Gaussian", "TGraph2D"
        int dimension;     ///< 1 or 2
        std::vector<std::string> paramNames; ///< Parameter names.
        std::vector<std::vector<double>> newPriorValues; ///< new [mean, sigma] pairs
        std::vector<std::vector<double>> oldPriorValues; ///< new [mean, sigma] pairs
        std::vector<bool> flatPrior;

        std::string weightBranchName; ///< Output weight branch name.
        bool enabled;

        // For TGraph 1D or 2D
        std::string fileName;         ///< ROOT file containing graph data.
        std::string graphName;        ///< Graph name in the ROOT file.

        // For TGraph1D
        std::unique_ptr<TGraph> graph_1D; ///< 1D interpolation graph.

        // For TGraph2D
        std::string hierarchyType; ///< "NO", "IO", or "auto"
        std::unique_ptr<TGraph2D> graph_NO; ///< Normal Ordering graph.
        std::unique_ptr<TGraph2D> graph_IO; ///< Inverted Ordering graph.
    };

  /// @class ReweightModule
  /// @brief Module for performing reweighting analysis on MCMC chains
  ///
  /// This module provides diagnostics such as convergence tests, autocorrelation
  /// analysis, and batch means calculations to assess the quality of MCMC samples.
  class ReweightModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~ReweightModule();

      /// @brief Function to interpolate 2D graph for Normal Ordering
      double Graph_interpolateNO(TGraph2D* graph, double theta13, double dm32);

      /// @brief Function to interpolate 2D graph for Inverted Ordering  
      double Graph_interpolateIO(TGraph2D* graph, double theta13, double dm32);

      /// @brief Function to interpolate 1D graph
      double Graph_interpolate1D(TGraph* graph, double theta13);

      /// @brief Get parameter information from MCMCProcessor
      bool GetParameterInfo(MCMCProcessor* processor, const std::string& paramName,
                            double& mean, double& sigma, bool& isFlat);

      /// @brief Load reweighting setting like 1D or 2D from YAML config
      void LoadReweightingSettings(std::vector<ReweightConfig>& reweightConfigs, const YAML::Node& reweight_settings);

      /// @brief Calculate 1D weight
      [[nodiscard]] double Get1DWeight(const ReweightConfig& rwConfig,
                                      const std::map<std::string, double>& paramValues);

      /// @brief Calculate 2D weight
      [[nodiscard]] double Get2DWeight(const ReweightConfig& rwConfig,
                                      const std::map<std::string, double>& paramValues);

      /// @brief Main executable responsible for reweighting MCMC chains
      /// @param inputFile MCMC Chain file path
      /// @param configFile Config file with reweighting settings
      /// @author David Riley
      /// @author Evan Goodman
      /// @todo add a generic 2D reweight that is not dm32 and theta13 specific DWR
      void ReweightMCMC(const std::string& configFile, const std::string& inputFile);

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC reweighting analysis
      /// @return Exit code (0 on success)
      int Run() override;
  };
}
