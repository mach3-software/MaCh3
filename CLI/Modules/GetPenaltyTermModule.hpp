/// @file GetPenaltyTermModule.hpp
/// @brief Module for extracting penalty terms from systematic chains
///
/// This module is designed to retrieve penalty terms from various sources,
/// such as flux and cross-section systematic chains. Since flux and cross-section
/// uncertainties are handled systematically, the penalty term cannot be taken
/// directly from the chain.
///
/// @todo This should really be moved to MCMC Processor
/// @ingroup MaCh3DiagnosticProcessing
/// @author Kamil Skwarczynski

#pragma once
#include "yaml-cpp/yaml.h"
#include "CLI/API/plugin.hpp"

namespace M3{

  /// @class GetPenaltyTermModule
  /// @brief Module for extracting specific penalty terms from MCMC chains
  ///
  /// This module retrieves penalty terms from flux and cross-section systematic
  /// chains, since these systematic uncertainties are handled separately and the
  /// penalty term cannot be taken directly from the chain.
  class GetPenaltyTermModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~GetPenaltyTermModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the penalty term extraction
      /// @return Exit code (0 on success)
      int Run() override;

    private:
      /// @brief Read covariance matrix and parameter information from file
      /// @param inputFile Path to the ROOT file containing the covariance
      /// @param Prior Vector to store prior parameter values
      /// @param isFlat Vector to store flat prior flags
      /// @param ParamNames Vector to store parameter names
      /// @param invCovMatrix 2D vector to store inverted covariance matrix
      /// @param nParams Reference to store the number of parameters
      void ReadCovFile(const std::string& inputFile,
                       std::vector <double>& Prior,
                       std::vector <bool>& isFlat,
                       std::vector<std::string>& ParamNames,
                       std::vector<std::vector<double>>& invCovMatrix,
                       int& nParams);

      /// @brief Load penalty term sets from YAML configuration
      /// @param Settings YAML configuration node
      /// @param SetsNames Vector to store set names
      /// @param FancyTitle Vector to store fancy titles for each set
      /// @param isRelevantParam 2D vector indicating which parameters belong to each set
      /// @param ParamNames Vector of all parameter names
      /// @param nParams Total number of parameters
      void LoadSettings(YAML::Node& Settings,
                        std::vector<std::string>& SetsNames,
                        std::vector<std::string>& FancyTitle,
                        std::vector<std::vector<bool>>& isRelevantParam,
                        const std::vector<std::string>& ParamNames,
                        const int nParams);

      /// @brief Calculate and plot penalty terms for parameter sets
      /// @param inputFile Path to the MCMC chain ROOT file
      /// @param configFile Path to the YAML configuration file
      void GetPenaltyTerm(const std::string& inputFile, const std::string& configFile);

  };
}
