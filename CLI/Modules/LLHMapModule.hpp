/// @file LLHMapModule.hpp
/// @brief Module for processing LLHMap outputs

#pragma once
#include "CLI/API/plugin.hpp"
#include "Manager/Core.h"

_MaCh3_Safe_Include_Start_ //{
// yaml Includes
#include "yaml-cpp/yaml.h"
// ROOT includes
#include <ROOT/RDataFrame.hxx>
_MaCh3_Safe_Include_End_ //}

namespace M3{

  /// @class LLHMapModule
  /// @brief Module for processing LLHMap outputs
  ///
  /// This module provides functionality for processing LLHMap outputs, including
  /// generating 1D and 2D profiled likelihoods as defined in the configuration.
  class LLHMapModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~LLHMapModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;

      /// @brief Process the parameters as given by the config file and present in the input LLHMap
      /// @param PoIs Vector of parameters of interest
      /// @param Map The input LLHMap as a ROOT RDataFrame
      /// @return Vector of parameters that will be processed
      std::vector<std::string> GetParams(std::vector<std::string> &PoIs, ROOT::RDataFrame Map);

      /// @brief Extract the binning for a given parameter from the LLHMap
      /// @param param The parameter for which to extract the binning
      /// @param Settings The YAML configuration node
      /// @param Map The input LLHMap as a ROOT RDataFrame
      /// @return A pair containing a boolean indicating success and a vector of bin edges
      std::pair<bool,std::vector<double>> ExtractBinning(std::string param, YAML::Node Settings, ROOT::RDataFrame Map);
  };
}
