/// @file UmbrellaPlotModule.hpp
/// @brief Module for Umbrella sampling plot analysis
/// @author David Riley
#pragma once
#include "CLI/API/plugin.hpp"
#include <string>

namespace M3{
/// @class UmbrellaPlotModule
/// @brief Module for Umbrella sampling plot analysis
///
/// This module provides functionality for generating plots related to Umbrella sampling,
/// such as free energy profiles and convergence diagnostics.
class UmbrellaPlotModule: public IModuleBase{
 public:
  /// @brief Destructor
  virtual ~UmbrellaPlotModule();

  /// @brief Get the argument parser for this module
  /// @return Pointer to the configured MaCh3ArgumentParser
  MaCh3ArgumentParser* get_parser() override;

  /// @brief Execute the Umbrella sampling plot analysis
  /// @return Exit code (0 on success)
  int Run() override;

 private:
  void UmbrellaPlotting(std::string InputFile, std::string outputDir);
};
}
