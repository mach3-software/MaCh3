/// @file RHatModule.hpp
/// @brief Module for computing the R-hat diagnostic for MCMC chains

#pragma once
#include "CLI/API/plugin.hpp"

namespace M3{

  /// @class RHatModule
  /// @brief Module for computing the R-hat diagnostic for MCMC chains
  ///
  /// This module calculates the R-hat statistic to assess the convergence of MCMC chains.
  class RHatModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~RHatModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;
  };
}
