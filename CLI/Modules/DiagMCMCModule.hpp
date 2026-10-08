/// @file DiagMCMCModule.hpp
/// @brief Module for MCMC diagnostic analysis

#pragma once
#include "CLI/API/plugin.hpp"
#include "CLI/Modules/DiagMCMCPlotModule.hpp"

namespace M3{

  /// @class DiagMCMCModule
  /// @brief Module for performing diagnostic analysis on MCMC chains
  ///
  /// This module provides diagnostics such as convergence tests, autocorrelation
  /// analysis, and batch means calculations to assess the quality of MCMC samples.
  class DiagMCMCModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~DiagMCMCModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;
    private:
      std::unique_ptr<MaCh3ArgumentParser> m_default_parser;
      DiagMCMCPlotModule m_plotModule;
  };
}
