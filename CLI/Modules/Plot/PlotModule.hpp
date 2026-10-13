/// @file PlotModule.hpp
/// @brief Module for plotting MCMC results

#pragma once
#include "CLI/API/plugin.hpp"
#include "CLI/Modules/Plot/LLHMapModule.hpp"
#include "CLI/Modules/Plot/LLHModule.hpp"
#include "CLI/Modules/Plot/PostFitParamModule.hpp"

namespace M3{

  /// @class PlotModule
  /// @brief Module for plotting MCMC results
  ///
  /// This module provides diagnostics such as convergence tests, autocorrelation
  /// analysis, and batch means calculations to assess the quality of MCMC samples.
  class PlotModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~PlotModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;
    private:
      LLHMapModule m_llhMapModule;
      LLHModule m_llhModule;
      PostFitParamModule m_postFitParamModule;
  };
}
