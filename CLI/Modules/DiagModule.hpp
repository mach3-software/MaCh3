#pragma once
#include "CLI/API/plugin.hpp"
#include "CLI/Modules/DiagMCMCModule.hpp"
#include "CLI/Modules/LLHMapModule.hpp"


namespace M3{
/// @class DiagModule
/// @brief Module for MCMC diagnostics
///
/// This module provides diagnostics such as convergence tests, autocorrelation
/// analysis, and batch means calculations to assess the quality of MCMC samples.
class DiagModule: public IModuleBase{
 public:
  /// @brief Destructor
  virtual ~DiagModule();

  /// @brief Get the argument parser for this module
  /// @return Pointer to the configured MaCh3ArgumentParser
  MaCh3ArgumentParser* get_parser() override;

  /// @brief Execute the MCMC diagnostics
  /// @return Exit code (0 on success)
  int Run() override;
 private:
  LLHMapModule m_llhMapModule;
  DiagMCMCModule m_diagMCMCModule;
};
}
