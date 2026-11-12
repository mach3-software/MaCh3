#pragma once
#include "CLI/API/plugin.hpp"



// MACH3 PLOTTING
#include "Utils/Plotting/PlottingManager.h"

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file LLHModule.hpp
/// @ingroup MaCh3Plotting
/// @author Ewan Miller

namespace M3{

  /// @class LLHModule
  /// @brief Module for plotting log-likelihood (LLH) scans and comparisons
  ///
  /// This module provides functionality to create and manage plots for log-likelihood
  /// scans, including split-sample comparisons and ratio plots.
  class LLHModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~LLHModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;

    private:

      /// @brief Rebuild the arguments used to invoke the PlotLLH module
      /// @return The reconstructed command-line arguments
      std::vector<std::string> rebuildCommandLine();

      void validateCommandLine();

      void SetTPads(TPad*& LLHPad, TPad*& ratioPad);

      void getSplitSampleStack(int fileIdx, const std::string& parameterName, const TH1D& LLH_allSams,
                               std::vector<float> &cumSums, std::vector<bool> &drawLabel,
                               THStack *sampleStack, TLegend *splitSamplesLegend,
                               float baselineLLH_main = 0.00001);

      // handy dandy helper function for drawing and nicely formatting a stack of ratio plots
      void drawRatioStack(THStack *ratioCompStack);

      void makeLLHScanComparisons(const std::string& paramName,
                                  const std::string& LLHType,
                                  const std::string& outputFileName,
                                  const std::unique_ptr<TCanvas>& canv);

      void makeSplitSampleLLHScanComparisons(const std::string& paramName,
                                       const std::string& outputFileName,
                                       const std::unique_ptr<TCanvas>& canv);
      
      int PlotLLH();

      /// @warning KS: keep raw pointer or ensure manual delete of PlotMan. If spdlog in automatically deleted before PlotMan then destructor has some spdlog and this could cause segfault
      Plotting::PlottingManager *PlotMan;
  };
}
