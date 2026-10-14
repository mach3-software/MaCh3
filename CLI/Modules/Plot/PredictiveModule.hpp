#pragma once
#include "CLI/API/plugin.hpp"

//MaCh3 Includes
#include "Utils/Plotting/PlottingUtils.h"
#include "Utils/Plotting/PlottingManager.h"
#include <numeric>

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file PredictiveModule.hpp
/// @ingroup MaCh3Plotting
/// @author Kamil Skwarczynski
/// @brief Module for predictive plotting analysis


namespace M3{

  /// @class PredictiveModule
  /// @brief Module for performing predictive plotting analysis
  ///
  /// This module provides functionality for generating predictive plots,
  /// including overlayed 1D histograms, violin plots, and mode-specific predictive plots.
  class PredictiveModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~PredictiveModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;
    private:
    std::vector<std::string> FindSamples(const std::string& File);
    std::vector<int> FindDimensions(const std::string& File, const std::vector<std::string>& Samples);
    std::vector<std::vector<std::string>> FindModes(const std::string& File,
                                                    const std::vector<std::string>& SampleNames);
    void PretifyHistogram(TH1* Hist, const std::string& SampleName);
    double GetPValue(const TH2D* hist);
    void PrintPosteriorPValue(const YAML::Node& Settings,
                              const std::vector<TFile*>& InputFiles,
                              const std::vector<std::string>& SampleNames);
    void OverlayViolin(const YAML::Node& Settings,
                       const std::vector<TFile*>& InputFiles,
                       const std::vector<std::string>& SampleNames,
                       const std::vector<int>& SampleDimension,
                       const std::unique_ptr<TCanvas>& canv);
    void OverlayPredicitve(const YAML::Node& Settings,
                           const std::vector<TFile*>& InputFiles,
                           const std::vector<std::string>& SampleNames,
                           const std::vector<int>& SampleDimension,
                           const std::unique_ptr<TCanvas>& canv);
    void OverlayPredicitveByMode(const YAML::Node& Settings,
                                 const std::vector<TFile*>& InputFiles,
                                 const std::vector<std::string>& SampleNames,
                                 const std::vector<int>& SampleDimension,
                                 const std::vector<std::vector<std::string>>& Modes,
                                 const std::unique_ptr<TCanvas>& canv);
    /// @brief KS: Get mean and error from gaussian fit to event distribution
    void GetMeanError(TH1D* hist, double &Mean, double &Error);
    /// @brief KS Print event rates in Latex like table
    void PrintPosteriorEventRates(const std::vector<TFile*>& InputFiles,
                                  const std::vector<std::string>& SampleNames);
    /// @brief KS: Print Fractional Uncertainties into Latex table format
    void PrintPosteriorFractionalUncertainties(const std::vector<TFile*>& InputFiles,
                                               const std::vector<std::string>& SampleNames);
    double GetLLH(TH1* hist);
    /// @brief KS Print Predictive LLH into Latex table format
    void PrintPredictiveLLH(const std::vector<TFile*>& InputFiles,
                            const std::vector<std::string>& SampleNames);
    void PredictivePlotting(const std::string& ConfigName,
                            const std::vector<std::string>& FileNames);
        /// @warning KS: keep raw pointer or ensure manual delete of PlotMan. If spdlog in automatically deleted before PlotMan then destructor has some spdlog and this could cause segfault
        Plotting::PlottingManager* PlotMan = nullptr;
  };
}
