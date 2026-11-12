/// @file UmbrellaSolver.cpp
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author David Riley
#pragma once
#include "CLI/API/plugin.hpp"
#include "CLI/Modules/UmbrellaPlotModule.hpp"
#include <vector>
#include <string>
#include "Fitters/Algorithms/MulticanonicalMCMCHandler.h"


_MaCh3_Safe_Include_Start_ //{
#include "TSystem.h"
#include "TChain.h"
#include "TSystemDirectory.h"
_MaCh3_Safe_Include_End_ //}

namespace M3{
bool debug_mode = false;
/// Structure to hold each window configuration
/// These are individual window settings, function type & params plus filename
struct WindowConfig {
  std::string name;
  M3::BiasFunction umbrellaBiasFunction;
  double center;
  double width;
  std::string input_file;
  double multicanonical_beta;
  double vonMises_kappa;
};

/// Structure to hold umbrella configuration
/// This is the global settings of the umbrella fit eg number of windows and solve settings
struct UmbrellaConfig {
  std::vector<WindowConfig> windows;
  std::string output_file;
  std::string variable_of_interest;
  bool dynamic_files;
  std::string dynamic_pattern;
  int dynamic_n_windows;
  int max_iterations;
  double tolerance;
  int print_frequency;
};

/// @class UmbrellaModule
/// @brief Module for performing Umbrella sampling
///
/// This module provides functionality for performing Umbrella sampling, a technique
/// used to enhance sampling of rare events in MCMC simulations.
class UmbrellaModule: public IModuleBase{
 public:
  /// @brief Destructor
  virtual ~UmbrellaModule();

  /// @brief Get the argument parser for this module
  /// @return Pointer to the configured MaCh3ArgumentParser
  MaCh3ArgumentParser* get_parser() override;

  /// @brief Execute the Umbrella sampling
  /// @return Exit code (0 on success)
  int Run() override;
 private:
  /// YAML-based config parser using yaml-cpp library
  UmbrellaConfig parseYAMLConfig(const std::string &filename);
  // Gaussian window function
  double gaussianWindow(double x, double center, double width);
  // von Mises window function (circular analog of Gaussian)
  double vonMisesWindow(double x, double center, double kappa);
  // Generalised Gaussian window allows harsher constraints on the tails of the distribution,
  // better for controlling windows in highly disfavoured regions
  double generalisedGaussian2(double x, double mean, double width);
  /// @todo modify this to instead use atan2 implementation for wrapping
  double GetMulticanonicalWeightGenGaussian(double deltacp, double mean, double width, double multicanonicalBeta);
  // A sub calculation for the overlap matrix
  // Sum of all windows weighted by z values
  double summedWindowsWeighted(double x,
                                const std::vector<WindowConfig> &windows,
                                const std::vector<double> &z_values);
  size_t calculateWindowCacheSizeBytes(const std::vector<WindowConfig> &windows,
                                        const std::vector<std::vector<double>> &samples);
  // Precompute all window evaluations once: cache[i][j][s] = window_j evaluated
  // at samples[i][s]
  // Memory heavy depending on number of steps/cores/windows
  std::vector<std::vector<std::vector<double>>> buildWindowCache(const std::vector<WindowConfig> &windows,
                                                                  const std::vector<std::vector<double>> &samples);
  // Main function to calculate the F matrix from a given set of samples plus weights
  std::vector<std::vector<double>> calcFmatrix(const std::vector<double> &z_current,
                                                const std::vector<WindowConfig> &windows,
                                                const std::vector<std::vector<double>> &samples,
                                                const std::vector<std::vector<std::vector<double>>> &window_cache);
  // Z-solver function implementing the fixed point matrix iteration algorithm
  std::vector<double> zSolver(const std::vector<double> &z_current,
                              const std::vector<WindowConfig> &windows,
                              const std::vector<std::vector<double>> &samples,
                              const std::vector<std::vector<std::vector<double>>> &window_cache,
                              [[maybe_unused]] bool verbose = false,
                              [[maybe_unused]] int *total_lines = nullptr);
  /// A few different convergence checks
  // Check convergence
  bool checkConvergence(const std::vector<double> &z_current, const std::vector<double> &z_prev, double tolerance);
  std::vector<double> getZDiffs(const std::vector<double> &z_current, const std::vector<double> &z_prev);
  // This implements a check on stalling of the evolution of the window weights. Seems to be a
  // better definition for convergence that above
  // moving-average stalled-convergence check on z values
  bool checkConvergenceStalled(const std::vector<double> &z_current,
                                const std::vector<double> &z_prev,
                                double tolerance);
  void LoadFiles(std::vector<std::vector<double>>& samples,
                  UmbrellaConfig& config,
                  std::vector<TFile *>& input_files,
                  std::vector<TTree *>& input_trees);
  void DiagnoseUmbrella(const UmbrellaConfig& config,
                        TFile * output_file,
                        const std::vector<std::vector<double>>& z_evolution,
                        TTree *combined_tree);
  // Main function to run the umbrella sampling solver
  void UmbrellaSolver(const std::string &config_file);

  std::unique_ptr<MaCh3ArgumentParser> m_default_parser;
  UmbrellaPlotModule m_plotModule;
};
}
