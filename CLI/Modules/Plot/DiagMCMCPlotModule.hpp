/// @file DiagMCMCPlotModule.hpp
/// @brief KS: This script is used to analyse output form DiagMCMC.
/// @warning This script support comparing up to 4 files, there is easy way to expand it up to five or six,
/// @todo this need serious refactor
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author Henry Wallace
/// @author Kamil Skwarczynski
#pragma once
#include <vector>
#include <memory>
#include <utility>
#include "CLI/API/plugin.hpp"

// MaCh3 includes
#include "Manager/Core.h"
_MaCh3_Safe_Include_Start_ //{
#include "TH1D.h"
#include "TString.h"
#include "TGraph.h"
#include "TDirectoryFile.h"
#include "TGraphAsymmErrors.h"
#include <TSystem.h>
_MaCh3_Safe_Include_End_ //}


namespace M3{

  /// @class DiagMCMCPlotModule
  /// @brief KS: This script is used to analyse output form DiagMCMC.
  class DiagMCMCPlotModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~DiagMCMCPlotModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics plotting
      /// @return Exit code (0 on success)
      int Run() override;

    private:
      /// @brief KS: function which looks for minimum in given range
      double GetMinimumInRange(TH1D *hist, const double minRange, const double maxRange);
      /// @brief HW: Check if histogram is flat within a given tolerance.
      bool IsHistogramAllOnes(TH1D *hist, double tolerance = 0.001, int max_failures = 100);
      /// @brief function which loops over Diag MCMC output and prints everything to PDF
      void MakeDiagPlot(std::vector<TString> input_files, std::vector<TString> hist_labels);
      /// @brief Plot autocorrelation for all params into single plot with color codding helping whether on average it is ok or not
      void PlotAutoCorr(const std::vector<TString>& fname);
      /// @brief HW: Create a band of minimum and maximum values from a histogram.
      std::pair<std::unique_ptr<TGraph>, std::unique_ptr<TGraph>> CreateMinMaxBand(TH1D *hist, Color_t color);
      /// @brief AR: Calculate a band of minimum and maximum values from a collection of histograms.
      /// @param histograms Vector of unique pointers to TH1D histograms
      /// @param color Color of the band
      /// @return Unique pointer to TGraphAsymmErrors representing the min-max band
      std::unique_ptr<TGraphAsymmErrors> CalculateMinMaxBand(const std::vector<std::unique_ptr<TH1D>> &histograms,
                                                             Color_t color);
      /// @brief Loop over directory get histograms into vector and add to averaged hist
      void ProcessAutoCorrelationDirectory(TDirectoryFile *autocor_dir,
                                           std::unique_ptr<TH1D>& average_hist,
                                           int &parameter_count,
                                           std::vector<std::unique_ptr<TH1D>> &histograms);

      void ProcessDiagnosticFile(const TString &file_path,
                                 std::unique_ptr<TH1D>& average_hist,
                                 int &parameter_count,
                                 std::vector<std::unique_ptr<TH1D>> &histograms);
      std::unique_ptr<TH1D> AutocorrProcessInputs(const TString &input_file,
                                                  std::vector<std::unique_ptr<TH1D>> &histograms);

      void CompareAverageAC(const std::vector<std::vector<std::unique_ptr<TH1D>>> &histograms,
                            const std::vector<std::unique_ptr<TH1D>> &averages,
                            const std::vector<TString> &hist_labels,
                            const TString &output_name,
                            bool draw_min_max = true,
                            bool draw_all = false,
                            bool draw_errors = true);

      /// @brief Calculate mean AC based on all parameters with error band. Great for comparing AC between different chains
      void PlotAverageACMult(std::vector<TString> input_files,
                             std::vector<TString> hist_labels,
                             const TString &output_name,
                             bool draw_min_max = true,
                             bool draw_all = false,
                             bool draw_errors = true);
  };
}
