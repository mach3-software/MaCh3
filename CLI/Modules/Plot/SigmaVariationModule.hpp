#pragma once
#include "CLI/API/plugin.hpp"

//MaCh3 Includes
#include "Utils/Plotting/PlottingUtils.h"
#include "Utils/Plotting/PlottingManager.h"

#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file Sigmavariation.hpp
/// @ingroup MaCh3Plotting
/// @author Kamil Skwarczynski
/// @todo add maybe ratio for PlotSigVar1D

namespace M3{

  /// @class SigmaVariationModule
  /// @brief Module for performing sigma variation analysis
  ///
  /// This module provides functionality for analyzing sigma variations in MCMC samples.
  class SigmaVariationModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~SigmaVariationModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;
    private:

        /// @brief Histograms have name like ND_CC0pi_1DProj0_Norm_Param_0_sig_n3.00_val_0.25. This code is trying to extract sigma names
        void FindKnot(std::vector<double>& SigmaValues,
                    const std::string& dirname,
                    const std::string& subdirname,
                    const std::string& ProjName,
                    std::string histname);

        std::unique_ptr<TLegend> MakeLegend(double x1, double y1, double x2, double y2,
                                            double textSize = 0.04);
        /// @brief Scan inputs to figure out dial name and used sample names
        void ScanInput(std::vector<std::string>& DialNameVecr,
                    std::vector<std::string>& SampleNameVec,
                    std::vector<int>& SampleDimVec,
                    std::vector<double>& SigmaValues,
                    const std::string& filename);

        /// @brief Check whether to skip directory or not based on defined strings
        bool SkipDirectory(const std::vector<std::string>& ExcludeString, const std::vector<std::string>& IncludeString, const std::string& dirname);

        /// @brief Extracts dial value for from histogram title
        std::vector<double> GetDialValues(const std::vector<std::unique_ptr<TH1D>>& Poly);
        void InitializePads(TCanvas* canv, TPad*& pad1, TPad*& pad2, double Pad1Bottom = 0.25, double Pad2Top = 0.25);
        void MakeRatio(const std::vector<std::unique_ptr<TH1D>>& Poly,
                       std::vector<std::unique_ptr<TH1D>>& Ratio);
        void PlotRatio(const std::vector<std::unique_ptr<TH1D>>& Poly,
                       const std::unique_ptr<TCanvas>& canv,
                       const std::string& Title,
                       const std::string& outfilename);
        void CompareSigVar1D(const std::string& filename, const YAML::Node& Settings);
        void PlotRatio2D(const std::vector<std::unique_ptr<TH2>>& Poly,
                         const std::unique_ptr<TCanvas>& canv,
                         const std::string& Title,
                         const std::string& outfilename);
        void CompareSigVar2D(const std::string& filename, const YAML::Node& Settings);
        void PlotEventRate(const std::vector<std::vector<std::unique_ptr<TH1D>>>& Poly,
                           const std::unique_ptr<TCanvas>& canv,
                           const std::string& Title,
                           const std::string& outfilename);
        void MakeEventRatePlot(const std::string& filename, const YAML::Node& Settings);
        void PlotSigVar1D(const std::vector<std::vector<std::unique_ptr<TH1D>>>& Projection,
                          const std::unique_ptr<TCanvas>& canv,
                          const std::string& Title,
                          const std::string& outfilename,
                          const std::vector<std::string>& ParamNames,
                          const std::vector<int>& ParamColour);
        void OverlaySigVar1D(const std::string& filename, const YAML::Node& Settings);
        std::vector<std::string> DialNameVector;
        std::vector<std::string> SampleNameVector;
        std::vector<int> SampleMaxDim;
        std::vector<double> sigmaArray;

        int PriorKnot = _BAD_INT_;

        static constexpr int NVars = 5;
        static constexpr Color_t Colours[NVars] = {kRed, kGreen+1, kBlack, kBlue+1, kOrange+1};
        static constexpr ELineStyle Style[NVars] = {kDotted, kDashed, kSolid, kDashDotted, kDashDotted};

        /// @warning KS: keep raw pointer or ensure manual delete of PlotMan. If spdlog in automatically deleted before PlotMan then destructor has some spdlog and this could cause segfault
        Plotting::PlottingManager* PlotMan = nullptr;

  };
}
