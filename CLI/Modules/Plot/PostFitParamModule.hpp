#pragma once
#include "CLI/API/plugin.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "Utils/Plotting/PlottingUtils.h"
#include "Utils/Plotting/PlottingManager.h"

_MaCh3_Safe_Include_Start_ //{
#include "TROOT.h"
#include "TGaxis.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TCandle.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"
_MaCh3_Safe_Include_End_ //}

/// @file GetPostfitParamPlots.hpp
/// This module generates post-fit parameter plots. The central postfit value is
/// taken as the Highest Posterior Density (HPD), but can be easily changed to
/// another method such as Gaussian. Be cautious as parameter names and the number
/// of parameters per plot are currently hardcoded.
///
/// @details
/// Usage:
/// ```
/// ./mach3 plot postfitparam ProcessMCMC_Output1.root <ProcessMCMC_Output2.root> <ProcessMCMC_Output3.root>
/// ```
///
/// @ingroup MaCh3Plotting
///
/// @author Clarence Wret
/// @author Will Parker
/// @author Kamil Skwarczynski
/// @author Ewan Miller

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

namespace M3{

  /// @class GetPostFitParamModule
  /// @brief Module for generating post-fit parameter plots
  ///
  /// This module provides diagnostics such as convergence tests, autocorrelation
  /// analysis, and batch means calculations to assess the quality of MCMC samples.
  class PostFitParamModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~PostFitParamModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;
    private:

        std::vector<std::string> rebuildCommandLine();
        void validateCommandLine();

        void copyParToBlockHist(const int localBin, const std::string& paramName, TH1D* blockHist,
                                const std::string& type, const int fileId, const bool setLabels = true);

        void InitializePads(const TCanvas* canvas, TPad*& pad3, TPad*& pad4);

        void CopyViolinToBlock(TH2D* FullViolin, TH2D* ReducedViolin, const std::vector<std::string>& ParamNames);

        template <typename HistType>
        void PrettifyTitles(HistType* hist);

        bool ReadSettings(const std::shared_ptr<TFile>& File1);
        std::unique_ptr<TH1D> makeRatio(TH1D *PrefitCopy, TH1D *PostfitCopy, bool setAxes);
        bool IsInvalidHist(const TH1D* hist, double invalid = M3::_BAD_DOUBLE_);
        void DrawPlots(TCanvas *plotCanv, TH1D* PrefitCopy, const std::vector<std::unique_ptr<TH1D>>& PostfitVec,
                       TPad *mainPad, TPad *ratioPad, const std::string& OutName);
        void MakeParameterPlots(TCanvas* canv);
        void MakeFluxPlots(TCanvas* canv);
        /// @warning This is legacy functions and will become deprecated
        void MakeNDDetPlots();
        void MakeRidgePlots();
        void GetPostfitParamPlots();
        std::unique_ptr<TGraphAsymmErrors> MakeTGraphAsymmErrors(const std::shared_ptr<TFile>& File,  std::vector<int> Index = {});
        /// @brief KS: Make fancy violin plots
        void GetViolinPlots();

        /// @brief KS: Make comparison of 2D Posteriors
        void Get2DComparison(const std::string& FileName1, const std::string& FileName2);

        Plotting::PlottingManager *PlotMan = nullptr;

        int NDParameters;
        int NDParametersStartingPos;

        std::vector<int> NDSamplesBins;
        std::vector<std::string> NDSamplesNames;

        std::string SaveName;
        std::string plotType;

  };
}
