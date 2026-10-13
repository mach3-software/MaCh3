#pragma once
#include "CLI/API/plugin.hpp"

//MaCh3 Includes
#include "Utils/Plotting/PlottingUtils.h"
#include "Utils/Plotting/PlottingManager.h"

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file MatrixModule.hpp
/// @ingroup MaCh3Plotting
/// @author Kamil Skwarczynski

namespace M3{

  /// @class MatrixModule
  /// @brief Module for plotting and analyzing matrices
  ///
  /// This module provides functionality for extracting submatrices, setting up
  /// canvases, and plotting and comparing matrices.
  class MatrixModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~MatrixModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;
    private:
    /// @brief Grab large Matrix only only extract submatrix based on label naming
        std::unique_ptr<TH2D> GetSubMatrix(TH2D *MatrixFull,
                                        const std::string& Title,
                                        const std::vector<std::string>& Params,
                                        const std::unique_ptr<Plotting::PlottingManager>& man);
        void SetupCanvas(TCanvas* canv);
        void DynamicLabelSize(TH2D* Hist);
        void SetupInfo(const std::string& Config, std::vector<std::string>& Title, std::vector<std::vector<std::string>>& Params);
        void PlotMatrix(const std::unique_ptr<Plotting::PlottingManager>& man, const std::string& Config, const std::string& File);
        void CompareMatrices(const std::unique_ptr<Plotting::PlottingManager>& man,
                             const std::string& Config, const std::string& File1, const std::string& Title1,
                             const std::string& File2, const std::string& Title2);
        std::unique_ptr<Plotting::PlottingManager> m_plottingManager;
  };
}
