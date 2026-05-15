//MaCh3 includes
#pragma once
// yaml Includes
#include "yaml-cpp/yaml.h"
#include "Fitters/MCMCProcessor.h"
#include "api/plugin.hpp"

namespace M3{

/// @file ProcessMCMC.cpp
/// @brief Main exectable responsible for different types of MCMC processing like drawing posteriors, triangle plots etc. Actual implantation of methods is in MCMCProcessor
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author Kamil Skwarczynski

/// @brief Main function processing MCMC and Producing plots

  class ProcessMCMCPlugin: public IPlugin{

    public:
      virtual ~ProcessMCMCPlugin();
      MaCh3ArgumentParser* get_parser() override;
      int run() override;


    private:
      std::map<std::string, std::pair<double, double>> GetCustomBinning(const YAML::Node& Settings);
      void ProcessMCMC(const std::string& inputFile);
      /// @brief Function producing comparison of posterior and more betwen a few MCMC chains
      void MultipleProcessMCMC();
      void CalcBayesFactor(MCMCProcessor* Processor);
      void CalcSavageDickey(MCMCProcessor* Processor);
      void CalcParameterEvolution(MCMCProcessor* Processor);
      void CalcBipolarPlot(MCMCProcessor* Processor);
      void GetTrianglePlot(MCMCProcessor* Processor);
      void DiagnoseCovarianceMatrix(MCMCProcessor* Processor, const std::string& inputFile);
      /// @brief KS: Convert TMatrix to TH2D, mostly useful for making fancy plots
      TH2D* TMatrixIntoTH2D(TMatrixDSym* Matrix, const std::string& title);
      /// @brief KS: Perform KS test to check if two posteriors for the same parameter came from the same distribution
      void KolmogorovSmirnovTest(const std::vector<std::unique_ptr<MCMCProcessor>>& Processor,
                                 const std::unique_ptr<TCanvas>& Posterior,
                                 const TString& canvasname);
      int nFiles;
      std::vector <std::string> FileNames;
      std::vector <std::string> TitleNames;
      std::string config;
      MaCh3ArgumentParser* m_parser;
  };
};
