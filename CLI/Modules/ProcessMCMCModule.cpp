/// @file ProcessMCMCModule.cpp
/// @brief Implementation of the ProcessMCMCModule class

//MaCh3 includes
#include "Fitters/OscProcessor.h"
#include "Manager/Manager.h"
#include "CLI/Modules/ProcessMCMCModule.hpp"


namespace M3{

  ProcessMCMCModule::~ProcessMCMCModule() = default;

  MaCh3ArgumentParser* ProcessMCMCModule::get_parser(){
    m_parser = std::make_unique<MaCh3ArgumentParser>("process", "1.0", argparse::default_arguments::help);
    m_parser->add_description("Main exectable responsible for different types of MCMC processing like drawing posteriors, triangle plots etc.");
    m_parser->add_epilog("""\
  ProcessMCMC The main application for analysing the ND280 chain.\n\
  It prints posterior distributions after the burn-in cut and allows comparison of two or three different chains.\n\
  Several options can be configured directly in the app, such as selection, burn-in cut, and whether to plot xsec+flux or only flux.\n\
  \n\
  Additional functionality includes:\n\
    Produce a covariance matrix with multithreading (RAM intensive due to caching)\n\
    Violin plots\n\
    Credible intervals and regions\n\
    Calculation of Bayes factors with significance based on the Jeffreys scale\n\
    Triangle plots\n\
    Study of covariance matrix stability\n""");
    m_parser->add_argument("--corr")
      .help("plot correlation - Same as PlotCorr option in config.")
      .flag();
    m_parser->add_argument("--MakeCredibleIntervals")
      .help("Same as MakeCredibleIntervals option in config.")
      .flag();
    m_parser->add_argument("--CalcBayesFactor")
      .help("Same as CalcBayesFactor option in config.")
      .flag();
    m_parser->add_argument("--CalcSavageDickey")
      .help("Same as CalcSavageDickey option in config.")
      .flag();
    m_parser->add_argument("--CalcBipolarPlot")
      .help("Same as CalcBipolarPlot option in config.")
      .flag();
    m_parser->add_argument("--CalcParameterEvolution")
      .help("Same as CalcParameterEvolution option in config.")
      .flag();
    m_parser->add_argument("config")
      .help("Config file.")
      .metavar("CONFIG")
      .required();
    m_parser->add_argument("mcmc-chain")
      .help("MCMC chain root files and titles.\n"
            "single chain mode: MCMC_CHAIN1\n"
            "two chain mode   : MCMC_CHAIN1 TITLE1 MCMC_CHAIN2 TITLE2\n"
            "three chain mode : MCMC_CHAIN1 TITLE1 MCMC_CHAIN2 TITLE2 MCMC_CHAIN3 TITLE3")
      .metavar("MCMC_CHAIN1 [TITLE1 MCMC_CHAIN2 TITLE2 [MCMC_CHAIN3 TITLE3]]")
      .nargs(1, 6)
      .required();
    return m_parser.get();
  }


  int ProcessMCMCModule::Run() {
    SetMaCh3LoggerFormat();
    nFiles = 0;
    config = m_parser->get<std::string>("config");
    auto mcmc_chain_args = m_parser->get<std::vector<std::string>>("mcmc-chain");

    int nargs = mcmc_chain_args.size();
    if (nargs != 1 && nargs !=4 && nargs != 6)
    {
      MACH3LOG_ERROR("invalid number of arguments: {}", nargs);
      std::cerr << *m_parser;
      MACH3LOG_ERROR("invalid number of arguments: {}", nargs);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    YAML::Node card_yaml = M3OpenConfig(config);
    if (!CheckNodeExists(card_yaml, "ProcessMCMC")) {
      MACH3LOG_ERROR("The 'ProcessMCMC' node is not defined in the YAML configuration.");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    if (mcmc_chain_args.size() == 1)
    {
      MACH3LOG_INFO("Producing single fit output");
      std::string filename = mcmc_chain_args[0];
      this->ProcessMCMC(filename);
    }
    // If we want to compare two or more fits (e.g. binning changes or introducing new params/priors)
    else if (mcmc_chain_args.size() > 1)
    {
      for (std::size_t i = 0; i < mcmc_chain_args.size(); i += 2) {
          FileNames.push_back(mcmc_chain_args[i]);
          TitleNames.push_back(mcmc_chain_args[i + 1]);
      }
      // MACH3LOG_INFO("Producing two fit comparison");
      // FileNames.push_back(files[0]);
      // TitleNames.push_back("ONE"); // todo fix

      // FileNames.push_back(files[1]);
      // TitleNames.push_back("TWO");
      // //KS: If there is third file add it
      // if(files.size() == 3)
      // {
      //   FileNames.push_back(files[2]);
      //   TitleNames.push_back("THREE");
      // }

      this->MultipleProcessMCMC();
    }
    return 0;
  }

  /// @brief Parse custom binning edges from a YAML configuration node.
  ///
  /// This function reads the `CustomBinEdges` section of the YAML config and
  /// returns a mapping of parameter names to their lower and upper edges.
  ///
  /// The expected YAML syntax is:
  /// @code{.yaml}
  /// CustomBinEdges:
  ///   delta_cp: [-3.141592, 3.141592]
  ///   another_param: [min, max]
  /// @endcode
  ///
  /// @param Settings YAML configuration node containing the optional
  ///        `CustomBinEdges` section.
  std::map<std::string, std::pair<double, double>> ProcessMCMCModule::GetCustomBinning(const YAML::Node& Settings)
  {
    std::map<std::string, std::pair<double, double>> CustomBinning;
    if (Settings["CustomBinEdges"]) {
      const YAML::Node& edges = Settings["CustomBinEdges"];

      for (const auto& node : edges) {
        std::string key = node.first.as<std::string>();
        auto values = node.second.as<std::vector<double>>();

        if (values.size() == 2) {
          CustomBinning[key] = std::make_pair(values[0], values[1]);
          MACH3LOG_DEBUG("Adding custom binning {} with {:.4f}, {:.4f}", key, values[0], values[1]);
        } else {
          MACH3LOG_ERROR("Invalid number of values for key: {}", key);
          throw MaCh3Exception(__FILE__ , __LINE__ );
        }
      }
    }
    return CustomBinning;
  }

  /// @brief Process a single MCMC chain file
  ///
  /// Loads the MCMC chain, applies configuration settings, and produces
  /// various diagnostic plots and analyses including posteriors, correlations,
  /// credible intervals, and Bayes factors.
  ///
  /// @param inputFile Path to the MCMC chain ROOT file
  void ProcessMCMCModule::ProcessMCMC(const std::string& inputFile)
  {
    MACH3LOG_INFO("File for study: {} with config  {}", inputFile, config);
    // Make the processor)
    auto Processor = std::make_unique<OscProcessor>(inputFile);

    YAML::Node card_yaml = M3OpenConfig(config.c_str());
    YAML::Node Settings = card_yaml["ProcessMCMC"];
    
    const bool PlotCorr = m_parser->get<bool>("--corr") || GetFromManager<bool>(Settings["PlotCorr"], false, __FILE__, __LINE__);

    Processor->SetExcludedTypes(GetFromManager<std::vector<std::string>>(Settings["ExcludedTypes"], {}, __FILE__, __LINE__));
    Processor->SetExcludedNames(GetFromManager<std::vector<std::string>>(Settings["ExcludedNames"], {}, __FILE__, __LINE__));
    Processor->SetExcludedGroups(GetFromManager<std::vector<std::string>>(Settings["ExcludedGroups"], {}, __FILE__, __LINE__));

    //Apply additional cuts to 1D posterior
    Processor->SetPosterior1DCut(GetFromManager<std::string>(Settings["Posterior1DCut"], "", __FILE__, __LINE__));

    if(PlotCorr) Processor->SetOutputSuffix("_drawCorr");
    //KS:Turn off plotting detector and some other setting, should be via some config
    Processor->SetPlotRelativeToPrior(GetFromManager<bool>(Settings["PlotRelativeToPrior"], false, __FILE__, __LINE__));
    Processor->SetPrintToPDF(GetFromManager<bool>(Settings["PrintToPDF"], true, __FILE__, __LINE__));

    //KS: Whether you want prior error bands for parameters with flat prior or not
    Processor->SetPlotErrorForFlatPrior(GetFromManager<bool>(Settings["PlotErrorForFlatPrior"], true, __FILE__, __LINE__));
    Processor->SetFancyNames(GetFromManager<bool>(Settings["FancyNames"], true, __FILE__, __LINE__));
    Processor->SetPlotBinValue(GetFromManager<bool>(Settings["PlotBinValue"], false, __FILE__, __LINE__));
    //KS: Plot only 2D posteriors with correlations greater than 0.2
    Processor->SetPost2DPlotThreshold(GetFromManager<double>(Settings["Post2DPlotThreshold"], 0.2, __FILE__, __LINE__));

    Processor->Initialise();

    if(Settings["BurnInSteps"]) {
      Processor->SetStepCut(Settings["BurnInSteps"].as<int>());
    } else {
      MACH3LOG_WARN("BurnInSteps not set, defaulting to 20%");
      Processor->SetStepCut(static_cast<int>(Processor->GetnSteps()/5));
    }
    if(Settings["MaxEntries"]) {
      Processor->SetEntries(Get<int>(Settings["MaxEntries"], __FILE__, __LINE__));
    }
    if(Settings["NBins"]) {
      Processor->SetNBins(Get<int>(Settings["NBins"], __FILE__, __LINE__));
    }
    if(Settings["Thinning"])
    {
      if(Settings["Thinning"][0].as<bool>()){
        Processor->ThinMCMC(Settings["Thinning"][1].as<int>());
      }
    }
    // Make the postfit
    Processor->MakePostfit(this->GetCustomBinning(Settings));
    Processor->DrawPostfit();
    //KS: Should set via config whether you want below or not
    if(m_parser->get<bool>("--MakeCredibleIntervals") || GetFromManager<bool>(Settings["MakeCredibleIntervals"], true, __FILE__, __LINE__)) {
      Processor->MakeCredibleIntervals(GetFromManager<std::vector<double>>(Settings["CredibleIntervals"], {0.99, 0.90, 0.68}, __FILE__, __LINE__),
                                      GetFromManager<std::vector<short int>>(Settings["CredibleIntervalsColours"], {436, 430, 422}, __FILE__, __LINE__),
                                      GetFromManager<bool>(Settings["CredibleInSigmas"], false, __FILE__, __LINE__));
    }
    if(m_parser->get<bool>("--CalcBayesFactor") || GetFromManager<bool>(Settings["CalcBayesFactor"], true, __FILE__, __LINE__))  this->CalcBayesFactor(Processor.get());
    if(m_parser->get<bool>("--CalcSavageDickey") || GetFromManager<bool>(Settings["CalcSavageDickey"], true, __FILE__, __LINE__)) this->CalcSavageDickey(Processor.get());
    if(m_parser->get<bool>("--CalcBipolarPlot") || GetFromManager<bool>(Settings["CalcBipolarPlot"], false, __FILE__, __LINE__)) this->CalcBipolarPlot(Processor.get());
    if(m_parser->get<bool>("--CalcParameterEvolution") || GetFromManager<bool>(Settings["CalcParameterEvolution"], false, __FILE__, __LINE__)) this->CalcParameterEvolution(Processor.get());

    if(PlotCorr)
    {
      Processor->SetSmoothing(GetFromManager<bool>(Settings["Smoothing"], true, __FILE__, __LINE__));
      // Make the covariance matrix
      //We have different treatment for multithread
      Processor->CacheSteps();
      //KS: Since we cached let's make fancy violins :)
      if(GetFromManager<bool>(Settings["MakeViolin"], true, __FILE__, __LINE__)) Processor->MakeViolin();
      Processor->MakeCovariance_MP();

      Processor->DrawCovariance();
      if(GetFromManager<bool>(Settings["MakeCovarianceYAML"], true, __FILE__, __LINE__)) Processor->MakeCovarianceYAML(GetFromManager<std::string>(Settings["CovarianceYAMLOutName"], "UpdatedCorrelationMatrix.yaml", __FILE__, __LINE__), GetFromManager<std::string>(Settings["CovarianceYAMLMeansMethod"], "HPD", __FILE__, __LINE__));

      auto const &MakeSubOptimality = Settings["MakeSubOptimality"];
      if(MakeSubOptimality[0].as<bool>()) Processor->MakeSubOptimality(MakeSubOptimality[1].as<int>());

      if(GetFromManager<bool>(Settings["MakeCredibleRegions"], false, __FILE__, __LINE__)) {
        Processor->MakeCredibleRegions(GetFromManager<std::vector<double>>(Settings["CredibleRegions"], {0.99, 0.90, 0.68}, __FILE__, __LINE__),
                                      GetFromManager<std::vector<short int>>(Settings["CredibleRegionStyle"], {2, 1, 3}, __FILE__, __LINE__),
                                      GetFromManager<std::vector<short int>>(Settings["CredibleRegionColor"], {413, 406, 416}, __FILE__, __LINE__),
                                      GetFromManager<bool>(Settings["CredibleInSigmas"], false, __FILE__, __LINE__),
                                      GetFromManager<bool>(Settings["Draw2DPosterior"], true, __FILE__, __LINE__),
                                      GetFromManager<bool>(Settings["DrawBestFit"], true, __FILE__, __LINE__));
      }
      if(GetFromManager<bool>(Settings["GetTrianglePlot"], true, __FILE__, __LINE__)) this->GetTrianglePlot(Processor.get());

      //KS: When creating covariance matrix longest time is spend on caching every step, since we already cached we can run some fancy covariance stability diagnostic
      if(GetFromManager<bool>(Settings["DiagnoseCovarianceMatrix"], false, __FILE__, __LINE__)) this->DiagnoseCovarianceMatrix(Processor.get(), inputFile);
    }
    if(GetFromManager<bool>(Settings["JarlskogAnalysis"], true, __FILE__, __LINE__)) Processor->PerformJarlskogAnalysis();
    if(GetFromManager<bool>(Settings["MakePiePlot"], true, __FILE__, __LINE__))      Processor->MakePiePlot();
  }

  /// @brief Compare multiple MCMC chains
  ///
  /// Processes multiple MCMC chains simultaneously, compares their posterior
  /// distributions, and performs Kolmogorov-Smirnov tests to check if posteriors
  /// are consistent across chains.
  void ProcessMCMCModule::MultipleProcessMCMC()
  {
    YAML::Node card_yaml = M3OpenConfig(config.c_str());
    YAML::Node Settings = card_yaml["ProcessMCMC"];

    constexpr Color_t PosteriorColor[] = {kBlue-1, kRed, kGreen+2};
    //constexpr Style_t PosteriorStyle[] = {kSolid, kDashed, kDotted};
    nFiles = int(FileNames.size());
    std::vector<std::unique_ptr<MCMCProcessor>> Processor(nFiles);

    if(!Settings["BurnInSteps"]) {
      MACH3LOG_WARN("BurnInSteps not set, defaulting to 20%");
    }

    for (int ik = 0; ik < nFiles; ik++)
    {
      MACH3LOG_INFO("File for study: {}", FileNames[ik]);
      // Make the processor
      Processor[ik] = std::make_unique<MCMCProcessor>(FileNames[ik]);
      Processor[ik]->SetOutputSuffix(("_" + std::to_string(ik)).c_str());

      Processor[ik]->SetExcludedTypes(GetFromManager<std::vector<std::string>>(Settings["ExcludedTypes"], {}, __FILE__, __LINE__));
      Processor[ik]->SetExcludedNames(GetFromManager<std::vector<std::string>>(Settings["ExcludedNames"], {}, __FILE__, __LINE__));
      Processor[ik]->SetExcludedGroups(GetFromManager<std::vector<std::string>>(Settings["ExcludedGroups"], {}, __FILE__, __LINE__));

      //Apply additional cuts to 1D posterior
      Processor[ik]->SetPosterior1DCut(GetFromManager<std::string>(Settings["Posterior1DCut"], "", __FILE__, __LINE__));

      Processor[ik]->SetPlotRelativeToPrior(GetFromManager<bool>(Settings["PlotRelativeToPrior"], false, __FILE__, __LINE__));
      Processor[ik]->SetFancyNames(GetFromManager<bool>(Settings["FancyNames"], true, __FILE__, __LINE__));
      Processor[ik]->Initialise();

      if(Settings["BurnInSteps"]) {
        Processor[ik]->SetStepCut(Settings["BurnInSteps"].as<int>());
      }else {
        Processor[ik]->SetStepCut(static_cast<int>(Processor[ik]->GetnSteps()/5));
      }

      if(Settings["MaxEntries"]) {
        Processor[ik]->SetEntries(Get<int>(Settings["MaxEntries"], __FILE__, __LINE__));
      }
      if(Settings["NBins"]) {
        Processor[ik]->SetNBins(Get<int>(Settings["NBins"], __FILE__, __LINE__));
      }
    }

    Processor[0]->MakePostfit(this->GetCustomBinning(Settings));
    Processor[0]->DrawPostfit();
    // Get edges from first histogram to ensure all params use same binning
    std::map<std::string, std::pair<double, double>> ParamEdges;
    for(int i = 0; i < Processor[0]->GetNParams(); ++i) {
      // Get the histogram for the i-th parameter
      TH1D* hist = Processor[0]->GetHpost(i);
      if (!hist) {
        MACH3LOG_DEBUG("Histogram for parameter {} is null.", i);
        continue;
      }

      // Get the parameter name (title of the histogram)
      std::string paramName = hist->GetTitle();

      // Get the axis limits (edges)
      TAxis* axis = hist->GetXaxis();
      double xmin = axis->GetXmin();
      double xmax = axis->GetXmax();

      MACH3LOG_DEBUG("Adding bin edges for {} equal to {:.4f}, {:.4f}",paramName, xmin, xmax);
      // Insert into the map
      ParamEdges[paramName] = std::make_pair(xmin, xmax);
    }

    //KS: Multithreading here is very tempting but there are some issues with root that need to be resovled :(
    for (int ik = 1; ik < nFiles; ik++)
    {
      // Make the postfit
      Processor[ik]->MakePostfit(ParamEdges);
      Processor[ik]->DrawPostfit();
    }

    // Open a TCanvas to write the posterior onto
    auto Posterior = std::make_unique<TCanvas>("PosteriorMulti", "PosteriorMulti", 0, 0, 1024, 1024);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    Posterior->SetGrid();
    Posterior->SetBottomMargin(0.1f);
    Posterior->SetTopMargin(0.05f);
    Posterior->SetRightMargin(0.03f);
    Posterior->SetLeftMargin(0.15f);
    
    // First filename: keep path, just remove ".root"
    // Would be nice to specify outpath in a later update
    size_t pos = FileNames[0].rfind(".root");
    std::string base = (pos == std::string::npos) ? FileNames[0] : FileNames[0].substr(0, pos);
    TString canvasname = base;

    // Remaining filenames: strip path and ".root"
    // So if you have /path/to/file1.root and /path/to/file2.root or /another/path/to/file2.root, canvasname = /path/to/file1_file2.root
    for (int ik = 1; ik < nFiles; ik++) {
      pos = FileNames[ik].find_last_of('/');
      base = (pos == std::string::npos) ? FileNames[ik] : FileNames[ik].substr(pos + 1);
      pos = base.rfind(".root");
      if (pos != std::string::npos) base = base.substr(0, pos);
      canvasname += "_" + TString(base);
    }
    
    canvasname = canvasname +".pdf[";

    Posterior->Print(canvasname);
    // Once the pdf file is open no longer need to bracket
    canvasname.ReplaceAll("[","");

    for(int i = 0; i < Processor[0]->GetNParams(); ++i) 
    {
      // This holds the posterior density
      std::vector<std::unique_ptr<TH1D>> hpost(nFiles);
      std::vector<std::unique_ptr<TLine>> hpd(nFiles);
      hpost[0] = M3::Clone(Processor[0]->GetHpost(i));
      hpost[0]->GetYaxis()->SetTitle("Posterior Density");
      bool Skip = false;
      for (int ik = 1 ; ik < nFiles;  ik++)
      {
        // KS: If somehow this chain doesn't given params we skip it
        const int Index = Processor[ik]->GetParamIndexFromName(hpost[0]->GetTitle());
        if(Index == M3::_BAD_INT_)
        {
          Skip = true;
          break;
        }
        hpost[ik] = M3::Clone(Processor[ik]->GetHpost(Index));
      }

      // Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
      if(hpost[0]->GetMaximum() == hpost[0]->Integral()*1.5 || Skip) {
        continue;
      }
      for (int ik = 0; ik < nFiles;  ik++)
      {
        RemoveFitter(hpost[ik].get(), "Gauss");

        // Set some nice colours
        hpost[ik]->SetLineColor(PosteriorColor[ik]);
        //hpost[ik]->SetLineStyle(PosteriorStyle[ik]);
        hpost[ik]->SetLineWidth(2);

        // Area normalise the distributions
        hpost[ik]->Scale(1./hpost[ik]->Integral());
      }
      TString Title;
      double Prior = 1.0;
      double PriorError = 1.0;

      Processor[0]->GetNthParameter(i, Prior, PriorError, Title);

      // Now make the TLine for the Asimov
      auto Asimov = std::make_unique<TLine>(Prior, hpost[0]->GetMinimum(), Prior, hpost[0]->GetMaximum());
      Asimov->SetLineColor(kRed-3);
      Asimov->SetLineWidth(2);
      Asimov->SetLineStyle(kDashed);

      // Make a nice little TLegend
      auto leg = std::make_unique<TLegend>(0.20, 0.7, 0.6, 0.97);
      leg->SetTextSize(0.03f);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetLineColor(0);
      leg->SetLineStyle(0);
      TString asimovLeg = Form("#splitline{Prior}{x = %.2f , #sigma = %.2f}", Prior, PriorError);
      leg->AddEntry(Asimov.get(), asimovLeg, "l");

      for (int ik = 0; ik < nFiles; ik++)
      {
        TString rebinLeg = Form("#splitline{%s}{#mu = %.2f, #sigma = %.2f}", TitleNames[ik].c_str(), hpost[ik]->GetMean(), hpost[ik]->GetRMS());
        leg->AddEntry(hpost[ik].get(),  rebinLeg, "l");

        hpd[ik] = std::make_unique<TLine>(hpost[ik]->GetBinCenter(hpost[ik]->GetMaximumBin()), hpost[ik]->GetMinimum(),
                                          hpost[ik]->GetBinCenter(hpost[ik]->GetMaximumBin()), hpost[ik]->GetMaximum());
        hpd[ik]->SetLineColor(hpost[ik]->GetLineColor());
        hpd[ik]->SetLineWidth(2);
        hpd[ik]->SetLineStyle(kSolid);
      }

      // Find the maximum value to nicely resize hist
      double maximum = 0;
      for (int ik = 0; ik < nFiles; ik++) maximum = std::max(maximum, hpost[ik]->GetMaximum());
      for (int ik = 0; ik < nFiles; ik++) hpost[ik]->SetMaximum(1.3*maximum);

      hpost[0]->Draw("hist");
      for (int ik = 1; ik < nFiles; ik++) hpost[ik]->Draw("hist same");
      Asimov->Draw("same");
      for (int ik = 0; ik < nFiles; ik++) hpd[ik]->Draw("same");
      leg->Draw("same");
      Posterior->cd();
      Posterior->Print(canvasname);
    }//End loop over parameters
      
    // Finally draw the parameter plot onto the PDF
    // Close the .pdf file with all the posteriors
    Posterior->cd();
    Posterior->Clear();

    if(GetFromManager<bool>(Settings["PerformKStest"], true, __FILE__, __LINE__)) this->KolmogorovSmirnovTest(Processor, Posterior, canvasname);
    
    // Close the pdf file
    MACH3LOG_INFO("Closing pdf {}", canvasname);
    canvasname+="]";
    Posterior->Print(canvasname);
  }

  /// @brief Calculate Bayes factors for hypothesis testing
  ///
  /// Computes Bayes factors to compare different models or hypotheses,
  /// particularly useful for oscillation parameters and switch parameters.
  /// Configuration is read from the YAML file under BayesFactor section.
  ///
  /// @param Processor Pointer to the MCMCProcessor instance
  void ProcessMCMCModule::CalcBayesFactor(MCMCProcessor* Processor)
  {
    YAML::Node card_yaml = M3OpenConfig(config.c_str());
    YAML::Node Settings = card_yaml["ProcessMCMC"];

    std::vector<std::string> ParNames;
    std::vector<std::vector<double>> Model1Bounds;
    std::vector<std::vector<double>> Model2Bounds;
    std::vector<std::vector<std::string>> ModelNames;
    for (const auto& dg : Settings["BayesFactor"])
    {
      ParNames.push_back(dg[0].as<std::string>());
      ModelNames.push_back(dg[1].as<std::vector<std::string>>());
      Model1Bounds.push_back(dg[2].as<std::vector<double>>());
      Model2Bounds.push_back(dg[3].as<std::vector<double>>());
    }

    Processor->GetBayesFactor(ParNames, Model1Bounds, Model2Bounds, ModelNames);
  }

  /// @brief Calculate Savage-Dickey ratios
  ///
  /// Computes Savage-Dickey ratios for Bayes factor estimation at specific
  /// parameter values. Configuration is read from the YAML file under
  /// SavageDickey section.
  ///
  /// @param Processor Pointer to the MCMCProcessor instance
  void ProcessMCMCModule::CalcSavageDickey(MCMCProcessor* Processor)
  {
    YAML::Node card_yaml = M3OpenConfig(config.c_str());
    YAML::Node Settings = card_yaml["ProcessMCMC"];

    std::vector<std::string> ParNames;
    std::vector<double> EvaluationPoint;
    std::vector<std::vector<double>> Bounds;
    
    for (const auto& d : Settings["SavageDickey"])
    {
      ParNames.push_back(d[0].as<std::string>());
      EvaluationPoint.push_back(d[1].as<double>());
      Bounds.push_back(d[2].as<std::vector<double>>());
    }
    Processor->GetSavageDickey(ParNames, EvaluationPoint, Bounds);
  }

  /// @brief Calculate parameter evolution over MCMC steps
  ///
  /// Tracks how parameters evolve during the MCMC chain, useful for
  /// diagnosing convergence and burn-in. Configuration is read from
  /// the YAML file under ParameterEvolution section.
  ///
  /// @param Processor Pointer to the MCMCProcessor instance
  void ProcessMCMCModule::CalcParameterEvolution(MCMCProcessor* Processor)
  {
    YAML::Node card_yaml = M3OpenConfig(config.c_str());
    YAML::Node Settings = card_yaml["ProcessMCMC"];

    std::vector<std::string> ParNames;
    std::vector<int> Intervals;
    for (const auto& d : Settings["ParameterEvolution"])
    {
      ParNames.push_back(d[0].as<std::string>());
      Intervals.push_back(d[1].as<int>());
    }
    Processor->ParameterEvolution(ParNames, Intervals);
  }

  /// @brief Create bipolar plots for parameter visualization
  ///
  /// Generates bipolar plots to visualize parameter distributions in a
  /// circular representation. Configuration is read from the YAML file
  /// under BipolarPlot section.
  ///
  /// @param Processor Pointer to the MCMCProcessor instance
  void ProcessMCMCModule::CalcBipolarPlot(MCMCProcessor* Processor)
  {
    YAML::Node card_yaml = M3OpenConfig(config.c_str());
    YAML::Node Settings = card_yaml["ProcessMCMC"];

    std::vector<std::string> ParNames;
    for (const auto& d : Settings["BipolarPlot"])
    {
      ParNames.push_back(d[0].as<std::string>());
    }
    Processor->GetPolarPlot(ParNames);
  }

  /// @brief Generate triangle plots showing parameter correlations
  ///
  /// Creates triangle plots displaying 1D and 2D posterior distributions
  /// for sets of correlated parameters. Configuration is read from the
  /// YAML file under TrianglePlot section.
  ///
  /// @param Processor Pointer to the MCMCProcessor instance
  void ProcessMCMCModule::GetTrianglePlot(MCMCProcessor* Processor) {
    YAML::Node card_yaml = M3OpenConfig(config.c_str());
    YAML::Node Settings = card_yaml["ProcessMCMC"];

    for (const auto& dg : Settings["TrianglePlot"])
    {
      std::string ParName = dg[0].as<std::string>();

      std::vector<std::string> NameVec = dg[1].as<std::vector<std::string>>();
      Processor->MakeTrianglePlot(NameVec,
                                  GetFromManager<std::vector<double>>(Settings["CredibleIntervals"], {0.99, 0.90, 0.68}, __FILE__, __LINE__),
                                  GetFromManager<std::vector<short int>>(Settings["CredibleIntervalsColours"], {436, 430, 422}, __FILE__, __LINE__),
                                  GetFromManager<std::vector<double>>(Settings["CredibleRegions"], {0.99, 0.90, 0.68}, __FILE__, __LINE__),
                                  GetFromManager<std::vector<short int>>(Settings["CredibleRegionStyle"], {2, 1, 3}, __FILE__, __LINE__),
                                  GetFromManager<std::vector<short int>>(Settings["CredibleRegionColor"], {413, 406, 416}, __FILE__, __LINE__),
                                  GetFromManager<bool>(Settings["CredibleInSigmas"], false, __FILE__, __LINE__));
    }
  }

  /// @brief Diagnose covariance matrix stability
  ///
  /// Validates the stability of the posterior covariance matrix by computing
  /// it at different burn-in cuts and comparing successive matrices. This helps
  /// determine when the matrix becomes stable and the appropriate burn-in length.
  ///
  /// @param Processor Pointer to the MCMCProcessor instance
  /// @param inputFile Path to the input file for naming output files
  void ProcessMCMCModule::DiagnoseCovarianceMatrix(MCMCProcessor* Processor, const std::string& inputFile)
  {
    //Turn of plots from Processor
    Processor->SetPrintToPDF(false);
    // Open a TCanvas to write the posterior onto
    auto Canvas = std::make_unique<TCanvas>("Canvas", "Canvas", 0, 0, 1024, 1024);
    Canvas->SetGrid();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    Canvas->SetTickx();
    Canvas->SetTicky();
    Canvas->SetBottomMargin(0.1f);
    Canvas->SetTopMargin(0.05f);
    Canvas->SetRightMargin(0.15f);
    Canvas->SetLeftMargin(0.10f);
    
    //KS: Fancy colours
    const int NRGBs = 10;
    TColor::InitializeColors();
    Double_t stops[NRGBs] = { 0.00, 0.10, 0.25, 0.35, 0.50, 0.60, 0.65, 0.75, 0.90, 1.00 };
    Double_t red[NRGBs]   = { 0.50, 1.00, 1.00, 0.25, 0.00, 0.10, 0.50, 1.00, 0.75, 0.55 };
    Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00, 0.60, 0.90, 1.00, 0.75, 0.75 };
    Double_t blue[NRGBs]  = { 0.00, 0.25, 1.00, 1.00, 0.50, 0.60, 0.90, 1.00, 0.05, 0.05 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 255);
    gStyle->SetNumberContours(255);
    
    std::string OutName = inputFile;
    OutName = OutName.substr(0, OutName.find(".root"));
    Canvas->Print(Form("Correlation_%s.pdf[", OutName.c_str()), "pdf");
    Canvas->Print(Form("Covariance_%s.pdf[", OutName.c_str()), "pdf");
    
    YAML::Node card_yaml = M3OpenConfig(config.c_str());
    YAML::Node Settings = card_yaml["ProcessMCMC"];

    const int entries = int(Processor->GetnSteps());
    const int NIntervals = GetFromManager<int>(Settings["NIntervals"], 5, __FILE__, __LINE__);
    const int IntervalsSize = entries/NIntervals;
    //We start with burn from 0 (no burn in at all)
    int BurnIn = 0;
    MACH3LOG_INFO("Diagnosing matrices with entries={}, NIntervals={} and IntervalsSize={}", entries, NIntervals, IntervalsSize);

    TMatrixDSym *Covariance = nullptr;
    TMatrixDSym *Correlation = nullptr;
    
    TH2D *CovariancePreviousHist = nullptr;
    TH2D *CorrelationPreviousHist = nullptr;

    TH2D *CovarianceHist = nullptr;
    TH2D *CorrelationHist = nullptr;

    //KS: Get first covariances, we need two for comparison...
    Processor->SetStepCut(BurnIn);
    Processor->GetCovariance(Covariance, Correlation);
    
    CovariancePreviousHist = this->TMatrixIntoTH2D(Covariance, "Covariance"); 
    CorrelationPreviousHist = this->TMatrixIntoTH2D(Correlation, "Correlation");
        
    delete Covariance;
    Covariance = nullptr;
    delete Correlation;
    Correlation = nullptr;
    
    //KS: Loop over all desired cuts
    for(int k = 1; k < NIntervals; ++k)
    {
      BurnIn = k*IntervalsSize;
      Processor->SetStepCut(BurnIn);
      Processor->GetCovariance(Covariance, Correlation);
      Processor->Reset2DPosteriors();
      
      CovarianceHist = this->TMatrixIntoTH2D(Covariance, "Covariance"); 
      CorrelationHist = this->TMatrixIntoTH2D(Correlation, "Correlation");
              
      TH2D *CovarianceDiff = static_cast<TH2D*>(CovarianceHist->Clone("Covariance_Ratio"));
      TH2D *CorrelationDiff = static_cast<TH2D*>(CorrelationHist->Clone("Correlation_Ratio"));
      
      //KS: Bit messy but quite often covariance is 0 is divided by 0 is problematic so
      #ifdef MULTITHREAD
      #pragma omp parallel for
      #endif
      for (int j = 1; j < CovarianceDiff->GetXaxis()->GetNbins()+1; ++j)
      {
        for (int i = 1; i < CovarianceDiff->GetYaxis()->GetNbins()+1; ++i)
        {
          if( std::fabs (CovarianceDiff->GetBinContent(j, i)) < 1.e-5 && std::fabs (CovariancePreviousHist->GetBinContent(j, i)) < 1.e-5)
          {
            CovarianceDiff->SetBinContent(j, i, M3::_BAD_DOUBLE_);
            CovariancePreviousHist->SetBinContent(j, i, M3::_BAD_DOUBLE_);
          }
          if( std::fabs (CorrelationDiff->GetBinContent(j, i)) < 1.e-5 && std::fabs (CorrelationPreviousHist->GetBinContent(j, i)) < 1.e-5)
          {
            CorrelationDiff->SetBinContent(j, i, M3::_BAD_DOUBLE_);
            CorrelationPreviousHist->SetBinContent(j, i, M3::_BAD_DOUBLE_);
          }
        }
      }
      //Divide matrices
      CovarianceDiff->Divide(CovariancePreviousHist);
      CorrelationDiff->Divide(CorrelationPreviousHist);
      
      //Now it is time for fancy names etc.
      for (int j = 0; j < CovarianceDiff->GetXaxis()->GetNbins(); ++j)
      {
        TString Title = "";
        double Prior = 1.0;
        double PriorError = 1.0;
    
        Processor->GetNthParameter(j, Prior, PriorError, Title);
        
        CovarianceDiff->GetXaxis()->SetBinLabel(j+1, Title);
        CovarianceDiff->GetYaxis()->SetBinLabel(j+1, Title);
        CorrelationDiff->GetXaxis()->SetBinLabel(j+1, Title);
        CorrelationDiff->GetYaxis()->SetBinLabel(j+1, Title);
      }
      CovarianceDiff->GetXaxis()->SetLabelSize(0.015f);
      CovarianceDiff->GetYaxis()->SetLabelSize(0.015f);
      CorrelationDiff->GetXaxis()->SetLabelSize(0.015f);
      CorrelationDiff->GetYaxis()->SetLabelSize(0.015f);
      
      std::stringstream ss;
      ss << "BCut_";
      ss << BurnIn;
      ss << "/";
      ss << "BCut_";
      ss << (k-1)*IntervalsSize;
      std::string str = ss.str();
      
      TString Title = "Cov " + str;
      CovarianceDiff->GetZaxis()->SetTitle( Title );
      Title = "Corr " + str;
      CorrelationDiff->GetZaxis()->SetTitle(Title);
      
      CovarianceDiff->SetMinimum(-2);
      CovarianceDiff->SetMaximum(2);
      CorrelationDiff->SetMinimum(-2);
      CorrelationDiff->SetMaximum(2);

      Canvas->cd();
      CovarianceDiff->Draw("colz");
      Canvas->Print(Form("Covariance_%s.pdf", OutName.c_str()), "pdf");

      CorrelationDiff->Draw("colz");
      Canvas->Print(Form("Correlation_%s.pdf", OutName.c_str()), "pdf");
          
      //KS: Current hist become previous as we need it for further comparison
      delete CovariancePreviousHist;
      CovariancePreviousHist = static_cast<TH2D*>(CovarianceHist->Clone());
      delete CorrelationPreviousHist;
      CorrelationPreviousHist = static_cast<TH2D*>(CorrelationHist->Clone());
      
      delete CovarianceHist;
      CovarianceHist = nullptr;
      delete CorrelationHist;
      CorrelationHist = nullptr;
      
      delete CovarianceDiff;
      delete CorrelationDiff;
      delete Covariance;
      Covariance = nullptr;
      delete Correlation;
      Correlation = nullptr;
    }
    Canvas->cd();
    Canvas->Print(Form("Covariance_%s.pdf]", OutName.c_str()), "pdf");
    Canvas->Print(Form("Correlation_%s.pdf]", OutName.c_str()), "pdf");
    
    Processor->SetPrintToPDF(true);
    if(Covariance != nullptr)              delete Covariance;
    if(Correlation != nullptr)             delete Correlation;
    if(CovariancePreviousHist != nullptr)  delete CovariancePreviousHist;
    if(CorrelationPreviousHist != nullptr) delete CorrelationPreviousHist;
    if(CovarianceHist != nullptr)          delete CovarianceHist;
    if(CorrelationHist != nullptr)         delete CorrelationHist;
  }

  /// @brief Convert TMatrixDSym to TH2D histogram
  ///
  /// Converts a ROOT symmetric matrix to a 2D histogram for visualization
  /// and easier manipulation with ROOT plotting tools.
  ///
  /// @param Matrix Pointer to the symmetric matrix
  /// @param title Title for the resulting histogram
  /// @return Pointer to the newly created TH2D histogram
  TH2D* ProcessMCMCModule::TMatrixIntoTH2D(TMatrixDSym* Matrix, const std::string& title)
  {
    TH2D* hMatrix = new TH2D(title.c_str(), title.c_str(), Matrix->GetNrows(), 0.0, Matrix->GetNrows(), Matrix->GetNcols(), 0.0, Matrix->GetNcols());
    for(int i = 0; i < Matrix->GetNrows(); i++)
    {
      for(int j = 0; j < Matrix->GetNcols(); j++)
      {
        //KS: +1 because there is offset in histogram relative to TMatrix
        hMatrix->SetBinContent(i+1,j+1, (*Matrix)(i,j));
      }
    }
    return hMatrix;
  }

  /// @brief Perform Kolmogorov-Smirnov test between posterior distributions
  ///
  /// Tests whether posterior distributions from different MCMC chains are
  /// consistent by computing the KS test statistic for each parameter.
  /// Results are visualized showing cumulative distributions and D-statistic.
  ///
  /// @param Processor Vector of MCMCProcessor instances (one per chain)
  /// @param Posterior Canvas for drawing the comparison plots
  /// @param canvasname Name for the output PDF file
  void ProcessMCMCModule::KolmogorovSmirnovTest(const std::vector<std::unique_ptr<MCMCProcessor>>& Processor,
                                          const std::unique_ptr<TCanvas>& Posterior,
                                          const TString& canvasname)
  {
    constexpr Color_t CumulativeColor[] = {kBlue-1, kRed, kGreen+2};
    constexpr Style_t CumulativeStyle[] = {kSolid, kDashed, kDotted};

    for(int i = 0; i < Processor[0]->GetNParams(); ++i) 
    {
      // This holds the posterior density
      std::vector<std::unique_ptr<TH1D>> hpost(nFiles);
      std::vector<std::unique_ptr<TH1D>> CumulativeDistribution(nFiles);
        
      TString Title;
      double Prior = 1.0;
      double PriorError = 1.0;

      Processor[0]->GetNthParameter(i, Prior, PriorError, Title);
      bool Skip = false;
      for (int ik = 0 ; ik < nFiles;  ik++)
      {
        int Index = 0;
        if(ik == 0 ) Index = i;
        else
        {
          // KS: If somehow this chain doesn't given params we skip it
          Index = Processor[ik]->GetParamIndexFromName(hpost[0]->GetTitle());
          if(Index == M3::_BAD_INT_)
          {
            Skip = true;
            break;
          }
        }
        hpost[ik] = M3::Clone(Processor[ik]->GetHpost(Index));
        CumulativeDistribution[ik] = M3::Clone(Processor[ik]->GetHpost(Index));
        CumulativeDistribution[ik]->Fill(0., 0.);
        CumulativeDistribution[ik]->Reset();
        CumulativeDistribution[ik]->SetMaximum(1.);
        TString TempTitle = Title+" Kolmogorov Smirnov";
        CumulativeDistribution[ik]->SetTitle(TempTitle);
        
        TempTitle = Title+" Value";
        CumulativeDistribution[ik]->GetXaxis()->SetTitle(TempTitle);
        CumulativeDistribution[ik]->GetYaxis()->SetTitle("Cumulative Probability");
        
        CumulativeDistribution[ik]->SetLineWidth(2);
        CumulativeDistribution[ik]->SetLineColor(CumulativeColor[ik]);
        CumulativeDistribution[ik]->SetLineStyle(CumulativeStyle[ik]);
      }

      // Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
      if(hpost[0]->GetMaximum() == hpost[0]->Integral()*1.5 || Skip) {
        continue;
      }
      
      for (int ik = 0 ; ik < nFiles;  ik++)
      {
        const int NumberOfBins = hpost[ik]->GetXaxis()->GetNbins();
        double Cumulative = 0;
        const double Integral = hpost[ik]->Integral();
        for (int j = 1; j < NumberOfBins+1; ++j)
        {
          Cumulative += hpost[ik]->GetBinContent(j)/Integral;
          CumulativeDistribution[ik]->SetBinContent(j, Cumulative);
        }
        //KS: Set overflow to 1 just in case
        CumulativeDistribution[ik]->SetBinContent(NumberOfBins+1, 1.);
      }
      
      std::vector<int> TestStatBin(nFiles, 0);
      std::vector<double> TestStatD(nFiles, -999);
      std::vector<std::unique_ptr<TLine>> LineD(nFiles);
      //Find KS statistic
      for (int ik = 1 ; ik < nFiles;  ik++)
      {
        const int NumberOfBins = CumulativeDistribution[0]->GetXaxis()->GetNbins();
        for (int j = 1; j < NumberOfBins+1; ++j)
        {
          const double BinValue = CumulativeDistribution[0]->GetBinCenter(j);
          const int BinNumber = CumulativeDistribution[ik]->FindBin(BinValue);
          //KS: Calculate D statistic for this bin, only save it if it's bigger than previously found value
          double TempDstat = std::fabs(CumulativeDistribution[0]->GetBinContent(j) - CumulativeDistribution[ik]->GetBinContent(BinNumber));
          if(TempDstat > TestStatD[ik])
          {
            TestStatD[ik] = TempDstat;
            TestStatBin[ik] = j;
          }
        }
      }

      for (int ik = 0 ; ik < nFiles;  ik++)
      {
        LineD[ik] = std::make_unique<TLine>(CumulativeDistribution[0]->GetBinCenter(TestStatBin[ik]), 0, CumulativeDistribution[0]->GetBinCenter(TestStatBin[ik]), CumulativeDistribution[0]->GetBinContent(TestStatBin[ik]));
        LineD[ik]->SetLineColor(CumulativeColor[ik]);
        LineD[ik]->SetLineWidth(2.0);
      }
      CumulativeDistribution[0]->Draw();
      for (int ik = 0 ; ik < nFiles;  ik++)
        CumulativeDistribution[ik]->Draw("SAME");
      
      auto leg = std::make_unique<TLegend>(0.15, 0.7, 0.5, 0.90);
      leg->SetTextSize(0.04f);
      for (int ik = 0; ik < nFiles;  ik++)
        leg->AddEntry(CumulativeDistribution[ik].get(), TitleNames[ik].c_str(), "l");
      for (int ik = 1; ik < nFiles;  ik++)
        leg->AddEntry(LineD[ik].get(), Form("#Delta D = %.4f", TestStatD[ik]), "l");
      
      leg->SetLineColor(0);
      leg->SetLineStyle(0);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->Draw("SAME");

      for (int ik = 1; ik < nFiles;  ik++)
        LineD[ik]->Draw("sam");
        
      Posterior->cd();
      Posterior->Print(canvasname);
    } //End loop over parameter
  }
};
