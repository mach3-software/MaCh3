/// @file DiagMCMCModule.cpp
/// @ingroup MaCh3DiagnosticProcessing
/// @brief Implementation of the DiagMCMCModule class
///
/// @author Clarence Wret
/// @author Kamil Skwarczynski
#include "Fitters/Processing/MCMCProcessor.h"
#include "Manager/Manager.h"
#include "CLI/Modules/DiagMCMCModule.hpp"

namespace M3{

  DiagMCMCModule::~DiagMCMCModule() = default;

  MaCh3ArgumentParser* DiagMCMCModule::get_parser(){
    m_parser = std::make_unique<MaCh3ArgumentParser>("diag", "1.0", argparse::default_arguments::help);
    m_parser->add_description("Tool for MCMC diagnostic like autocorrelations.");
    m_parser->add_argument("mcmc-output")
      .help("MCMC chain root file.")
      .metavar("MCMC_CHAIN")
      .required();
    m_parser->add_argument("config")
      .help("Config file.")
      .metavar("CONFIG")
      .required();
    return m_parser.get();
  }


  /// @brief Main function  creating MCMCProcessor and calling MCMC Diagnostic
  /// @param inputFile MCMC Chain
  /// @param config Config file with settings
  int DiagMCMCModule::Run() {
    SetMaCh3LoggerFormat();
    MACH3LOG_INFO("Producing single fit output");
    std::string inputFile = m_parser->get<std::string>("mcmc-output");
    std::string config = m_parser->get<std::string>("config");

    MACH3LOG_INFO("File for study: {}", inputFile);

    YAML::Node Settings = M3OpenConfig(config);

    // Make the processor
    auto Processor = std::make_unique<MCMCProcessor>(inputFile);
    Processor->SetOutputSuffix("_MCMC_Diag");
    //KS:Turn off plotting detector and some other setting
    Processor->SetExcludedTypes(GetFromManager<std::vector<std::string>>(Settings["DiagMCMC"]["ExcludedTypes"], {}, __FILE__, __LINE__));
    Processor->SetExcludedNames(GetFromManager<std::vector<std::string>>(Settings["DiagMCMC"]["ExcludedNames"], {}, __FILE__, __LINE__));
    Processor->SetExcludedGroups(GetFromManager<std::vector<std::string>>(Settings["DiagMCMC"]["ExcludedGroups"], {}, __FILE__, __LINE__));
    Processor->SetPlotRelativeToPrior(GetFromManager<bool>(Settings["DiagMCMC"]["PlotRelativeToPrior"], false, __FILE__, __LINE__));
    //KS: Use 20 batches for batched means
    Processor->SetnBatches(GetFromManager<int>(Settings["DiagMCMC"]["nBatches"], 20, __FILE__, __LINE__));
    Processor->SetnLags(GetFromManager<int>(Settings["DiagMCMC"]["nLags"], 25000, __FILE__, __LINE__));
    Processor->SetPrintToPDF(GetFromManager<bool>(Settings["PrintToPDF"], true, __FILE__, __LINE__));
    Processor->Initialise();
    if(Settings["MaxEntries"]) {
        Processor->SetEntries(Get<int>(Settings["MaxEntries"], __FILE__, __LINE__));
    }
    //KS: finally call main method
    Processor->DiagMCMC();

    return 0;
  }
}
