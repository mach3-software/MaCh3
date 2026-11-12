#include "Fitters/Processing/MCMCProcessor.h"
#include "Samples/HistogramUtils.h"
#include "Manager/Manager.h"
#include "CLI/Modules/SmearModule.hpp"

/// @file SmearModule.cpp
/// @brief Allows you to smear contour. For example after performing sets of study one finds out that used sets of uncertainty doesn't fully cover analysis need. Then one can smear additionally contour.
/// @ingroup MaCh3DiagnosticProcessing


namespace M3{
    
    SmearModule::~SmearModule() = default;
    
    MaCh3ArgumentParser* SmearModule::get_parser(){
        m_parser = std::make_unique<MaCh3ArgumentParser>("smear", "1.0", argparse::default_arguments::help);
        m_parser->add_description("Tool for smearing MCMC chains.");
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
    
    void SmearModule::SmearChain(const std::string& inputFile, const std::string& config){
        MACH3LOG_INFO("File for study: {}", inputFile);
        
        YAML::Node Settings = M3OpenConfig(config);
        
        // Make the processor
        auto Processor = std::make_unique<MCMCProcessor>(inputFile);
        Processor->SetOutputSuffix("_Smear_MCMC");
        Processor->Initialise();
        
        const auto& Smear = Settings["SmearChain"];
        
        std::vector<std::string> Names = Get<std::vector<std::string>>(Smear["Smear"][0], __FILE__, __LINE__);
        std::vector<double> ErrorValue = Get<std::vector<double>>(Smear["Smear"][1], __FILE__, __LINE__);
        
        bool SaveUnsmearedBranch = GetFromManager<bool>(Smear["SaveUnsmearedBranch"], false, __FILE__ , __LINE__);
        Processor->SmearChain(Names, ErrorValue, SaveUnsmearedBranch);
    }
    
    int SmearModule::Run() {
        SetMaCh3LoggerFormat();
        MACH3LOG_INFO("Producing single fit output");
        std::string filename = m_parser->get<std::string>("mcmc-output");
        std::string config = m_parser->get<std::string>("config");
        SmearChain(filename, config);
        
        return 0;
    }
}


















void SmearChain(const std::string& inputFile, const std::string& config)
{
    MACH3LOG_INFO("File for study: {}", inputFile);
    
    YAML::Node Settings = M3OpenConfig(config);
    
    // Make the processor
    auto Processor = std::make_unique<MCMCProcessor>(inputFile);
    Processor->SetOutputSuffix("_Smear_MCMC");
    Processor->Initialise();
    
    const auto& Smear = Settings["SmearChain"];
    
    std::vector<std::string> Names = Get<std::vector<std::string>>(Smear["Smear"][0], __FILE__, __LINE__);
    std::vector<double> ErrorValue = Get<std::vector<double>>(Smear["Smear"][1], __FILE__, __LINE__);
    
    bool SaveUnsmearedBranch = GetFromManager<bool>(Smear["SaveUnsmearedBranch"], false, __FILE__ , __LINE__);
    Processor->SmearChain(Names, ErrorValue, SaveUnsmearedBranch);
}

int main(int argc, char *argv[]) {
    SetMaCh3LoggerFormat();
    if (argc != 3)
    {
        MACH3LOG_ERROR("How to use: {} MCMC_Output.root config", argv[0]);
        throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    MACH3LOG_INFO("Producing single fit output");
    std::string filename = argv[1];
    std::string config = argv[2];
    SmearChain(filename, config);
    
    return 0;
}
