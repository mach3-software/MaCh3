#include "Fitters/MCMCProcessor.h"

#include "Manager/Manager.h"

/// @file SmearChain.cpp
/// @brief Allows you to smear contour. For example after performing sets of study one finds out that used sets of uncertainty doesn't fully cover analysis need. Then one can smear additionally contour.

/// @brief Main function  creating MCMCProcessor and calling Smear Chain
/// @param inputFile MCMC Chain
/// @param config Config file with settings
void DiagMCMC(const std::string& inputFile, const std::string& config)
{
  MACH3LOG_INFO("File for study: {}", inputFile);

  YAML::Node Settings = M3OpenConfig(config);

  // Make the processor
  auto Processor = std::make_unique<MCMCProcessor>(inputFile);
  Processor->SetOutputSuffix("_Smear_MCMC");
  Processor->Initialise();

  const auto& Smear = Settings["SmearChain"];
  std::vector<std::string> Names = Smear["Smear"][0].as<std::vector<std::string>>();
  std::vector<double> ErrorValue = Smear["Smear"][1].as<std::vector<double>>();

  bool SaveUnsmearedBranch = GetFromManager<bool>(Smear["SaveUnsmearedBranch"], false);
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
  DiagMCMC(filename, config);

  return 0;
}
