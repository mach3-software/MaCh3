#include "mcmc/MCMCProcessor.h"

#include "manager/manager.h"


void DiagMCMC(std::string inputFile, std::string config);

int main(int argc, char *argv[]) {
  SetMaCh3LoggerFormat();
  if (argc != 3)
  {
    MACH3LOG_ERROR("How to use: DiagMCMC MCMC_Output.root config");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  MACH3LOG_INFO("Producing single fit output");
  std::string filename = argv[1];
  std::string config = argv[2];
  DiagMCMC(filename, config);
    
  return 0;
}


void DiagMCMC(std::string inputFile, std::string config)
{
  MACH3LOG_INFO("File for study: {}", inputFile);

  YAML::Node Settings = YAML::LoadFile(config);

  // Make the processor
  MCMCProcessor* Processor = new MCMCProcessor(inputFile);

  Processor->SetOutputSuffix("_MCMC_Diag");
  //KS:Turn off plotting detector and some other setting
  Processor->SetExcludedTypes(GetFromManager<std::vector<std::string>>(Settings["DiagMCMC"]["ExcludedTypes"], {""}));
  Processor->SetExcludedNames(GetFromManager<std::vector<std::string>>(Settings["DiagMCMC"]["ExcludedNames"], {""}));
  Processor->SetPlotRelativeToPrior(GetFromManager<bool>(Settings["DiagMCMC"]["PlotRelativeToPrior"], false));
  //KS: Use 20 batches for batched means
  Processor->SetnBatches(GetFromManager<int>(Settings["DiagMCMC"]["nBatches"], 20));
  Processor->SetnLags(GetFromManager<int>(Settings["DiagMCMC"]["nLags"], 25000));
  Processor->SetPrintToPDF(GetFromManager<bool>(Settings["PrintToPDF"], true));
  Processor->Initialise();

  //KS: finally call main method
  Processor->DiagMCMC();

  delete Processor;
}


