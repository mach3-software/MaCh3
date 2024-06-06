#include "mcmc/MCMCProcessor.h"

#include "manager/manager.h"


void DiagMCMC(std::string inputFile, std::string config);

int main(int argc, char *argv[]) {

  if (argc != 3)
  {
    std::cerr << "How to use:   DiagMCMC MCMC_Output.root config" << std::endl;
    exit(-1);
  }
  
  std::cout << "Producing single fit output" << std::endl;
  std::string filename = argv[1];
  std::string config = argv[2];
  DiagMCMC(filename, config);
    
  return 0;
}


void DiagMCMC(std::string inputFile, std::string config)
{
    std::cout << "File for study:       " << inputFile << std::endl;
      
    YAML::Node Settings = YAML::LoadFile(config);

   // Make the processor
    MCMCProcessor* Processor = new MCMCProcessor(inputFile, false);
    
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


