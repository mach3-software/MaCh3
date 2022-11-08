#include "mcmc/MCMCProcessor.h"

void DiagMCMC(std::string inputFile);

int main(int argc, char *argv[]) {

  if (argc != 2) 
  {
    std::cerr << "How to use:   DiagMCMC MCMC_Output.root" << std::endl;
    exit(-1);
  }
  
    if (argc == 2) 
    {
        std::cout << "Producing single fit output" << std::endl;
        std::string filename = argv[1];
        DiagMCMC(filename);
  }
    
  return 0;
}


void DiagMCMC(std::string inputFile)
{
    std::cout << "File for study:       " << inputFile << std::endl;
      
   // Make the processor
    MCMCProcessor* Processor = new MCMCProcessor(inputFile, false);
    
    Processor->SetOutputSuffix("_MCMC_Diag");
    //KS:Turn off plotting detector and some other setting, should be via some config
    Processor->SetPlotDet(false);
    Processor->SetPlotRelativeToPrior(true);
    Processor->Initialise();
    //KS: Use 20 batches for batched means
    Processor->SetnBatches(20);
    
    //KS: finally call main method
    Processor->DiagMCMC();
    
    delete Processor;
}


