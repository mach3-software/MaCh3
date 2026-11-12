// MaCh3 includes
#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/Plot/DiagMCMCPlotModule.hpp"

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 plot diagmcmc' as a direct replacement instead.");
  argv[0] = "diagmcmc";
  M3::DiagMCMCPlotModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();
  parser->parse_args(argc, argv);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 plot diagmcmc' as a direct replacement instead.");
  return 0;
}
