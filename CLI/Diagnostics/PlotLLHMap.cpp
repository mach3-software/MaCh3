/// @file PlotLLHMap.cpp
/// @brief Processing n-dimensional LLHMap outputs generating 1D and 2D profiled likelihoods as defined in config
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author Tomas Nosek

#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/LLHMapModule.hpp"

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 diagnostics llhmap' as a direct replacement instead.");
  argv[0] = "llhmap";
  M3::LLHMapModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();
  parser->parse_args(argc, argv);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 diagnostics llhmap' as a direct replacement instead.");
  return 0;
}
