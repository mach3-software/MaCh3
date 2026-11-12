#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/CombineModule.hpp"

/// @file CombineMaCh3Chains.cpp
/// @brief Combine chains files produced by **MCMC**, enforcing the condition that all the files to combine were made using the exact same software versions and config files.
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author Ewan Miller
/// @author Kamil Skwarczynski

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 combine' as a direct replacement instead.");
  argv[0] = "combine";
  M3::CombineModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();
  parser->parse_args(argc, argv);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 combine' as a direct replacement instead.");
  return 0;
}
