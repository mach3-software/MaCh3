/// @file UmbrellaSolver.cpp
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author David Riley

#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/UmbrellaModule.hpp"

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 umbrellasolver' as a direct replacement instead.");
  argv[0] = "umbrellasolver";
  M3::UmbrellaModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();
  std::vector<std::string> args(argv, argv + argc);
  parser->insert_default_subcommands(args);  // The UmbrellaModule parser doesn't automatically insert default subcommand so treats the first argument as a subcommand
  parser->parse_args(args);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 umbrellasolver' as a direct replacement instead.");
  return 0;
}

