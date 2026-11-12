#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/DiagMCMCModule.hpp"

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 diagnostics mcmc' as a direct replacement instead.");
  argv[0] = "mcmc";
  M3::DiagMCMCModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();

  // add compatibility layer for old binaries which have config file after input
  // ////////////////////////////////////////////////////////////////////////////
  if(argc > 2){
    std::string arg1 = argv[1];
    std::string arg2 = argv[2];
    if(arg1.size() > 5 && arg1.substr(arg1.size() - 5) == ".root" &&
      ((arg2.size() > 5 && arg2.substr(arg2.size() - 5) == ".yaml") ||
      (arg2.size() > 4 && arg2.substr(arg2.size() - 4) == ".yml"))){
      std::swap(argv[1], argv[2]);
    }
  }
  // ////////////////////////////////////////////////////////////////////////////

  std::vector<std::string> args(argv, argv + argc);
  parser->insert_default_subcommands(args);  // The DiagMCMCModule parser doesn't automatically insert default subcommand so treats the first argument as a subcommand
  parser->parse_args(args);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 diagnostics mcmc' as a direct replacement instead.");
  return 0;
}
