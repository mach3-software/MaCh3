/// @file UmbrellaPlotting.cpp
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author David Riley

#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/UmbrellaPlotModule.hpp"

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 umbrellasolver plot' as a direct replacement instead.");
  argv[0] = "plot";
  M3::UmbrellaPlotModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();
  parser->parse_args(argc, argv);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 umbrellasolver plot' as a direct replacement instead.");
  return 0;
}
