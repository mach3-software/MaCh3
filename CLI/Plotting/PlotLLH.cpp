// C++
#include <algorithm>
#include <iomanip>

// MACH3 PLOTTING
#include "Utils/Plotting/PlottingManager.h"

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file PlotLLH.cpp
/// @ingroup MaCh3Plotting
/// @author Ewan Miller

#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/Plot/LLHModule.hpp"

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 plot llhscan' as a direct replacement instead.");
  argv[0] = "llhscan";
  M3::LLHModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();
  parser->parse_args(argc, argv);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 plot llhscan' as a direct replacement instead.");
  return 0;
}
