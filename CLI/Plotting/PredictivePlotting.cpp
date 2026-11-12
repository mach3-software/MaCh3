/// @file PredictivePlotting.cpp
/// @ingroup MaCh3Plotting
/// @author Kamil Skwarczynski

/// @warning KS: keep raw pointer or ensure manual delete of PlotMan. If spdlog in automatically deleted before PlotMan then destructor has some spdlog and this could cause segfault
#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/Plot/PredictiveModule.hpp"

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 plot predictive' as a direct replacement instead.");
  argv[0] = "predictive";
  M3::PredictiveModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();
  parser->parse_args(argc, argv);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 plot predictive' as a direct replacement instead.");
  return 0;
}
