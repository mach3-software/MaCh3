/// @file GetPostfitParamPlots.cpp
/// This script generates post-fit parameter plots. The central postfit value is
/// taken as the Highest Posterior Density (HPD), but can be easily changed to
/// another method such as Gaussian. Be cautious as parameter names and the number
/// of parameters per plot are currently hardcoded.
///
/// @details
/// Usage:
/// ```
/// ./GetPostfitParamPlots ProcessMCMC_Output1.root <ProcessMCMC_Output2.root> <ProcessMCMC_Output3.root>
/// ```
///
/// @ingroup MaCh3Plotting
///
/// @author Clarence Wret
/// @author Will Parker
/// @author Kamil Skwarczynski
/// @author Ewan Miller


#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/Plot/PostFitParamModule.hpp"

int main(int argc, char const* argv[]){
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 plot postfitparam' as a direct replacement instead.");
  argv[0] = "postfitparam";
  M3::PostFitParamModule proc;
  M3::MaCh3ArgumentParser* parser = proc.get_parser();
  parser->parse_args(argc, argv);
  proc.Run();
  MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
  MACH3LOG_WARN("                   : you can use 'mach3 plot postfitparam' as a direct replacement instead.");
  return 0;
}
