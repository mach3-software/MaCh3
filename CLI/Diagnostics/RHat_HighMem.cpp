// MaCh3 includes
#include "Fitters/RHatCalculator.h"

/// @file RHat_HighMem.cpp
/// @brief This executable calculates the \f$ \hat{R} \f$ estimator for Markov Chain Monte Carlo (MCMC) convergence.
///
/// KS: This exe is meant to calculate the \f$ \hat{R} \f$ estimator. For a well-converged chain, this distribution
/// should be centered at one. The \f$ \hat{R} \f$ statistic is used to assess the convergence of MCMC simulations
/// and helps determine whether the chains have reached a stable distribution.
///
/// @cite gelman2019.
///
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author Kamil Skwarczynski
/// @author Michael Reh

// *******************

// *******************
int main(int argc, char *argv[]) {
// *******************
  SetMaCh3LoggerFormat();
  M3::Utils::MaCh3Welcome();

  std::vector<std::string> MCMCFile;
  if (argc == 1 || argc == 2)
  {
    MACH3LOG_ERROR("Wrong arguments");
    MACH3LOG_ERROR("./RHat Ntoys MCMCchain_1.root MCMCchain_2.root MCMCchain_3.root ... [how many you like]");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  int Ntoys = atoi(argv[1]);
  int Nchains = 0;
  //KS Gelman suggests to diagnose on more than one chain
  for (int i = 2; i < argc; i++)
  {
    MCMCFile.push_back(std::string(argv[i]));
    MACH3LOG_INFO("Adding file: {}", MCMCFile.back());
    Nchains++;
  }

  if(Ntoys < 1)
  {
    MACH3LOG_ERROR("You specified {} specify larger greater than 0", Ntoys);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if(Nchains == 1)
  {
    MACH3LOG_WARN("Gelman is going to be sad :(. He suggested you should use more than one chain (at least 4). Code works fine for one chain, however, estimator might be biased.");
    MACH3LOG_WARN("Multiple chains are more likely to reveal multimodality and poor adaptation or mixing:");
  }
  MACH3LOG_INFO("Diagnosing {} chains, with {} toys", Nchains, Ntoys);

  auto RHatCalc = std::make_unique<RHatCalculator>(true, MCMCFile, Ntoys);
  RHatCalc->RunDiagnostic();
  return 0;
}
