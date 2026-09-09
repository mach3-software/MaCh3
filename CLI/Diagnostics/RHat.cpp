// MaCh3 includes
#include "Fitters/Processing/RHatCalculator.h"

/// @file RHat.cpp
/// @brief This executable calculates the \f$ \hat{R} \f$ estimator for Markov Chain Monte Carlo (MCMC) convergence.
///
/// KS: This exe is meant to calculate the \f$ \hat{R} \f$ estimator. For a well-converged chain, this distribution
/// should be centered at one. The \f$ \hat{R} \f$ statistic is used to assess the convergence of MCMC simulations
/// and helps determine whether the chains have reached a stable distribution.
///
/// @cite gelman2019.
///
/// MJR: Update -- Improved memory usage so that whole chains can be quickly loaded without requiring copious amounts
///      of RAM. This comes at the cost of not being able to calculate the Folded RHat since finding the median
///      requires the loading of full chains at a time. The method has been validated to give identical results to the
///      "High Memory" (original) version at a fraction of the runtime and resources.
///
///      The input format is also slightly altered; since we can now load entire chains, there's less need to
///      specify how many toys are desired for a sub-sample, so the Ntoys input has been removed.
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
  int Nchains = 0;
  if (argc < 2)
  {
    MACH3LOG_ERROR("Wrong arguments");
    MACH3LOG_ERROR("./RHat NThin MCMCchain_1.root MCMCchain_2.root MCMCchain_3.root ... [how many you like]");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  int NThin = atoi(argv[1]);

  //KS Gelman suggests to diagnose on more than one chain
  for (int i = 2; i < argc; i++)
  {
    MCMCFile.push_back(std::string(argv[i]));
    MACH3LOG_INFO("Adding file: {}", MCMCFile.back());
    Nchains++;
  }

  if(Nchains == 1)
  {
    MACH3LOG_WARN("Gelman is going to be sad :(. He suggested you should use more than one chain (at least 4). Code works fine for one chain, however, estimator might be biased.");
    MACH3LOG_WARN("Multiple chains are more likely to reveal multimodality and poor adaptation or mixing:");
  }
  MACH3LOG_INFO("Diagnosing {} chains", Nchains);
  auto RHatCalc = std::make_unique<RHatCalculator>(false, MCMCFile, NThin);
  RHatCalc->RunDiagnostic();

  return 0;
}
