/// @file RHatModule.cpp
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
#include "Fitters/Processing/RHatCalculator.h"
#include "CLI/Modules/RHatModule.hpp"

namespace M3{

  RHatModule::~RHatModule() = default;

  MaCh3ArgumentParser* RHatModule::get_parser(){
    m_parser = std::make_unique<MaCh3ArgumentParser>("rhat", "1.0", argparse::default_arguments::help);
    m_parser->add_description("Tool for computing the R-hat diagnostic for MCMC chains.");
    m_parser->add_argument("n-thin")
      .help("if in high mem mode this is number of toys per chain, for low mem thinning setting")
      .metavar("NTHINTOYS")
      .scan<'i', int>()
      .required();
    m_parser->add_argument("--high-mem")
      .help("Enable high memory mode.")
      .flag();
    m_parser->add_argument("mcmc-chain")
      .help("MCMC chain root files")
      .metavar("MCMC_CHAIN1 [MCMC_CHAIN2 [...]]")
      .nargs(argparse::nargs_pattern::at_least_one)
      .required();
    return m_parser.get();
  }

  int RHatModule::Run() {
    SetMaCh3LoggerFormat();
    M3::Utils::MaCh3Welcome();

    std::vector<std::string> MCMCFile;
    auto mcmc_chain_args = m_parser->get<std::vector<std::string>>("mcmc-chain");

    if (mcmc_chain_args.empty()){  // No MCMC chain files provided, shouldn't be possible with the argparser req.
        MACH3LOG_ERROR("Wrong arguments");
        MACH3LOG_ERROR("./RHat NThin MCMCchain_1.root MCMCchain_2.root MCMCchain_3.root ... [how many you like]");
        throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    int NThinToys = m_parser->get<int>("n-thin");  // Number of thinning for MCMC chains, Number of toys for high memory mode
    std::size_t Nchains = mcmc_chain_args.size();

    //KS Gelman suggests to diagnose on more than one chain
    for (const auto& file : mcmc_chain_args){
        MACH3LOG_INFO("Adding file: {}", file);
    }

    if(Nchains == 1){
        MACH3LOG_WARN("Gelman is going to be sad :(. He suggested you should use more than one chain (at least 4). Code works fine for one chain, however, estimator might be biased.");
        MACH3LOG_WARN("Multiple chains are more likely to reveal multimodality and poor adaptation or mixing:");
    }


    MACH3LOG_INFO("Diagnosing {} chains", Nchains);
    bool high_mem = m_parser->get<bool>("--high-mem");
    if (high_mem){
        MACH3LOG_INFO("With {} toys", NThinToys);
        if (NThinToys < 1){
            MACH3LOG_ERROR("You specified {} for NToys, specify larger greater than 0", NThinToys);
            throw MaCh3Exception(__FILE__ , __LINE__ );
        }
    }

    auto RHatCalc = std::make_unique<RHatCalculator>(high_mem, mcmc_chain_args, NThinToys);
    RHatCalc->RunDiagnostic();

    return 0;
  }

}