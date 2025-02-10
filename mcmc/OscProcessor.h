#pragma once


// MaCh3 includes
#include "mcmc/MCMCProcessor.h"
#include "samplePDF/HistogramUtils.h"

/// @author Clarence Wret
/// @author Kamil Skwarczynski
class OscProcessor : public MCMCProcessor {
  public:
    /// @brief Constructs an OscProcessor object with the specified input file and options.
    /// @param InputFile The path to the input file containing MCMC data.
    OscProcessor(const std::string &InputFile);
    /// @brief Destroys the OscProcessor object.
    virtual ~OscProcessor();
};
