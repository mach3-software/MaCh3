#pragma once

// MaCh3 includes
#include "Fitters/MCMCProcessor.h"
#include "Samples/HistogramUtils.h"

/// @author Clarence Wret
/// @author Kamil Skwarczynski
/// @brief This class extends MCMC and allow specialised for Oscillation parameters analysis which require specialised hardcoding
class OscProcessor : public MCMCProcessor {
  public:
    /// @brief Constructs an OscProcessor object with the specified input file and options.
    /// @param InputFile The path to the input file containing MCMC data.
    OscProcessor(const std::string &InputFile);
    /// @brief Destroys the OscProcessor object.
    virtual ~OscProcessor();

  protected:
    /// @brief Read the Osc cov file and get the input central values and errors
    /// Here we allow Jarlskog Shenanigans
    virtual void ReadOSCFile() override;

    /// Will plot Jarlskog Invariant using information in the chain
    bool PlotJarlskog;
};
