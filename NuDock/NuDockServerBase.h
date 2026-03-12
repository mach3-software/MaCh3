#pragma once

/// @file NuDockServerBase.h
/// @brief Base class for the MaCh3 NuDock server, exposing fit internals via a
///        JSON-based request/response protocol.
/// @author Hank Hua

#include "Fitters/FitterBase.h"
#include "NuDockFactory.h"
#include "Samples/SampleHandlerFD.h"
_MaCh3_Safe_Include_Start_ //{
#include <nlohmann/json.hpp>
_MaCh3_Safe_Include_End_ //}

/// @brief A class that implements a NuDock server for cross-fitter communication.
///
/// @details NuDockServerBase acts as the server-side adapter in the NuDock protocol.
/// It inherits from FitterBase and translates incoming JSON requests into
/// operations on MaCh3 systematic & sample objects (set/get parameters,
/// compute likelihoods, set Asimov data, retrieve predicted spectra, etc.).
///
/// Experiment-specific servers should derive from this class and override the
/// virtual handlers where custom behaviour is needed.
///
/// @author Hank Hua
class NuDockServerBase : public FitterBase {
public:
  /// @brief Construct the NuDock server.
  /// @param fitMan Pointer to the MaCh3 manager that owns fit configuration.
  NuDockServerBase(Manager* const fitMan);

  /// @brief Destructor.
  virtual ~NuDockServerBase();

  /// @brief Return a human-readable identifier for this fitter.
  inline std::string GetName()const {return "NuDockServer";};

  /// @brief Perform one-time server setup.
  ///
  /// Resizes internal likelihood storage vectors and reads the "Verbose" and
  /// "AddPriorLLH" flags from the NuDock configuration block.
  void setup();

  /// @brief Compute the total log-likelihood across all systematics and samples.
  ///
  /// The systematic prior likelihoods are optionally included (controlled by
  /// `add_prior_llh`). Sample likelihoods are always summed.
  ///
  /// @return The combined log-likelihood value.
  virtual double getLogLikelihood();

  /// @brief Handle a "set_parameters" request.
  ///
  /// Expects the JSON payload to contain "sys_pars" and "osc_pars" maps.
  /// Oscillation parameters are converted from NuDock convention to MaCh3
  /// convention before being applied. Samples are reweighted after all
  /// parameters are updated.
  ///
  /// @param request JSON object with "sys_pars" and "osc_pars" entries.
  /// @return JSON response with a "status" field.
  virtual nlohmann::json setParameters(const nlohmann::json &request);

  /// @brief Handle a "get_parameter_names" request.
  ///
  /// Returns the names of all unfixed, non-oscillation systematic parameters.
  ///
  /// @param request JSON object (unused).
  /// @return JSON response with a "sys_pars" array of parameter name strings.
  virtual nlohmann::json getParametersNames(const nlohmann::json &request);

  /// @brief Handle a "get_parameters" request.
  ///
  /// Returns current values of oscillation and systematic parameters.
  ///
  /// @param request JSON object (unused).
  /// @return JSON response with "osc_pars" and "sys_pars" maps.
  virtual nlohmann::json getParameters(const nlohmann::json &request);

  /// @brief Handle a "log_likelihood" request.
  ///
  /// Returns twice the log-likelihood (2NLL) to match NuDock convention.
  ///
  /// @param request JSON object (unused).
  /// @return JSON response with a "log_likelihood" field containing 2*NLL.
  virtual nlohmann::json get2NLL(const nlohmann::json &request);

  /// @brief Handle a "set_asimov" request.
  ///
  /// Reweights all samples then clones the MC prediction as the "data".
  /// Prefit parameter central values are updated to match the current proposed values.
  ///
  /// @param request JSON object (unused).
  /// @return JSON response with a "status" field.
  virtual nlohmann::json setAsimovPoint(const nlohmann::json &request);

  /// @brief Handle a "get_spectrum" request.
  ///
  /// Retrieves predicted bin contents and bin edges from every registered
  /// SampleHandlerFD sample, supporting 1-D and 2-D histograms.
  ///
  /// @param request JSON object (unused).
  /// @return JSON response containing sample names, dimensions, axis titles,
  ///         bin edges, and bin values.
  /// @throw MaCh3Exception if a sample does not derive from SampleHandlerFD.
  virtual nlohmann::json getSpectrum(const nlohmann::json &request);

  /// @brief No-op implementation --- the server does not run own MCMC.
  void RunMCMC() override {(void)0; /* do nothing */};

protected:
  /// @brief Flag controlling verbose logging of server operations.
  bool verbose;
  /// @brief If true, systematic prior likelihoods are added to the total LLH.
  bool add_prior_llh;
};
