#include "NuDock/SampleHandlerNuDockBase.h"
#include <unordered_map>
#include <algorithm>
#include <NuDock/NuDockFactory.h>

// ***************************************************************************
/// @copydoc SampleHandlerNuDockBase::SampleHandlerNuDockBase
// ***************************************************************************
SampleHandlerNuDockBase::SampleHandlerNuDockBase(std::string configFile, ParameterHandlerGeneric* xsec_cov)
: SampleHandlerInterface() {
  MACH3LOG_INFO("Creating SampleHandlerNuDock object..");
  MACH3LOG_INFO("- Using NuDock sample config in this file {}", configFile);
  ParHandler = xsec_cov;
  SampleManager = std::make_unique<Manager>(configFile.c_str());
  verbose = GetFromManager(SampleManager->raw()["NuDockClient"]["Verbose"], false, __FILE__, __LINE__);
  checkServerOscParams = GetFromManager(SampleManager->raw()["NuDockClient"]["CheckServerOscParams"], true, __FILE__, __LINE__);

  // ReadConfig();
  // SetupReweightArrays();
  Init();
}

// ***************************************************************************
/// @copydoc SampleHandlerNuDockBase::~SampleHandlerNuDockBase
// ***************************************************************************
SampleHandlerNuDockBase::~SampleHandlerNuDockBase() {
}

// ***************************************************************************
/// @copydoc SampleHandlerNuDockBase::Init
// ***************************************************************************
void SampleHandlerNuDockBase::Init() {
  InitialiseNuDockObj(SampleManager.get(), nudock_ptr);

  // Gather the indices of NuDock parameters from the ParameterHandler, so that we can easily retrieve their values in Reweight()
  auto nudockParamInds_func = ParHandler->GetParsIndexFromSampleName(kNuDockSampleTag, SystType::kFunc);
  // HH: We should only be using kFunc for NuDock, but for convenience of testing (e.g. when validating
  // M3 client with M3 server we'd want to test the LLH of all the params) we include all param types here.
  auto nudockParamInds_norm = ParHandler->GetParsIndexFromSampleName(kNuDockSampleTag, SystType::kNorm);
  auto nudockParamInds_spline = ParHandler->GetParsIndexFromSampleName(kNuDockSampleTag, SystType::kSpline);
  nudockParamInds = nudockParamInds_func;
  nudockParamInds.insert(nudockParamInds.end(), nudockParamInds_norm.begin(), nudockParamInds_norm.end());
  nudockParamInds.insert(nudockParamInds.end(), nudockParamInds_spline.begin(), nudockParamInds_spline.end());

  // Oscillation parameters now should be tagged in the covariance YAML
  AssertNuDockOscParamsTagged(ParHandler);
  nudockOscParamInds = GetNuDockOscParIndices(ParHandler);
  MACH3LOG_INFO("NuDock oscillation-parameter set contains {} parameter(s)", nudockOscParamInds.size());
  for (const auto &iParam : nudockOscParamInds) {
    MACH3LOG_INFO("  - {}", ParHandler->GetParFancyName(iParam));
  }

  nudock_ptr->start_client();

  if (checkServerOscParams) {
    CheckServerOscParams();
  } else {
    MACH3LOG_WARN("NuDockClient:CheckServerOscParams is off. A client/server disagreement on the");
    MACH3LOG_WARN("oscillation parameter set will not be detected, and any parameter the server");
    MACH3LOG_WARN("does not recognise will silently keep its previous value.");
  }
}

// ***************************************************************************
/// @copydoc SampleHandlerNuDockBase::CheckServerOscParams
// ***************************************************************************
void SampleHandlerNuDockBase::CheckServerOscParams() {
  nlohmann::json request = "";
  nlohmann::json response;
  try {
    response = nudock_ptr->send_request("/get_parameter_names", request);
  } catch (const std::exception &e) {
    MACH3LOG_WARN("Could not query server parameter names ({}); skipping oscillation", e.what());
    MACH3LOG_WARN("parameter check. A client/server mismatch will not be detected.");
    return;
  }

  if (!response.contains("osc_pars") || !response["osc_pars"].is_array() || response["osc_pars"].empty()) {
    MACH3LOG_WARN("Server reported no osc_pars list; skipping oscillation parameter");
    MACH3LOG_WARN("check. A client/server mismatch will not be detected.");
    return;
  }

  const auto serverNames = response["osc_pars"].get<std::vector<std::string>>();

  std::vector<std::string> missingOnServer;
  for (const auto &iParam : nudockOscParamInds) {
    const std::string nameM3 = ParHandler->GetParFancyName(iParam);
    const auto it = NuDockOscNameMap_r.find(nameM3);
    const std::string nameNuDock = (it != NuDockOscNameMap_r.end()) ? it->second : nameM3;
    if (std::find(serverNames.begin(), serverNames.end(), nameNuDock) == serverNames.end()) {
      missingOnServer.push_back(nameNuDock);
    }
  }

  if (!missingOnServer.empty()) {
    MACH3LOG_ERROR("Client and server disagree on the oscillation parameter set.");
    MACH3LOG_ERROR("The server does not recognise the following, so it would silently keep them");
    MACH3LOG_ERROR("at their previous values rather than applying what we send:");
    for (const auto &name : missingOnServer) MACH3LOG_ERROR("  - {}", name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  MACH3LOG_INFO("Oscillation parameter set checked against server ({} parameters)", serverNames.size());
}

// ***************************************************************************
/// @copydoc SampleHandlerNuDockBase::Reweight
// ***************************************************************************
void SampleHandlerNuDockBase::Reweight() {
  nlohmann::json request;
  std::unordered_map<std::string, double> osc_params;
  std::unordered_map<std::string, double> xsec_params;
  
  // Loop for systs
  for (const auto& iParam : nudockParamInds) {
    std::string paramName = ParHandler->GetParFancyName(iParam);
    double paramValue = ParHandler->GetParProp(iParam);
    xsec_params[paramName] = paramValue;
  }

  // Loop over oscillation params tagged with nudock sample
  for (const auto& iParam : nudockOscParamInds) {
    const std::string paramNameM3 = ParHandler->GetParFancyName(iParam);
    double paramValue = ParHandler->GetParProp(iParam);
    // Convert sin2_theta to theta
    FormatOscParsForNuDock(paramNameM3, paramValue);
    // Translate MaCh3's name to NuDock's name
    const auto it = NuDockOscNameMap_r.find(paramNameM3);
    osc_params[(it != NuDockOscNameMap_r.end()) ? it->second : paramNameM3] = paramValue;
  }

  request["osc_pars"] = osc_params;
  request["sys_pars"] = xsec_params;

  auto response = nudock_ptr->send_request("/set_parameters", request);
  if (verbose) {
    try {
      MACH3LOG_INFO("NuDock response: {}", response.dump());
    } catch (const std::exception &e) {
      MACH3LOG_ERROR("Error dumping NuDock response: {}", e.what());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
}

// ***************************************************************************
/// @copydoc SampleHandlerNuDockBase::GetLikelihood
// ***************************************************************************
double SampleHandlerNuDockBase::GetLikelihood() const {
  nlohmann::json request = "";
  double llh_value = 0.0;
  auto response = nudock_ptr->send_request("/log_likelihood", request);
  try {
    llh_value = response["log_likelihood"].get<double>();
    llh_value /= 2; // NuDock returns 2NLL, so we divide by 2 to be consistent with M3's definition of LLH.
    return llh_value;
  } catch (const std::exception &e) {
    MACH3LOG_ERROR("Error retrieving log-likelihood from NuDock response: {}", e.what());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}
