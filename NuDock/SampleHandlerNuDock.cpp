#include "NuDock/SampleHandlerNuDock.h"
#include <unordered_map>
#include <NuDock/NuDockFactory.h>

SampleHandlerNuDock::SampleHandlerNuDock(std::string configFile, ParameterHandlerGeneric* xsec_cov, const std::shared_ptr<OscillationHandler>& Oscillator)
: SampleHandlerFD(configFile, xsec_cov, Oscillator) {
  MACH3LOG_INFO("Creating SampleHandlerNuDock object..");
  MACH3LOG_INFO("- Using NuDock sample config in this file {}", configFile);
  ParHandler = xsec_cov;
  SampleManager = std::make_unique<manager>(configFile.c_str());
  verbose = GetFromManager(SampleManager->raw()["NuDockClient"]["Verbose"], false);

  ReadConfig();
  SetupReweightArrays();
  Init();
}

SampleHandlerNuDock::~SampleHandlerNuDock() {
}

void SampleHandlerNuDock::Init() {
  InitialiseNuDockObj(SampleManager.get(), nudock_ptr);

  // Gather the indices of NuDock parameters from the ParameterHandler, so that we can easily retrieve their values in Reweight()
  auto nudockParamInds_func = ParHandler->GetParsIndexFromSampleName("NuDock", SystType::kFunc);
  // HH: We should only be using kFunc for NuDock, but for convenience of testing (e.g. when validating
  // M3 client with M3 server we'd want to test the LLH of all the params) we include all param types here.
  auto nudockParamInds_norm = ParHandler->GetParsIndexFromSampleName("NuDock", SystType::kNorm);
  auto nudockParamInds_spline = ParHandler->GetParsIndexFromSampleName("NuDock", SystType::kSpline);
  nudockParamInds = nudockParamInds_func;
  nudockParamInds.insert(nudockParamInds.end(), nudockParamInds_norm.begin(), nudockParamInds_norm.end());
  nudockParamInds.insert(nudockParamInds.end(), nudockParamInds_spline.begin(), nudockParamInds_spline.end());

  nudock_ptr->start_client();
}

// HH: Instead of reweighting here, we send the parameters to the NuDock server
void SampleHandlerNuDock::Reweight() {
  nlohmann::json request;
  std::unordered_map<std::string, double> osc_params;
  std::unordered_map<std::string, double> xsec_params;
  
  // Loop for systs
  for (const auto& iParam : nudockParamInds) {
    std::string paramName = ParHandler->GetParFancyName(iParam);
    double paramValue = ParHandler->GetParProp(iParam);
    xsec_params[paramName] = paramValue;
  }

  // Loop over NuDockOscNameMap_r to get osc params
  for (auto const& [paramNameM3, paramNameNuDock] : NuDockOscNameMap_r) {
    int iParam = ParHandler->GetParIndex(paramNameM3);
    double paramValue = ParHandler->GetParProp(iParam);
    // Convert sin2_theta to theta
    FormatOscParsForNuDock(paramNameNuDock, paramValue);
    osc_params[paramNameNuDock] = paramValue;
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

// HH: Instead of calculating the likelihood from M3, we get it from NuDock server
double SampleHandlerNuDock::GetLikelihood() const {
  nlohmann::json request = "";
  double llh_value = 0.0;
  auto response = nudock_ptr->send_request("/log_likelihood", request);
  try {
    llh_value = response["log_likelihood"].get<double>();
    return llh_value;
  } catch (const std::exception &e) {
    MACH3LOG_ERROR("Error retrieving log-likelihood from NuDock response: {}", e.what());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}
