#include "NuDock/SampleHandlerNuDock.h"
#include <unordered_map>
#include <NuDock/NuDockFactory.h>
#include "NuDock/NuDockServer.h"

SampleHandlerNuDock::SampleHandlerNuDock(std::string configFile, ParameterHandlerGeneric* xsec_cov, const std::shared_ptr<OscillationHandler>& Oscillator)
: SampleHandlerFD(configFile, xsec_cov, Oscillator) {
  MACH3LOG_INFO("Creating SampleHandlerNuDock object..");
  MACH3LOG_INFO("- Using NuDockT2K sample config in this file {}", configFile);
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
  nudock_ptr->start_client();
}

void SampleHandlerNuDock::Reweight() {
  nlohmann::json request;
  std::unordered_map<std::string, double> osc_params;
  std::unordered_map<std::string, double> xsec_params;
  
  // Loop for systs
  auto nudockParamInds_func = ParHandler->GetParsIndexFromSampleName("NuDock", SystType::kFunc);
  // HH: We should only be using kFunc for NuDock, but for convenience of testing (e.g. when validation
  // M3 client with M3 server we'd want to test the LLH of all the params) we include all param types here.
  auto nudockParamInds_norm = ParHandler->GetParsIndexFromSampleName("NuDock", SystType::kNorm);
  auto nudockParamInds_spline = ParHandler->GetParsIndexFromSampleName("NuDock", SystType::kSpline);
  std::vector<int> nudockParamInds = nudockParamInds_func;
  nudockParamInds.insert(nudockParamInds.end(), nudockParamInds_norm.begin(), nudockParamInds_norm.end());
  nudockParamInds.insert(nudockParamInds.end(), nudockParamInds_spline.begin(), nudockParamInds_spline.end());
  for (const auto& iParam : nudockParamInds) {
    std::string paramName = ParHandler->GetParFancyName(iParam);
    double paramValue = ParHandler->GetParProp(iParam);
    xsec_params[paramName] = paramValue;
  }

  // Loop for osc
  // HH: will use SK so we don't need to change the osc config
  // TODO: this needs to be changed to be more general
  auto oscParamInds = ParHandler->GetParsIndexFromSampleName("SK", SystType::kOsc);
  for (const auto& iParam : oscParamInds) {
    std::string paramName = ParHandler->GetParFancyName(iParam);
    double paramValue = ParHandler->GetParProp(iParam);
    // Skip if not any of the osc params NuDock uses
    if (NuDockOscNameMap_r.find(paramName) == NuDockOscNameMap_r.end()) continue;
    paramName = NuDockOscNameMap_r.at(paramName);
    // Convert sin2_theta to theta
    if (paramName == "Theta12" || paramName == "Theta13" || paramName == "Theta23") {
      paramValue = asin(sqrt(paramValue));
    }
    osc_params[paramName] = paramValue;
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
  // if (response.contains("status") && response["status"] == "success") {
  //   MACH3LOG_INFO("Successfully set parameters in NuDock.");
  // } else {
  //   MACH3LOG_ERROR("Failed to set parameters in NuDock.");
  //   throw MaCh3Exception(__FILE__, __LINE__);
  // }
}

double SampleHandlerNuDock::GetLikelihood() const {
  MACH3LOG_INFO("Requesting log-likelihood from NuDock.");
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
