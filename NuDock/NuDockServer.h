#pragma once

#include "Fitters/FitterBase.h"
#include "Samples/SampleHandlerFD.h"
#include <nlohmann/json.hpp>

class NuDockServer : public FitterBase {
public:
  NuDockServer(manager * const fitMan);
  virtual ~NuDockServer();

  inline std::string GetName()const {return "NuDockServer";};

  void setup();

  double getLogLikelihood();

  void RunMCMC() override {(void)0; /* do nothing */};

  nlohmann::json setParameters(const nlohmann::json &request);
  nlohmann::json getParametersNames(const nlohmann::json &request);
  nlohmann::json getParameters(const nlohmann::json &request);
  nlohmann::json unnormaliseParameters(const nlohmann::json &request);
  nlohmann::json getLogLikelihood(const nlohmann::json &request);
  nlohmann::json setAsimovPoint(const nlohmann::json &request);

  std::unordered_map<std::string, std::string> NuDockOscNameMap = {
    {"Theta12", "sin2th_12"},
    {"Theta13", "sin2th_13"},
    {"Theta23", "sin2th_23"},
    {"DeltaCP", "delta_cp"},
    {"Deltam2_21", "delm2_12"},
    {"Deltam2_32", "delm2_23"},
  };

  std::unordered_map<std::string, std::string> NuDockOscNameMap_r = {
    {"sin2th_12", "Theta12"},
    {"sin2th_13", "Theta13"},
    {"sin2th_23", "Theta23"},
    {"delta_cp", "DeltaCP"},
    {"delm2_12", "Deltam2_21"},
    {"delm2_23", "Deltam2_32"},
  };

private:
  bool verbose;
  bool add_prior_llh;
};
