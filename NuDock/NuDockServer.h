#pragma once

#include "Fitters/FitterBase.h"
#include "NuDockFactory.h"
_MaCh3_Safe_Include_Start_ //{
#include <nlohmann/json.hpp>
_MaCh3_Safe_Include_End_ //}

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
  nlohmann::json getLogLikelihood(const nlohmann::json &request);
  nlohmann::json setAsimovPoint(const nlohmann::json &request);

private:
  bool verbose;
  bool add_prior_llh;
};
