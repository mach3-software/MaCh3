#pragma once

#include "Fitters/FitterBase.h"
#include "NuDockFactory.h"
#include "Samples/SampleHandlerFD.h"
_MaCh3_Safe_Include_Start_ //{
#include <nlohmann/json.hpp>
_MaCh3_Safe_Include_End_ //}

class NuDockServer : public FitterBase {
public:
  NuDockServer(manager * const fitMan);
  virtual ~NuDockServer();

  inline std::string GetName()const {return "NuDockServer";};

  void setup();

  virtual double getLogLikelihood();

  virtual nlohmann::json setParameters(const nlohmann::json &request);
  virtual nlohmann::json getParametersNames(const nlohmann::json &request);
  virtual nlohmann::json getParameters(const nlohmann::json &request);
  virtual nlohmann::json getLogLikelihood(const nlohmann::json &request);
  virtual nlohmann::json setAsimovPoint(const nlohmann::json &request);
  virtual nlohmann::json getSpectrum(const nlohmann::json &request);

  void RunMCMC() override {(void)0; /* do nothing */};

private:
  bool verbose;
  bool add_prior_llh;
};
