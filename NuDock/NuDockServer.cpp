#include "NuDockServer.h"

NuDockServer::NuDockServer(manager *man) : FitterBase(man) {
}

NuDockServer::~NuDockServer() {
}

void NuDockServer::setup() {
  // Resize the likelihood storage vectors
  syst_llh.resize(systematics.size(), M3::_BAD_DOUBLE_);
  sample_llh.resize(samples.size(), M3::_BAD_DOUBLE_);
  verbose = GetFromManager(fitMan->raw()["NuDock"]["Verbose"], false);
  add_prior_llh = GetFromManager(fitMan->raw()["NuDock"]["AddPriorLLH"], false);
}

double NuDockServer::getLogLikelihood() {
  // HH: Based on MR2T2::ProposeStep()
  // Initial likelihood
  double llh = 0.0;

  // Initiate to false
  bool reject = false;

  if (verbose) {
    std::cout << "==Before reweighting==" << std::endl;
    for (size_t i = 0; i < samples.size(); ++i) {
      // Get the sample likelihoods and add them
      sample_llh[i] = samples[i]->GetLikelihood();
      std::cout << "Total LLH in sample handler " << i << ": " << sample_llh[i] << std::endl;
      for (int iSample = 0; iSample < samples[i]->GetNsamples(); ++iSample) {
        double sample_llh_ind = samples[i]->GetSampleLikelihood(iSample);
        std::cout << "  Sample " << samples[i]->GetSampleTitle(iSample)
                  << " llh: " << sample_llh_ind << std::endl;
      }
    }
    std::cout << "======================" << std::endl;
  }

  // HH: Only add prior llh if requested
  if (add_prior_llh) {
    // Loop over the systematics and propose the initial step
    for (size_t s = 0; s < systematics.size(); ++s) {
      // Get the likelihood from the systematics
      syst_llh[s] = systematics[s]->GetLikelihood();
      if (verbose) std::cout << "Syst " << systematics[s]->GetName() << " llh: " << syst_llh[s] << std::endl;
      llh += syst_llh[s];

      if (verbose) std::cout << "LLH after " << systematics[s]->GetName() << " " << llh << std::endl;
    }

    // Check if we've hit a boundary in the systematics
    // In this case we can save time by not having to reconfigure the simulation
    if (llh >= M3::_LARGE_LOGL_) {
      reject = true;
      if (verbose) std::cout << "Rejecting based on boundary" << std::endl;
    }
  }

  // Only reweight when we have a good parameter configuration
  // This speeds things up considerably because for every bad parameter configuration we don't have to reweight the MC
  if (!reject) {
    // Could multi-thread this
    // But since sample reweight is multi-threaded it's probably better to do that
    for (size_t i = 0; i < samples.size(); ++i) {
      samples[i]->Reweight();
    }

    //DB for atmospheric event by event sample migration, need to fully reweight all samples to allow event passing prior to likelihood evaluation
    for (size_t i = 0; i < samples.size(); ++i) {
      // Get the sample likelihoods and add them
      sample_llh[i] = samples[i]->GetLikelihood();
      if (verbose) {
        std::cout << "Total LLH in sample handler " << i << ": " << sample_llh[i] << std::endl;
        for (int iSample = 0; iSample < samples[i]->GetNsamples(); ++iSample) {
          double sample_llh_ind = samples[i]->GetSampleLikelihood(iSample);
          std::cout << "  Sample " << samples[i]->GetSampleTitle(iSample)
                    << " llh: " << sample_llh_ind << std::endl;
        }
      }
      llh += sample_llh[i];
      if (verbose) std::cout << "LLH after sample " << i << " " << llh << std::endl;
    }

    // For when we don't have to reweight, set sample to madness
  } else {
    for (size_t i = 0; i < samples.size(); ++i) {
      // Set the sample_llh[i] to be madly high also to signify a step out of bounds
      sample_llh[i] = M3::_LARGE_LOGL_;
      if (verbose) std::cout << "LLH after REJECT sample " << i << " " << llh << std::endl;
    }
  }

  return llh;

}

nlohmann::json NuDockServer::setParameters(const nlohmann::json &request) {
  if (verbose) {
    MACH3LOG_INFO("Setting parameters from NuDock request: {}", request.dump());
  }
  std::map<std::string, double> syst_params = request["sys_pars"].get<std::map<std::string, double>>();
  std::map<std::string, double> osc_params = request["osc_pars"].get<std::map<std::string, double>>();
  // Loop over systematic objects
  for (size_t s = 0; s < systematics.size(); ++s) {
    int npars = systematics[s]->GetNumParams();
    // Loop over each parameter in the systematics object
    for (size_t p = 0; p < npars; ++p) {
      std::string param_name = systematics[s]->GetParFancyName(p);
      // Check if this is an oscillation parameter
      if (NuDockOscNameMap_r.find(param_name) != NuDockOscNameMap_r.end()) {
        param_name = NuDockOscNameMap_r.at(param_name);
        // Check if it exists in the request
        if (osc_params.find(param_name) != osc_params.end()) {
          double param_value = osc_params[param_name];
          // Convert to theta to sin2_theta
          if (param_name == "Theta12" || param_name == "Theta13" || param_name == "Theta23") {
            param_value = sin(param_value) * sin(param_value);
          }
          systematics[s]->SetParCurrProp(p, param_value);
          if (verbose) MACH3LOG_INFO("Setting osc param {} to value {}", param_name, param_value);
        } else {
          MACH3LOG_WARN("Oscillation parameter {} not found in request", param_name);
        }
      } else {
        if (!systematics[s]->IsParameterFixed(p)) {
          // Non-oscillation parameter
          // Check if it exists in the request
          if (syst_params.find(param_name) != syst_params.end()) {
            double param_value = syst_params[param_name];
            systematics[s]->SetParCurrProp(p, param_value);
            if (verbose) MACH3LOG_INFO("Setting syst param {} to value {}", param_name, param_value);
          } else {
            MACH3LOG_WARN("Systematic parameter {} not found in request", param_name);
          }
        }
      }
    }
  }
  nlohmann::json response;
  response["status"] = "success";
  return response;
}

nlohmann::json NuDockServer::setAsimovPoint(const nlohmann::json &request) {
  (void)request;
  nlohmann::json response;
  response["status"] = "not implemented";
  return response;
}

nlohmann::json NuDockServer::getParametersNames(const nlohmann::json &request) {
  (void)request;
  std::vector<std::string> syst_par_names;

  for (size_t s = 0; s < systematics.size(); ++s) {
    int npars = systematics[s]->GetNumParams();
    for (size_t p = 0; p < npars; ++p) {
      std::string param_name = systematics[s]->GetParFancyName(p);
      if (NuDockOscNameMap_r.find(param_name) != NuDockOscNameMap_r.end()) continue;
      if (!systematics[s]->IsParameterFixed(p)) {
        syst_par_names.push_back(param_name);
      }
    }
  }
  nlohmann::json response;
  response["sys_pars"] = syst_par_names;
  return response;
}

nlohmann::json NuDockServer::getParameters(const nlohmann::json &request) {
  std::map<std::string, double> osc_params;
  std::map<std::string, double> syst_params;

  for (size_t s = 0; s < systematics.size(); ++s) {
    int npars = systematics[s]->GetNumParams();
    for (size_t p = 0; p < npars; ++p) {
      std::string param_name = systematics[s]->GetParFancyName(p);
      double param_value = systematics[s]->GetParProp(p);
      if (NuDockOscNameMap_r.find(param_name) != NuDockOscNameMap_r.end()){
        param_name = NuDockOscNameMap_r.at(param_name);
        osc_params[param_name] = param_value;
      } else {
        if (!systematics[s]->IsParameterFixed(p)) {
          syst_params[param_name] = param_value;
        }
      }
    }
  }
  (void)request;
  // Convert maps to json
  nlohmann::json response;
  response["syst_params"] = syst_params;
  response["osc_params"] = osc_params;
  return response;
}

nlohmann::json NuDockServer::getLogLikelihood(const nlohmann::json &request) {
  (void)request;
  nlohmann::json response;
  response["log_likelihood"] = getLogLikelihood();
  return response;
}
