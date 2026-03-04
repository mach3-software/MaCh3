#include "NuDockServerBase.h"

// ***************************************************************************
/// @copydoc NuDockServerBase::NuDockServerBase
// ***************************************************************************
NuDockServerBase::NuDockServerBase(manager *man) : FitterBase(man) {
}

// ***************************************************************************
/// @copydoc NuDockServerBase::~NuDockServerBase
// ***************************************************************************
NuDockServerBase::~NuDockServerBase() {
}

// ***************************************************************************
/// @copydoc NuDockServerBase::setup
// ***************************************************************************
void NuDockServerBase::setup() {
  // Resize the likelihood storage vectors
  syst_llh.resize(systematics.size(), M3::_BAD_DOUBLE_);
  sample_llh.resize(samples.size(), M3::_BAD_DOUBLE_);
  verbose = GetFromManager(fitMan->raw()["NuDock"]["Verbose"], false);
  add_prior_llh = GetFromManager(fitMan->raw()["NuDock"]["AddPriorLLH"], false);
}

// ***************************************************************************
/// @copydoc NuDockServerBase::getLogLikelihood
// ***************************************************************************
double NuDockServerBase::getLogLikelihood() {
  // HH: Based on MR2T2::ProposeStep()
  // Initial likelihood
  double llh = 0.0;

  if (add_prior_llh) {
    // Loop over the systematics and propose the initial step
    for (size_t s = 0; s < systematics.size(); ++s) {
      // Get the likelihood from the systematics
      syst_llh[s] = systematics[s]->GetLikelihood();
      if (verbose) std::cout << "Syst " << systematics[s]->GetName() << " llh: " << syst_llh[s] << std::endl;
      llh += syst_llh[s];

      if (verbose) std::cout << "LLH after " << systematics[s]->GetName() << " " << llh << std::endl;
    }
  }

  // DB for atmospheric event by event sample migration, need to fully reweight all samples to allow event passing prior to likelihood evaluation
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

  return llh;
}

// ***************************************************************************
/// @copydoc NuDockServerBase::setParameters
// ***************************************************************************
nlohmann::json NuDockServerBase::setParameters(const nlohmann::json &request) {
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
        std::string param_name_nudock = NuDockOscNameMap_r.at(param_name);
        // Check if it exists in the request
        if (osc_params.find(param_name_nudock) != osc_params.end()) {
          double param_value = osc_params[param_name_nudock];
          FormatOscParsForMaCh3(param_name_nudock, param_value);
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
  // HH: Reweights here is probably a better idea than getLikelihood, thanks to Jude
  for (size_t i = 0; i < samples.size(); ++i) {
    samples[i]->Reweight();
  }

  nlohmann::json response;
  response["status"] = "success";
  return response;
}

// ***************************************************************************
/// @copydoc NuDockServerBase::setAsimovPoint
// ***************************************************************************
nlohmann::json NuDockServerBase::setAsimovPoint(const nlohmann::json &request) {
  (void)request;
  // reweight sample and add to data
  for (size_t ipdf=0; ipdf<samples.size(); ipdf++) {
    samples[ipdf]->Reweight();
    if (auto* fd_casted_sample = dynamic_cast<SampleHandlerFD*>(samples[ipdf])) {
      for (int iSample = 0; iSample < samples[ipdf]->GetNsamples(); ++iSample) {
        if (fd_casted_sample->GetNDim(iSample) == 1) {
          fd_casted_sample->AddData(iSample, (TH1D*)fd_casted_sample->GetMCHist(iSample, fd_casted_sample->GetNDim(iSample))->Clone());
        } else if (fd_casted_sample->GetNDim(iSample) == 2) {
          fd_casted_sample->AddData(iSample, (TH2D*)fd_casted_sample->GetMCHist(iSample, fd_casted_sample->GetNDim(iSample))->Clone());
        } else {
          MACH3LOG_ERROR("Unsupported histogram dimension for SampleHandlerFD: {}", fd_casted_sample->GetNDim(iSample));
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
    } else {
      MACH3LOG_ERROR("Sample object does not derive from SampleHandlerFD. Consider overloading setAsimovPoint for this sample type.");
      throw MaCh3Exception(__FILE__,__LINE__);
    } 
  }

  // set prefit parameter values to the current parameter values
  for (size_t s = 0; s < systematics.size(); ++s) {
    int npars = systematics[s]->GetNumParams();
    for (int p = 0; p < npars; ++p) {
      std::string param_name = systematics[s]->GetParFancyName(p);
      double param_value = systematics[s]->GetParCurr(p);
      systematics[s]->SetPar(p, param_value);
      if (verbose) MACH3LOG_INFO("Setting prefit param {} to current value {}", param_name, param_value);
    }
  }

  nlohmann::json response;
  response["status"] = "success";
  return response;
}

// ***************************************************************************
/// @copydoc NuDockServerBase::getSpectrum
// ***************************************************************************
nlohmann::json NuDockServerBase::getSpectrum(const nlohmann::json &request) {
  (void)request;
  nlohmann::json response;

  std::vector<std::string> sample_titles = {};

  for (size_t ipdf=0; ipdf<samples.size(); ipdf++) {
    if (auto* fd_casted_sample = dynamic_cast<SampleHandlerFD*>(samples[ipdf])) {
      for (int iSample = 0; iSample < samples[ipdf]->GetNsamples(); ++iSample) {
        std::string sample_title = samples[ipdf]->GetSampleTitle(iSample);
        sample_titles.push_back(sample_title);

        int dimension = fd_casted_sample->GetNDim(iSample);
        response["dimensions"][sample_title] = dimension;

        if (dimension == 1) {
          TH1D* mc_hist = (TH1D*)fd_casted_sample->GetMCHist(iSample, dimension)->Clone();

          TAxis *ax = mc_hist->GetXaxis();
          std::string xtitle = ax->GetTitle(); 
          int nxbins = ax->GetNbins();
          std::vector<double> xbins(nxbins+1); 
          std::vector<double> binvals(nxbins);

          xbins[0] = ax->GetBinLowEdge(1);
          for (int ix = 1; ix <= nxbins; ++ix) {
            xbins[ix] = ax->GetBinUpEdge(ix);
            binvals[ix-1] = mc_hist->GetBinContent(ix);
          }
          response["axis_titles"][sample_title] = {xtitle};
          response["bin_edges"][sample_title] = {xbins};
          response["bin_values"][sample_title] = binvals;
        } else if (dimension == 2) {
          TH2D* mc_hist = (TH2D*)fd_casted_sample->GetMCHist(iSample, dimension)->Clone();

          TAxis *x_axis = mc_hist->GetXaxis();
          std::string xtitle = x_axis->GetTitle();
          int nxbins = x_axis->GetNbins();
          std::vector<double> xbins(nxbins+1);

          TAxis *y_axis = mc_hist->GetYaxis();
          std::string ytitle = y_axis->GetTitle();
          int nybins = y_axis->GetNbins();
          std::vector<double> ybins(nybins+1);

          std::vector<std::vector<double>> binvals(nxbins, std::vector<double>(nybins));

          xbins[0] = x_axis->GetBinLowEdge(1);
          ybins[0] = y_axis->GetBinLowEdge(1);
          for (int ix = 1; ix <= nxbins; ++ix) {
            for (int iy = 1; iy <= nybins; ++iy) {
              xbins[ix] = x_axis->GetBinUpEdge(ix);
              ybins[iy] = y_axis->GetBinUpEdge(iy);
              binvals[ix-1][iy-1] = mc_hist->GetBinContent(ix, iy);
            }
          }
          response["axis_titles"][sample_title] = {xtitle, ytitle};
          response["bin_edges"][sample_title] = {xbins, ybins};
          response["bin_values"][sample_title] = binvals;
        }
      }
    } else {
      MACH3LOG_ERROR("Sample object does not derived from SampleHandlerFD. Consider overloading getSpectrum for this sample type.");
      throw MaCh3Exception(__FILE__,__LINE__);
    } 
  }
  response["sample_names"] = sample_titles;
  return response;
}

// ***************************************************************************
/// @copydoc NuDockServerBase::getParametersNames
// ***************************************************************************
nlohmann::json NuDockServerBase::getParametersNames(const nlohmann::json &request) {
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

// ***************************************************************************
/// @copydoc NuDockServerBase::getParameters
// ***************************************************************************
nlohmann::json NuDockServerBase::getParameters(const nlohmann::json &request) {
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
  response["sys_pars"] = syst_params;
  response["osc_pars"] = osc_params;
  return response;
}

// ***************************************************************************
/// @copydoc NuDockServerBase::get2NLL
// ***************************************************************************
nlohmann::json NuDockServerBase::get2NLL(const nlohmann::json &request) {
  (void)request;
  nlohmann::json response;
  response["log_likelihood"] = 2*getLogLikelihood();
  return response;
}
