#include "mcmc/MulticanonicalMCMCHandler.h"
#include "covariance/covarianceOsc.h"

MulticanonicalMCMCHandler::MulticanonicalMCMCHandler() {
  // Initialize member variables with defaults
  oscCovVar = -1;
  multicanonicalVar = -1;
  multicanonicalVar_dm23 = -1;
  multicanonicalSeparate = false;
  multicanonicalSpline = false;
  multicanonicalBeta = 1.0;
  multicanonicalSigma = 1.0;
  delta_cp_value = 0.0;
  delm23_value = 0.0;
  dcp_spline_IO = nullptr;
  dcp_spline_NO = nullptr;
  multicanonicalSeparateMean = 0.0;
  multicanonicalSeparateSigma = 1.0;
  umbrellaNumber = 5;
  umbrellaOverlapMode = false;
  umbrellaSigmaOverlap = 3.0;
  umbrellaAdjustStepScale = false;
  umbrellaStepScaleFactor = 1.0;
  flipWindow = false;
}

MulticanonicalMCMCHandler::~MulticanonicalMCMCHandler() {
  // Destructor - no dynamic memory to clean up (TSpline3 pointers are owned by TFile)
}

void MulticanonicalMCMCHandler::InitializeMulticanonicalHandlerConfig(manager* fitMan, std::vector<covarianceBase*>& systematics){
    
    bool foundDeltaCP = false;
    bool foundDelm23 = false;
    
    MACH3LOG_INFO("Looping over systematics to find delta_cp parameter");
    // Loop over the systematics and find the osc_cov systematic and the delta_cp parameter number
    for (size_t s = 0; s < systematics.size(); s++) {
      MACH3LOG_INFO("Systematic: {}", systematics[static_cast<int>(s)]->getName());
      if (systematics[static_cast<int>(s)]->getName() == "osc_cov") {
        oscCovVar = static_cast<int>(s);
        MACH3LOG_INFO("Found osc_cov systematic saving in variable {}", oscCovVar);
        for (int i = 0; i < systematics[static_cast<int>(s)]->GetNumParams(); i++) {
          MACH3LOG_INFO("Parameter: {}", systematics[static_cast<int>(s)]->GetParName(i));
          if (systematics[static_cast<int>(s)]->GetParName(i) == "delta_cp") {
            multicanonicalVar = i;
            MACH3LOG_INFO("Setting multicanonical weight on delta_cp parameter int {}",i);
            foundDeltaCP = true;
          }
          if (systematics[static_cast<int>(s)]->GetParName(i) == "delm2_23") {
            multicanonicalVar_dm23 = i;
            MACH3LOG_INFO("Setting delm2_23 parameter int {}",i);
            foundDelm23 = true;
          }
        }
      }
      if (!foundDeltaCP) {
        MACH3LOG_ERROR("Could not find delta_cp parameter in osc_cov systematic");
        throw std::runtime_error("Could not find delta_cp parameter in osc_cov systematic");
      }
      if (!foundDelm23) {
        MACH3LOG_ERROR("Could not find delm2_23 parameter in osc_cov systematic");
        throw std::runtime_error("Could not find delm2_23 parameter in osc_cov systematic");
      }
    }

    multicanonicalSpline = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["MulticanonicalSpline"],false);

    multicanonicalSeparate = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["MulticanonicalSeparate"],false);

    if (multicanonicalSpline && multicanonicalSeparate) {
      MACH3LOG_ERROR("Cannot use both multicanonical spline and separate method at the same time. Please choose one.");
      throw std::runtime_error("Cannot use both multicanonical spline and separate method at the same time.");
    }
    
    if (multicanonicalSpline){
      std::string splineFile = GetFromManager<std::string>(fitMan->raw()["General"]["MCMC"]["MulticanonicalSplineFile"],"nofile");
      TFile *file = new TFile(splineFile.c_str(), "READ");
      if (!file || file->IsZombie()) {
        MACH3LOG_ERROR("Could not open multicanonical spline file: {}", splineFile);
        throw std::runtime_error("Could not open multicanonical spline file");
      }
      dcp_spline_IO = static_cast<TSpline3*>(file->Get("dcp_spline_IO"));
      MACH3LOG_INFO("Using multicanonical spline from file {}", splineFile);
      dcp_spline_IO->Eval(0.0); // check that the spline is valid
      std::cout << "Spline evaluated at 0.0 gives value " << dcp_spline_IO->Eval(0.0) << std::endl;

      dcp_spline_NO = static_cast<TSpline3*>(file->Get("dcp_spline_NO"));
      MACH3LOG_INFO("Using multicanonical spline from file {}", splineFile);
      dcp_spline_NO->Eval(0.0); // check that the spline is valid
      std::cout << "Spline evaluated at 0.0 gives value " << dcp_spline_NO->Eval(0.0) << std::endl;

      // // check for empty TSpline
      // if (dcp_spline == nullptr) {
      //   MACH3LOG_ERROR("Spline not found in file {}", splineFile);
      //   throw std::runtime_error("Spline not found in file");
      // }

    } else if (multicanonicalSeparate) { 
      // If we are using the multicanonical method in separate chains, we need to get the separate mean and sigma values
      MACH3LOG_INFO("Using separate multicanonical method");
      multicanonicalSeparateMean = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["MulticanonicalSeparateMean"], -TMath::Pi());
      MACH3LOG_INFO("Setting multicanonical mean to {}", multicanonicalSeparateMean);
      umbrellaNumber = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["UmbrellaNumber"], 5);
      umbrellaOverlapMode = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["UmbrellaOverlapMode"], false);
      if (umbrellaOverlapMode) {
        MACH3LOG_INFO("Setting width based on # of sigma overlapping between umbrellas");
        umbrellaSigmaOverlap = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["UmbrellaSigmaOverlap"], 3.0);
        MACH3LOG_INFO("Setting umbrella number to {}", umbrellaNumber);
        multicanonicalSeparateSigma = TMath::Pi()/((umbrellaNumber - 1)*(umbrellaSigmaOverlap));
      } else {
        multicanonicalSeparateSigma = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["MulticanonicalSeparateSigma"],(2*TMath::Pi())/umbrellaNumber);
        MACH3LOG_INFO("Setting width based on value in config {}", multicanonicalSeparateSigma);
      } 
      // set individual step scale for dcp, so that the ratio of the step scale to the multicanonical sigma is stepscale/1sigmaerror = 1/2pi 
      umbrellaAdjustStepScale = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["UmbrellaAdjustStepScale"], false);
      umbrellaStepScaleFactor = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["UmbrellaStepScaleFactor"], 1.0);
      if (umbrellaAdjustStepScale){
        MACH3LOG_INFO("Adjusting umbrella step scale to keep ratio of step scale to multicanonical sigma constant");
        MACH3LOG_INFO("Setting umbrella step scale factor to {}", umbrellaStepScaleFactor);
        double stepScale = (multicanonicalSeparateSigma * umbrellaStepScaleFactor) / (2.0 * TMath::Pi());
        MACH3LOG_INFO("Setting individual step scale for multicanonical separate to {}", stepScale);
        // Set the individual step scale for the multicanonical variable
        for (auto& syst : systematics) {
          if (syst->getName() == "osc_cov") {
            syst->setIndivStepScale(multicanonicalVar, stepScale);
            MACH3LOG_INFO("Setting individual step scale for {} systematic to {}", multicanonicalVar, stepScale);
          }
        }       
      } else {
        MACH3LOG_INFO("Not adjusting umbrella step scale, using value in OscCov config");
      }

      // Set the flip window flag for the oscillation systematic
      flipWindow = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["FlipWindow"], false);
      MACH3LOG_INFO("Flip Window: {}", flipWindow);

    } else {
      MACH3LOG_INFO("Using multicanonical Gaussian umbrellas");
      // Get the multicanonical sigma values from the configuration file
      multicanonicalSigma = fitMan->raw()["General"]["MCMC"]["MulticanonicalSigma"].as<double>();
      MACH3LOG_INFO("Setting multicanonical sigma to {}", multicanonicalSigma);
    }

    // Get the multicanonical beta value from the configuration file
    multicanonicalBeta = fitMan->raw()["General"]["MCMC"]["MulticanonicalBeta"].as<double>();
    MACH3LOG_INFO("Setting multicanonical beta to {}", multicanonicalBeta);
}

void MulticanonicalMCMCHandler::InitializeMulticanonicalParams(std::vector<covarianceBase*>& systematics){
    if (multicanonicalSeparate){
    // set the starting point of the chain to the mean value of the multicanonical umbrella
    for (auto& syst : systematics) {
        if (syst->getName() == "osc_cov") {
            syst->printNominalCurrProp();
            syst->setParCurrProp(multicanonicalVar, multicanonicalSeparateMean);
            syst->setParProp(multicanonicalVar, multicanonicalSeparateMean);
            MACH3LOG_INFO("Setting starting point of chain to mean value for multicanonical separate: {}", multicanonicalSeparateMean);
            // pass the mean to the covarianceOsc object for parameter flipping
            if (flipWindow) {
            auto* oscCov = dynamic_cast<covarianceOsc*>(syst);
            if (oscCov) {
                oscCov->setFlipWindow(flipWindow);
                oscCov->setMulticanonicalSeparateMean(multicanonicalSeparateMean);
            }
            }
            syst->printNominalCurrProp();
        }
    }
  }
}

double MulticanonicalMCMCHandler::GetMulticanonicalWeightGaussian(double deltacp){
  // precalculate constants
  constexpr double inv_sqrt_2pi = 0.3989422804014337;
  double sigma = multicanonicalSigma;
  const double neg_half_sigma_sq = -1/(2*sigma*sigma);
  // three gaussians centered at -pi, 0, pi with sigma pre-defined above
  double exp1 = std::exp(neg_half_sigma_sq * (deltacp - TMath::Pi()) * (deltacp - TMath::Pi()));
  double exp2 = std::exp(neg_half_sigma_sq * (deltacp) * (deltacp));
  double exp3 = std::exp(neg_half_sigma_sq * (deltacp + TMath::Pi()) * (deltacp + TMath::Pi()));
  ///delta_cp_log_likelihood = -TMath::Log(TMath::Gaus(deltacp,TMath::Pi(),1,kTRUE)+TMath::Gaus(deltacp,0,1,kTRUE)+TMath::Gaus(deltacp,-TMath::Pi(),1,kTRUE));

  // return the log likelihood, ie the log of the normalised sum of the gaussians
  return -std::log(inv_sqrt_2pi * (1/sigma) * (exp1 + exp2 + exp3))*(multicanonicalBeta);
}

double MulticanonicalMCMCHandler::GetMulticanonicalWeightSpline(double deltacp, double delm23){
  double dcp_spline_val;

  if (delm23 < 0){
    dcp_spline_val = dcp_spline_IO->Eval(deltacp);
    return -(-std::log(dcp_spline_val)+std::log(dcp_spline_IO->Eval(-TMath::Pi()/2)))*(multicanonicalBeta); // do I want this offset?? does it matter?
  } else {
    dcp_spline_val = dcp_spline_NO->Eval(deltacp);
    return -(-std::log(dcp_spline_val)+std::log(dcp_spline_NO->Eval(-TMath::Pi()/2)))*(multicanonicalBeta);
  }
  // std::cout << "Evaluating spline at delta_cp = " << deltacp << " gives value " << dcp_spline_val << "with -log lh of :" << -log(dcp_spline_val) << std::endl;
}

double MulticanonicalMCMCHandler::GetMulticanonicalWeightSeparate(double deltacp){
  // precalculate constants
  constexpr double inv_sqrt_2pi = 0.3989422804014337;
  const double neg_half_sigma_sq = -1/(2*multicanonicalSeparateSigma*multicanonicalSeparateSigma); 
  // return the log likelihood, ie the log of the normalised gaussian
  return (-std::log(inv_sqrt_2pi * (1/multicanonicalSeparateSigma) * std::exp(neg_half_sigma_sq * (deltacp - multicanonicalSeparateMean) * (deltacp - multicanonicalSeparateMean))))*(multicanonicalBeta);
}

