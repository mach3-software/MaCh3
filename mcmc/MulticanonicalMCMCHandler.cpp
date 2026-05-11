#include "mcmc/MulticanonicalMCMCHandler.h"
#include "covariance/covarianceOsc.h"
#include <stdexcept>

MulticanonicalMCMCHandler::BiasFunction ParseBiasFunction(const std::string& biasFunctionName) {
  if (biasFunctionName == "gaussian") {
    return MulticanonicalMCMCHandler::BiasFunction::Gaussian;
  }
  if (biasFunctionName == "vonMises") {
    return MulticanonicalMCMCHandler::BiasFunction::VonMises;
  }
  if (biasFunctionName == "generalisedGaussian") {
    return MulticanonicalMCMCHandler::BiasFunction::GeneralisedGaussian;
  }

  throw std::runtime_error("Unknown multicanonical bias function: " + biasFunctionName);
}

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
  multicanonicalGenGaussianMean = 0.0;
  multicanonicalGenGaussianWidth = 0.0;
  umbrellaNumber = 5;
  umbrellaOverlapMode = false;
  umbrellaSigmaOverlap = 3.0;
  umbrellaAdjustStepScale = false;
  umbrellaStepScaleFactor = 1.0;
  flipWindow = false;

  vonMises_mode = false;
  vonMises_kappa = -1.0;
  vonMises_I0_kappa = -1.0;
  multicanonicalBiasFunction = BiasFunction::Gaussian;
  multicanonicalBiasFunctionName = "gaussian";
}

MulticanonicalMCMCHandler::~MulticanonicalMCMCHandler() {
  // Destructor 
}

void MulticanonicalMCMCHandler::InitializeMulticanonicalHandlerConfig(manager* fitMan, std::vector<covarianceBase*>& systematics){
    
    bool foundDeltaCP = false;
    bool foundDelm23 = false;
    
    MACH3LOG_INFO("Looping over systematics to find delta_cp parameter");
    // Loop over the systematics and find the osc_cov systematic and the delta_cp parameter number
    // TODO make this its own grab parameter functions, and use ENUMS instead
    for (size_t s = 0; s < systematics.size(); s++) {
      MACH3LOG_INFO("Systematic: {}", systematics[static_cast<int>(s)]->getName());
      
      // if we find the covariance object
      if (systematics[static_cast<int>(s)]->getName() == "osc_cov") {
        oscCovVar = static_cast<int>(s);
        MACH3LOG_INFO("Found osc_cov systematic saving in variable {}", oscCovVar);
        
	for (int i = 0; i < systematics[static_cast<int>(s)]->GetNumParams(); i++) {
          MACH3LOG_INFO("Parameter: {}", systematics[static_cast<int>(s)]->GetParName(i));
          
	  // check for params of interest
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

      // if we didn't find both parameters we need throw
      if (!foundDeltaCP) {
        MACH3LOG_ERROR("Could not find delta_cp parameter in osc_cov systematic");
        throw std::runtime_error("Could not find delta_cp parameter in osc_cov systematic");
      }
      if (!foundDelm23) {
        MACH3LOG_ERROR("Could not find delm2_23 parameter in osc_cov systematic");
        throw std::runtime_error("Could not find delm2_23 parameter in osc_cov systematic");
      }
    }
    const auto mcmcConfig = fitMan->raw()["General"]["MCMC"];
    
    multicanonicalSpline = GetFromManager<bool>(mcmcConfig["MulticanonicalSpline"],false);

    multicanonicalSeparate = GetFromManager<bool>(mcmcConfig["MulticanonicalSeparate"],false);

    const std::string biasFunctionName = GetFromManager<std::string>(mcmcConfig["MulticanonicalBiasFunction"], "");
    const bool hasBiasFunction = !biasFunctionName.empty();

    if (multicanonicalSpline && (multicanonicalSeparate || hasBiasFunction)) {
      MACH3LOG_ERROR("Cannot use multicanonical spline together with umbrella bias function selection.");
      throw std::runtime_error("Cannot use multicanonical spline together with umbrella bias function selection.");
    }

    if (hasBiasFunction) {
      multicanonicalBiasFunction = ParseBiasFunction(biasFunctionName);
      multicanonicalBiasFunctionName = biasFunctionName;
      multicanonicalSeparate = true;
    } else if (multicanonicalSeparate) {
      if (GetFromManager<bool>(mcmcConfig["VonMisesMode"], false)) {
        multicanonicalBiasFunction = BiasFunction::VonMises;
        multicanonicalBiasFunctionName = "vonMises";
      } else if (GetFromManager<bool>(mcmcConfig["MulticanonicalGenGaussian"], false)) {
        multicanonicalBiasFunction = BiasFunction::GeneralisedGaussian;
        multicanonicalBiasFunctionName = "genGaussian";
      } else {
        multicanonicalBiasFunction = BiasFunction::Gaussian;
        multicanonicalBiasFunctionName = "gaussian";
      }
    }

    MACH3LOG_INFO("Using multicanonical bias function {}", multicanonicalBiasFunctionName);

    // setup for spline bias mode
    if (multicanonicalSpline){

      std::string splineFile = GetFromManager<std::string>(fitMan->raw()["General"]["MCMC"]["MulticanonicalSplineFile"],"nofile");
      
      TFile *file = new TFile(splineFile.c_str(), "READ");
      if (!file || file->IsZombie()) {
        MACH3LOG_ERROR("Could not open multicanonical spline file: {}", splineFile);
        throw std::runtime_error("Could not open multicanonical spline file");
      }
      
      // grab the splines and do a quick check that they are evaluatable
      dcp_spline_IO = static_cast<TSpline3*>(file->Get("dcp_spline_IO"));
      MACH3LOG_INFO("Using multicanonical spline from file {}", splineFile);
      dcp_spline_IO->Eval(0.0); // check that the spline is valid
      MACH3LOG_INFO("Spline evaluated at 0.0 gives value: {}",dcp_spline_IO->Eval(0.0));

      dcp_spline_NO = static_cast<TSpline3*>(file->Get("dcp_spline_NO"));
      MACH3LOG_INFO("Using multicanonical spline from file {}", splineFile);
      dcp_spline_NO->Eval(0.0); // check that the spline is valid
      MACH3LOG_INFO("Spline evaluated at 0.0 gives value {}",dcp_spline_NO->Eval(0.0));

      // // check for empty TSpline
      // if (dcp_spline == nullptr) {
      //   MACH3LOG_ERROR("Spline not found in file {}", splineFile);
      //   throw std::runtime_error("Spline not found in file");
      // }

    } else if (multicanonicalSeparate) { 
      // If we are using the multicanonical method in separate chains, we need to get the separate mean and sigma values
      MACH3LOG_INFO("Using umbrella multicanonical method");
      multicanonicalSeparateMean = GetFromManager<double>(mcmcConfig["MulticanonicalMean"], GetFromManager<double>(mcmcConfig["MulticanonicalSeparateMean"], 0));
      multicanonicalGenGaussianMean = GetFromManager<double>(mcmcConfig["MulticanonicalMean"], multicanonicalSeparateMean);
      
      MACH3LOG_INFO("Setting multicanonical mean to {}", multicanonicalSeparateMean);
      umbrellaNumber = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["UmbrellaNumber"], 5);
      
      // dynamically adjust the width of the gaussians to ensure a certain level of overlap? be careful, not enough windows here will lead to posterior (rather than bias) dominated results
      umbrellaOverlapMode = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["UmbrellaOverlapMode"], false);
      if (umbrellaOverlapMode) {
        MACH3LOG_INFO("Setting width based on # of sigma overlapping between umbrellas");
        umbrellaSigmaOverlap = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["UmbrellaSigmaOverlap"], 3.0);
        MACH3LOG_INFO("Setting umbrella number to {}", umbrellaNumber);
        multicanonicalSeparateSigma = TMath::Pi()/((umbrellaNumber - 1)*(umbrellaSigmaOverlap));
        multicanonicalGenGaussianWidth = multicanonicalSeparateSigma;
      } else {
      	// just grab the width directly from the config
        multicanonicalSeparateSigma = GetFromManager<double>(mcmcConfig["MulticanonicalWidth"], GetFromManager<double>(mcmcConfig["MulticanonicalSeparateSigma"], (2*TMath::Pi())/umbrellaNumber));
        multicanonicalGenGaussianWidth = GetFromManager<double>(mcmcConfig["MulticanonicalWidth"], GetFromManager<double>(mcmcConfig["MulticanonicalGenGaussianWidth"], multicanonicalSeparateSigma));
        MACH3LOG_INFO("Setting width based on value in config {}", multicanonicalSeparateSigma);
      } 

      multicanonicalSigma = multicanonicalSeparateSigma;

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

      // Set the flip window flag for the oscillation systematic DO NOT USE THIS 
      // flipWindow = GetFromManager<bool>(fitMan->raw()["General"]["MCMC"]["FlipWindow"], false);
      // MACH3LOG_INFO("Flip Window: {}", flipWindow);

    } else {
      // this method is out of date and not much use, remove
      MACH3LOG_INFO("Using multicanonical Gaussian umbrellas");
      // Get the multicanonical sigma values from the configuration file
      multicanonicalSigma = fitMan->raw()["General"]["MCMC"]["MulticanonicalSigma"].as<double>();
      MACH3LOG_INFO("Setting multicanonical sigma to {}", multicanonicalSigma);
    }

    // initialize von Mises parameters if needed
    if (multicanonicalBiasFunction == BiasFunction::VonMises) {
      double temp_vonMises_sigma;
      temp_vonMises_sigma = GetFromManager<double>(mcmcConfig["MulticanonicalWidth"], GetFromManager<double>(mcmcConfig["VonMisesSigma"], -1.0));
      multicanonicalSigma = temp_vonMises_sigma;
      vonMises_kappa = 1.0 / (temp_vonMises_sigma * temp_vonMises_sigma);
      vonMises_I0_kappa = TMath::BesselI0(vonMises_kappa);
      MACH3LOG_INFO("Using von Mises distribution with kappa = {} and I0(kappa) = {}", vonMises_kappa, vonMises_I0_kappa);
    }

    if (multicanonicalBiasFunction == BiasFunction::GeneralisedGaussian) {
      multicanonicalGenGaussianWidth = GetFromManager<double>(mcmcConfig["MulticanonicalWidth"], GetFromManager<double>(mcmcConfig["MulticanonicalGenGaussianWidth"], multicanonicalSeparateSigma));
      multicanonicalGenGaussianMean = GetFromManager<double>(mcmcConfig["MulticanonicalMean"], multicanonicalSeparateMean);
      multicanonicalSigma = multicanonicalGenGaussianWidth;
      MACH3LOG_INFO("Using generalised Gaussian with mean {} and width {}", multicanonicalGenGaussianMean, multicanonicalGenGaussianWidth);
    }

    // Get the multicanonical beta value from the configuration file
    multicanonicalBeta = mcmcConfig["MulticanonicalBeta"].as<double>();
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
        // pass the mean to the covarianceOsc object for parameter flipping DO NOT USE
        //if (flipWindow) {
        //  auto* oscCov = dynamic_cast<covarianceOsc*>(syst);
        //  if (oscCov) {
        //    oscCov->setFlipWindow(flipWindow);
        //    oscCov->setMulticanonicalSeparateMean(multicanonicalSeparateMean);
        //  }
        //}
        syst->printNominalCurrProp();
      }
    }
  }
}

double MulticanonicalMCMCHandler::GetMulticanonicalWeightVonMises(double deltacp){
  // calculate the Log form of the von Mises instead to avoid numerical issues and return directly
  double log_vonMises = vonMises_kappa * std::cos(deltacp - multicanonicalSeparateMean) - std::log(2 * TMath::Pi() * vonMises_I0_kappa);
  // return the log likelihood, ie the log of the normalised von Mises
  return -log_vonMises * (multicanonicalBeta);
}

// this now sorts through the available bias functions in a single function
double MulticanonicalMCMCHandler::GetMulticanonicalWeight(double deltacp, double delm23){
  if (multicanonicalSpline) {
    return GetMulticanonicalWeightSpline(deltacp, delm23);
  }

  switch (multicanonicalBiasFunction) {
    case BiasFunction::Gaussian:
      return GetMulticanonicalWeightGaussian(deltacp);
    case BiasFunction::VonMises:
      return GetMulticanonicalWeightVonMises(deltacp);
    case BiasFunction::GeneralisedGaussian:
      return GetMulticanonicalWeightGenGaussian(deltacp);
  }

  return GetMulticanonicalWeightGaussian(deltacp);
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

double MulticanonicalMCMCHandler::GetMulticanonicalWeightGaussian(double deltacp){
  constexpr double inv_sqrt_2pi = 0.3989422804014337;
  const double neg_half_sigma_sq = -1/(2*multicanonicalSeparateSigma*multicanonicalSeparateSigma); 
  // return the log likelihood, ie the log of the normalised gaussian
  return (-std::log(inv_sqrt_2pi * (1/multicanonicalSeparateSigma) * std::exp(neg_half_sigma_sq * (deltacp - multicanonicalSeparateMean) * (deltacp - multicanonicalSeparateMean))))*(multicanonicalBeta);
}

double MulticanonicalMCMCHandler::generalisedGaussian2(double x, double mean, double width){
  constexpr int n = 2; // this controls the tightness of the gaussian fixed at 2 for now due to normalisation
  const double normFactor =1/((0.906402477055)*2*std::sqrt(2)*width); //the normalisation is a little ugly (uses gamma functions), im just going to hardcode them for now
  double likelihood = normFactor*std::exp(std::pow(-(std::pow(x-mean,2)/(2*std::pow(width,2))),n));
  return likelihood;
}

double MulticanonicalMCMCHandler::GetMulticanonicalWeightGenGaussian(double deltacp){
  // implemenetation of the generalised gaussian as a bias function
  // for now with a fixed n = 2 for simplicity

  double g0 = generalisedGaussian2(deltacp,multicanonicalGenGaussianMean,multicanonicalGenGaussianWidth);
  double g1 = generalisedGaussian2(deltacp,multicanonicalGenGaussianMean - 2*TMath::Pi(),multicanonicalGenGaussianWidth); // these two repeats are required for wrapping the gaussian around -pi and pi
  double g2 = generalisedGaussian2(deltacp,multicanonicalGenGaussianMean + 2*TMath::Pi(),multicanonicalGenGaussianWidth); 
  return -std::log(g0 + g1 + g2)*(multicanonicalBeta);
}

double MulticanonicalMCMCHandler::GetMulticanonicalWeightTripleGaussian(double deltacp){ // pretty much deprecated at this point, just here for testing
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

