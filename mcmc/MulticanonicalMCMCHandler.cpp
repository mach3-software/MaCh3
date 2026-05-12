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
  multicanonicalSpline = false;
  multicanonicalBeta = 1.0;
  multicanonicalSigma = 1.0;
  delta_cp_value = 0.0;
  delm23_value = 0.0;
  dcp_spline_IO = nullptr;
  dcp_spline_NO = nullptr;
  umbrellaMean = 0.0;
  umbrellaWidth = 1.0;
  umbrellaNumber = 5;
  umbrellaOverlapMode = false;
  umbrellaSigmaOverlap = 3.0;
  umbrellaAdjustStepScale = false;
  umbrellaStepScaleFactor = 1.0;
  flipWindow = false;

  vonMises_kappa = -1.0;
  vonMises_I0_kappa = -1.0;
  umbrellaBiasFunction = BiasFunction::Gaussian;
  umbrellaBiasFunctionName = "gaussian";
}

MulticanonicalMCMCHandler::~MulticanonicalMCMCHandler() {
  // Destructor 
}

void MulticanonicalMCMCHandler::FindOscCovParams(const std::vector<covarianceBase*>& systematics){
  MACH3LOG_INFO("Looping over systematics to find delta_cp parameter"); 
  
  bool foundDeltaCP = false;
  bool foundDelm23 = false;
  
  // Loop over the systematics and find the osc_cov systematic and the delta_cp parameter number
  for (size_t s = 0; s < systematics.size(); s++) {

    auto* syst = systematics[static_cast<int>(s)];

    MACH3LOG_INFO("Systematic: {}", syst->getName());
    
    // if we find the covariance object
    if (syst->getName() == "osc_cov") {
      oscCovVar = static_cast<int>(s);
      MACH3LOG_INFO("Found osc_cov systematic saving in variable {}", oscCovVar);
      
      for (int i = 0; i < syst->GetNumParams(); i++) {
        MACH3LOG_INFO("Parameter: {}", syst->GetParName(i));
              
        // check for params of interest
        if (syst->GetParName(i) == "delta_cp") {
          multicanonicalVar = i;
          MACH3LOG_INFO("Setting multicanonical weight on delta_cp parameter int {}",i);
          foundDeltaCP = true;
        }
        if (syst->GetParName(i) == "delm2_23") {
          multicanonicalVar_dm23 = i;
          MACH3LOG_INFO("Setting delm2_23 parameter int {}",i);
          foundDelm23 = true;
        }
      }
    }
  }

  // if we didn't find both parameters we need to throw
  if (!foundDeltaCP) {
    MACH3LOG_ERROR("Could not find delta_cp parameter in osc_cov systematic");
    throw std::runtime_error("Could not find delta_cp parameter in osc_cov systematic");
  }
  if (!foundDelm23) {
    MACH3LOG_ERROR("Could not find delm2_23 parameter in osc_cov systematic");
    throw std::runtime_error("Could not find delm2_23 parameter in osc_cov systematic");
  }

  return;
}

void MulticanonicalMCMCHandler::InitializeMulticanonicalHandlerConfig(manager* fitMan, std::vector<covarianceBase*>& systematics){
  
  FindOscCovParams(systematics);
  
  const auto mcmcConfig = fitMan->raw()["General"]["MCMC"];
  
  // TODO DR: I realised this will fail if you set spline to true but don't delete the umbrella section, thats not ideal
  multicanonicalSpline = GetFromManager<bool>(mcmcConfig["Multicanonical"]["Spline"]["SplineMode"],false);
  
  const std::string biasFunctionName = GetFromManager<std::string>(mcmcConfig["Multicanonical"]["Umbrella"]["UmbrellaBiasFunction"], "");
  const bool hasBiasFunction = !biasFunctionName.empty();

  // Spline and umbrella modes are mutually exclusive
  if (multicanonicalSpline && hasBiasFunction) {
    MACH3LOG_ERROR("Cannot use multicanonical spline together with umbrella bias function selection.");
    throw std::runtime_error("Cannot use multicanonical spline together with umbrella bias function selection.");
  }

  // Umbrella mode requires an explicit bias function selection
  if (!multicanonicalSpline && !hasBiasFunction) {
    MACH3LOG_ERROR("Multicanonical umbrella mode requires UmbrellaBiasFunction to be set (gaussian, vonMises, or generalisedGaussian).");
    throw std::runtime_error("Multicanonical umbrella mode requires UmbrellaBiasFunction to be set.");
  }

  // Parse and set the umbrella bias function enum
  if (hasBiasFunction) {
    umbrellaBiasFunction = ParseBiasFunction(biasFunctionName);
    umbrellaBiasFunctionName = biasFunctionName;
    MACH3LOG_INFO("Using umbrella bias function {}", umbrellaBiasFunctionName);
  }

  // setup for spline bias mode
  if (multicanonicalSpline){

    std::string splineFile = GetFromManager<std::string>(mcmcConfig["Multicanonical"]["Spline"]["SplineFile"],"nofile");
    
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

  } else {
    // Umbrella mode with explicit bias function selection
    MACH3LOG_INFO("Using umbrella multicanonical method with bias function {}", umbrellaBiasFunctionName);
    umbrellaMean = GetFromManager<double>(mcmcConfig["Multicanonical"]["Umbrella"]["UmbrellaMean"], 0);
    MACH3LOG_INFO("Setting multicanonical mean to {}", umbrellaMean);
    
    umbrellaNumber = GetFromManager<int>(mcmcConfig["Multicanonical"]["Umbrella"]["UmbrellaNumber"], 5);
    
    // dynamically adjust the width of the gaussians to ensure a certain level of overlap? be careful, not enough windows here will lead to posterior (rather than bias) dominated results
    umbrellaOverlapMode = GetFromManager<bool>(mcmcConfig["Multicanonical"]["Umbrella"]["AutoOverlapMode"], false);
    
    if (umbrellaOverlapMode) {
      MACH3LOG_INFO("Setting width based on # of sigma overlapping between umbrellas");
      umbrellaSigmaOverlap = GetFromManager<double>(mcmcConfig["Multicanonical"]["Umbrella"]["SigmaOverlap"], 3.0);
      MACH3LOG_INFO("Setting umbrella number to {}", umbrellaNumber);
      umbrellaWidth = TMath::Pi()/((umbrellaNumber - 1)*(umbrellaSigmaOverlap));
    } else {
      // just grab the width directly from the config
      umbrellaWidth = GetFromManager<double>(mcmcConfig["Multicanonical"]["Umbrella"]["UmbrellaWidth"], (2*TMath::Pi())/umbrellaNumber);
      MACH3LOG_INFO("Setting width based on value in config {}", umbrellaWidth);
    } 

    multicanonicalSigma = umbrellaWidth;

    // set individual step scale for dcp, so that the ratio of the step scale to the multicanonical sigma is stepscale/1sigmaerror = 1/2pi 
    umbrellaAdjustStepScale = GetFromManager<bool>(mcmcConfig["Multicanonical"]["Umbrella"]["UmbrellaAdjustStepScale"], false);
    umbrellaStepScaleFactor = GetFromManager<double>(mcmcConfig["Multicanonical"]["Umbrella"]["UmbrellaStepScaleFactor"], 1.0);
    AdjustUmbrellaStepScale(systematics);

    // Set the flip window flag for the oscillation systematic DO NOT USE THIS 
    // flipWindow = GetFromManager<bool>(mcmcConfig["Multicanonical"]["Umbrella"]["FlipWindow"], false);
    // MACH3LOG_INFO("Flip Window: {}", flipWindow);

  }

  // initialize von Mises parameters
  if (umbrellaBiasFunction == BiasFunction::VonMises) {
    double temp_vonMises_sigma;
    temp_vonMises_sigma = GetFromManager<double>(mcmcConfig["Multicanonical"]["Umbrella"]["UmbrellaWidth"], umbrellaWidth);
    multicanonicalSigma = temp_vonMises_sigma;
    vonMises_kappa = 1.0 / (temp_vonMises_sigma * temp_vonMises_sigma);
    vonMises_I0_kappa = TMath::BesselI0(vonMises_kappa);
    MACH3LOG_INFO("Using von Mises distribution with kappa = {} and I0(kappa) = {}", vonMises_kappa, vonMises_I0_kappa);
  }

  if (umbrellaBiasFunction == BiasFunction::GeneralisedGaussian) {
    multicanonicalSigma = umbrellaWidth;
    MACH3LOG_INFO("Using generalised Gaussian with mean {} and width {}", umbrellaMean, umbrellaWidth);
  }

  // Get the multicanonical beta value from the configuration file
  multicanonicalBeta = mcmcConfig["Multicanonical"]["Beta"].as<double>();
  MACH3LOG_INFO("Setting multicanonical beta to {}", multicanonicalBeta);
}

void MulticanonicalMCMCHandler::AdjustUmbrellaStepScale(const std::vector<covarianceBase*>& systematics){
  if (umbrellaAdjustStepScale){
    MACH3LOG_INFO("Adjusting umbrella step scale to keep ratio of step scale to multicanonical sigma constant");
    MACH3LOG_INFO("Setting umbrella step scale factor to {}", umbrellaStepScaleFactor);
    double stepScale = (umbrellaWidth * umbrellaStepScaleFactor) / (2.0 * TMath::Pi());
    MACH3LOG_INFO("Setting individual step scale for multicanonical separate to {}", stepScale);
    for (auto& syst : systematics) {
      if (syst->getName() == "osc_cov") {
        syst->setIndivStepScale(multicanonicalVar, stepScale);
        MACH3LOG_INFO("Setting individual step scale for {} systematic to {}", multicanonicalVar, stepScale);
      }
    }       
  } else {
    MACH3LOG_INFO("Not adjusting umbrella step scale, using value in OscCov config");
  }
}

void MulticanonicalMCMCHandler::InitializeMulticanonicalParams(std::vector<covarianceBase*>& systematics){
  // Set starting point of chain to umbrella center for non-spline umbrella modes
  if (!multicanonicalSpline) {
    for (auto& syst : systematics) {
      if (syst->getName() == "osc_cov") {
        syst->printNominalCurrProp();
        syst->setParCurrProp(multicanonicalVar, umbrellaMean);
        syst->setParProp(multicanonicalVar, umbrellaMean);
        MACH3LOG_INFO("Setting starting point of chain to mean value for multicanonical separate: {}", umbrellaMean);
        // pass the mean to the covarianceOsc object for parameter flipping DO NOT USE
        //if (flipWindow) {
        //  auto* oscCov = dynamic_cast<covarianceOsc*>(syst);
        //  if (oscCov) {
        //    oscCov->setFlipWindow(flipWindow);
        //    oscCov->setMulticanonicalSeparateMean(umbrellaMean);
        //  }
        //}
        syst->printNominalCurrProp();
        MACH3LOG_INFO("Setting starting point of chain to umbrella center: {}", umbrellaMean);
      }
    }
  }
}

double MulticanonicalMCMCHandler::GetMulticanonicalWeightVonMises(double deltacp){
  // calculate the Log form of the von Mises instead to avoid numerical issues and return directly
  double log_vonMises = vonMises_kappa * std::cos(deltacp - umbrellaMean) - std::log(2 * TMath::Pi() * vonMises_I0_kappa);
  // return the log likelihood, ie the log of the normalised von Mises
  return -log_vonMises * (multicanonicalBeta);
}

// this now sorts through the available bias functions in a single function
double MulticanonicalMCMCHandler::GetMulticanonicalWeight(double deltacp, double delm23){
  if (multicanonicalSpline) {
    return GetMulticanonicalWeightSpline(deltacp, delm23);
  }

  switch (umbrellaBiasFunction) {
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
  const double neg_half_sigma_sq = -1/(2*umbrellaWidth*umbrellaWidth); 
  // return the log likelihood, ie the log of the normalised gaussian
  return (-std::log(inv_sqrt_2pi * (1/umbrellaWidth) * std::exp(neg_half_sigma_sq * (deltacp - umbrellaMean) * (deltacp - umbrellaMean))))*(multicanonicalBeta);
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

  double g0 = generalisedGaussian2(deltacp,umbrellaMean,umbrellaWidth);
  double g1 = generalisedGaussian2(deltacp,umbrellaMean - 2*TMath::Pi(),umbrellaWidth); // these two repeats are required for wrapping the gaussian around -pi and pi
  double g2 = generalisedGaussian2(deltacp,umbrellaMean + 2*TMath::Pi(),umbrellaWidth); 
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

