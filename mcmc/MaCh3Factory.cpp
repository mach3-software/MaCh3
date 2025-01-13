// MaCh3 includes
#include "mcmc/MaCh3Factory.h"

// ********************************************
std::unique_ptr<FitterBase> MaCh3FitterFactory(manager *fitMan) {
// ********************************************
  std::unique_ptr<FitterBase> MaCh3Fitter = nullptr;

  auto Algorithm = GetFromManager<std::string>(fitMan->raw()["General"]["FittingAlgorithm"], "MCMC");

  if(Algorithm == "MCMC") {
    MaCh3Fitter = std::make_unique<mcmc>(fitMan);
  } else if (Algorithm == "PSO") {
    MaCh3Fitter = std::make_unique<PSO>(fitMan);
  } else if (Algorithm == "Minuit2") {
    #ifdef MaCh3_MINUIT2
    MaCh3Fitter = std::make_unique<MinuitFit>(fitMan);
    #else
    MACH3LOG_ERROR("Trying to use Minuit2 however MaCh3 was compiled without Minuit2 support");
    throw MaCh3Exception(__FILE__ , __LINE__ );
    #endif
  } else {
    MACH3LOG_ERROR("You want to use algorithm {}, I don't recognize it, sry", Algorithm);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return MaCh3Fitter;
}

// ********************************************
covarianceXsec* MaCh3CovarianceFactory(manager *FitManager, const std::string& PreFix) {
// ********************************************
  return MaCh3CovarianceFactory<covarianceXsec>(FitManager, PreFix);
}
