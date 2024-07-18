// MaCh3 includes
#include "mcmc/MaCh3Factory.h"

FitterBase* MaCh3FitterFactory(manager *fitMan, std::vector<samplePDFBase*>& Samples, std::vector<covarianceBase*>& Covariances)
{
  FitterBase* MaCh3Fitter = nullptr;

  std::string Algorithm = GetFromManager<std::string>(fitMan->raw()["General"]["FittingAlgorithm"], "MCMC");

  if(Algorithm == "MCMC") {
    MaCh3Fitter = new mcmc(fitMan);
  } else if (Algorithm == "PSO") {
    MaCh3Fitter = new PSO(fitMan);
  } else if (Algorithm == "Minuit2") {
  #ifdef MaCh3_MINUIT2
    MaCh3Fitter = new MinuitFit(fitMan);
  #else
    MACH3LOG_ERROR("Trying to use Minuit2 however MaCh3 was compiled without Minuit2 support");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  #endif
  } else {
    MACH3LOG_ERROR("You want to use algorithm {}, I don't recognize it, sry", Algorithm);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  //KS: Adding samples and covariances to the Fitter class could be in the factory
  for(unsigned int i = 0; Samples.size(); i++)
    MaCh3Fitter->addSamplePDF(Samples[i]);
  for(unsigned int i = 0; Covariances.size(); i++)
    MaCh3Fitter->addSystObj(Covariances[i]);

  return MaCh3Fitter;
}
