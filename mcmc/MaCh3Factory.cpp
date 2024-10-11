// MaCh3 includes
#include "mcmc/MaCh3Factory.h"


// ********************************************
std::unique_ptr<FitterBase> MaCh3FitterFactory(manager *fitMan, std::vector<samplePDFBase*>& Samples, std::vector<covarianceBase*>& Covariances) {
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

  //KS: Adding samples and covariances to the Fitter class could be in the factory
  for(unsigned int i = 0; Samples.size(); i++)
    MaCh3Fitter->addSamplePDF(Samples[i]);
  for(unsigned int i = 0; Covariances.size(); i++)
    MaCh3Fitter->addSystObj(Covariances[i]);

  return MaCh3Fitter;
}

// ********************************************
covarianceXsec* MaCh3CovarianceFactory(manager *fitMan, const std::string& name) {
// ********************************************
  MACH3LOG_INFO("Initialising {} matrix", name);

  // config for our matrix
  YAML::Node Settings = fitMan->raw()["General"]["Systematics"][name];
  // yaml files initialising out matrix
  auto CovMatrixFile = Settings["CovFile"].as<std::vector<std::string>>();

  // PCA threshold, -1 means no pca
  auto PCAThreshold = GetFromManager<double>(Settings["PCAThreshold"], -1);
  // do we pca whole matrix or only submatrix
  auto PCAParamRegion = GetFromManager<std::vector<double>>(Settings["PCAParams"], {-999, -999});

  // create our glorious matrix
  covarianceXsec* xsec = new covarianceXsec(CovMatrixFile, "xsec_cov", PCAThreshold, PCAParamRegion[0], PCAParamRegion[1]);

  // Fill the parameter values with their nominal values
  // should _ALWAYS_ be done before overriding with fix or flat
  xsec->setParameters();

  auto FixParams = GetFromManager<std::vector<std::string>>(Settings["FixParams"], {});

  // Fixed xsec parameters loop
  if (FixParams.size() == 1 && FixParams.at(0) == "All") {
    for (int j = 0; j < xsec->GetNumParams(); j++) {
      xsec->toggleFixParameter(j);
    }
  } else {
    for (unsigned int j = 0; j < FixParams.size(); j++) {
      xsec->toggleFixParameter(FixParams.at(j));
    }
  }
  //Global step scale for matrix
  auto StepScale = Settings["StepScale"].as<double>();

  MACH3LOG_INFO("With Step Scale {:.4f}", StepScale);
  xsec->setStepScale(StepScale);

  // TODO
  // Adaptive MCMC stuff
  //xsec->initialiseAdaption(FitManager->raw());

  MACH3LOG_INFO("Factory successful");

  return xsec;
}
