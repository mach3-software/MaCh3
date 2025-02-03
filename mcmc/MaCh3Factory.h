#pragma once

// MaCh3 includes
#include "mcmc/FitterBase.h"
#include "mcmc/mcmc.h"
#include "mcmc/PSO.h"
#include "mcmc/LikelihoodFit.h"
#ifdef MaCh3_MINUIT2
#include "mcmc/MinuitFit.h"
#endif

#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"

/// @file MaCh3Factory.h
/// @brief Factory methods for MaCh3 software which streamline initialisation of different objects
/// @author Kamil Skwarczynski

/// @brief MaCh3 Factory initiates one of implemented fitting algorithms
/// @param fitMan pointer to Manager class
/// @param Samples vector of Sample PDF which will be used in the fits
/// @param Covariances vector of Systematic objects which will be used in the fits
///
/// @note Example YAML configuration:
/// @code
/// General:
///   FittingAlgorithm: ["MCMC"]
std::unique_ptr<FitterBase> MaCh3FitterFactory(manager *fitMan);

/// @brief Initializes the config manager class and allows overriding settings via command-line arguments.
/// @param argc number of arguments
/// @param argv name of arguments
/// @return A unique pointer to the initialized `manager` instance with optional overrides applied.
/// @example Usage examples:
/// ```
/// ./bin/MCMCTutorial Inputs/FitterConfig.yaml General:OutputFile:blarb.root
/// ./bin/MCMCTutorial Inputs/FitterConfig.yaml General:OutputFile:blarb.root General:MCMC:NSteps:50000
/// ```
std::unique_ptr<manager> MaCh3ManagerFactory(int argc, char **argv);

/// @brief Factory function for creating a covariance class for systematic handling.
covarianceXsec* MaCh3CovarianceFactory(manager *FitManager, const std::string& PreFix);


// ********************************************
/// @brief Factory function for creating a covariance class for systematic handling.
///
/// @param fitMan Pointer to the manager class that holds the configuration settings.
/// @param name Prefix, for example Xsec, then code will look for XsecCovFile
/// @return Pointer to the initialized covarianceXsec matrix object.
///
/// @note Example YAML configuration:
/// @code
///MaCh3CovarianceFactory(, Xsec)
///
/// General:
///   Systematics:
///       XsecCovName: "xsec_cov"
///       XsecCovFile: ["inputs/blarb1.yaml",
///                     "inputs/blarb2.yaml"]
///       XsecFix: ["Param_0",
///                 "Param_1"]
///       XsecPCAThreshold: -1
///       #XsecPCAThreshold: 0.00001
///
///       XsecPCAParams: [-999, -999]
///       XsecStepScale: 0.0075
/// @endcode
///
/// @todo add adaptive stuff
template <typename CovType>
CovType* MaCh3CovarianceFactory(manager *FitManager, const std::string& PreFix){
// ********************************************
  // config for our matrix
  YAML::Node Settings = FitManager->raw()["General"]["Systematics"];
  auto CovMatrixName = Settings[std::string(PreFix) + "CovName"].as<std::string>();
  MACH3LOG_INFO("Initialising {} matrix", CovMatrixName);

  // yaml files initialising out matrix
  auto CovMatrixFile = Settings[std::string(PreFix) + "CovFile"].as<std::vector<std::string>>();

  // PCA threshold, -1 means no pca
  auto PCAThreshold = GetFromManager<int>(Settings[std::string(PreFix) + "PCAThreshold"], -1);
  // do we pca whole matrix or only submatrix
  auto PCAParamRegion = GetFromManager<std::vector<int>>(Settings[std::string(PreFix) + "PCAParams"], {-999, -999});

  CovType* CovObject = new CovType(CovMatrixFile, CovMatrixName, PCAThreshold, PCAParamRegion[0], PCAParamRegion[1]);

  // Fill the parameter values with their nominal values
  // should _ALWAYS_ be done before overriding with fix or flat
  CovObject->setParameters();

  auto FixParams = GetFromManager<std::vector<std::string>>(Settings[std::string(PreFix) + "Fix"], {});

  // Fixed CovObject parameters loop
  if (FixParams.size() == 1 && FixParams.at(0) == "All") {
    for (int j = 0; j < CovObject->GetNumParams(); j++) {
      CovObject->toggleFixParameter(j);
    }
  } else {
    for (unsigned int j = 0; j < FixParams.size(); j++) {
      CovObject->toggleFixParameter(FixParams.at(j));
    }
  }
  //Global step scale for matrix
  auto StepScale = Settings[std::string(PreFix) + "StepScale"].as<double>();

  CovObject->setStepScale(StepScale);

  // Adaptive MCMC stuff
  if(FitManager->raw()["AdaptionOptions"])
    CovObject->initialiseAdaption(FitManager->raw());

  MACH3LOG_INFO("Factory successful");

  return CovObject;
}

// ********************************************
/// @brief Factory function for creating SamplePDF and initialisation with systematic.
///
/// @tparam SampleType The class type of the sample to create, e.g., `samplePDFTutorial`.
/// @param SampleConfig Path to sample config.
/// @param xsec A pointer to a covarianceXsec object for cross-section systematic settings.
/// @param osc (Optional) A pointer to a covarianceOsc object for oscillation systematic settings.
/// @return Vector of SampleType object, initialized and ready for use.
///
/// @note Example
/// ```cpp
/// auto mySamples = MaCh3SamplePDFFactory<samplePDFTutorial>(SampleConfig, xsec, osc);
/// ```
template <typename SampleType>
std::vector<SampleType*> MaCh3SamplePDFFactory(const std::vector<std::string>& SampleConfig,
                                               covarianceXsec* xsec,
                                               covarianceOsc* osc = nullptr) {
// ********************************************
  std::vector<SampleType*> PDFs(SampleConfig.size());
  for (size_t i = 0; i < SampleConfig.size(); ++i)
  {
    // Instantiate the sample using the specified class type
    SampleType* Sample = new SampleType(SampleConfig[i], xsec, osc);
    Sample->reweight();

    // Obtain sample name and create a TString version for histogram naming
    std::string name = Sample->GetName();
    TString NameTString = TString(name.c_str());

    // Clone the 1D histogram with a modified name
    TH1D* SampleHistogramPrior = static_cast<TH1D*>(Sample->get1DHist()->Clone(NameTString + "_Prior"));
    Sample->addData(SampleHistogramPrior);
    PDFs[i] = Sample;
  }
  return PDFs;
}
