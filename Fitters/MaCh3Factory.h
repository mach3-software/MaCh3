#pragma once

// MaCh3 includes
#include "Fitters/FitterBase.h"
#include "Fitters/mcmc.h"
#include "Fitters/PSO.h"
#include "Fitters/LikelihoodFit.h"
#ifdef MaCh3_MINUIT2
#include "Fitters/MinuitFit.h"
#endif

#include "Parameters/ParameterHandlerGeneric.h"

/// @file MaCh3Factory.h
/// @brief Factory methods for MaCh3 software which streamline initialisation of different objects
/// @author Kamil Skwarczynski

/// @brief MaCh3 Factory initiates one of implemented fitting algorithms
/// @param fitMan pointer to Manager class
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

// ********************************************
/// @brief Factory function for creating a covariance class for systematic handling.
///
/// @param FitManager Pointer to the manager class that holds the configuration settings.
/// @param PreFix Prefix, for example Xsec, then code will look for XsecCovFile
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
std::unique_ptr<CovType> MaCh3CovarianceFactory(manager *FitManager, const std::string& PreFix){
// ********************************************
  // config for our matrix
  YAML::Node Settings = FitManager->raw()["General"]["Systematics"];
  auto CovMatrixName = Get<std::string>(Settings[std::string(PreFix) + "CovName"], __FILE__, __LINE__);
  MACH3LOG_INFO("Initialising {} matrix", CovMatrixName);

  // yaml files initialising out matrix
  auto CovMatrixFile = Get<std::vector<std::string>>(Settings[std::string(PreFix) + "CovFile"], __FILE__, __LINE__);

  // PCA threshold, -1 means no pca
  auto PCAThreshold = GetFromManager<int>(Settings[std::string(PreFix) + "PCAThreshold"], -1);
  // do we pca whole matrix or only submatrix
  auto PCAParamRegion = GetFromManager<std::vector<int>>(Settings[std::string(PreFix) + "PCAParams"], {-999, -999});

  auto CovObject = std::make_unique<CovType>(CovMatrixFile, CovMatrixName, PCAThreshold, PCAParamRegion[0], PCAParamRegion[1]);

  // Fill the parameter values with their nominal values
  // should _ALWAYS_ be done before overriding with fix or flat
  CovObject->SetParameters();

  auto FixParams = GetFromManager<std::vector<std::string>>(Settings[std::string(PreFix) + "Fix"], {});

  // Fixed CovObject parameters loop
  if (FixParams.size() == 1 && FixParams.at(0) == "All") {
    for (int j = 0; j < CovObject->GetNumParams(); j++) {
      CovObject->ToggleFixParameter(j);
    }
  } else {
    for (unsigned int j = 0; j < FixParams.size(); j++) {
      CovObject->ToggleFixParameter(FixParams.at(j));
    }
  }
  //Global step scale for matrix
  auto StepScale = Get<double>(Settings[std::string(PreFix) + "StepScale"], __FILE__, __LINE__);

  CovObject->SetStepScale(StepScale);

  // Adaptive MCMC stuff
  if (CheckNodeExists(FitManager->raw(), "AdaptionOptions")) {
    CovObject->InitialiseAdaption(FitManager->raw());
  }
  
  return CovObject;
}

// ********************************************
/// @brief Factory function for creating SampleHandler and initialisation with systematic.
///
/// @tparam SampleType The class type of the sample to create, e.g., `SampleHandlerTutorial`.
/// @param SampleConfig Path to sample config.
/// @param xsec A pointer to a ParameterHandlerGeneric object for cross-section systematic settings.
/// @return Vector of SampleType object, initialized and ready for use.
///
/// @note Example
/// ```cpp
/// auto mySamples = MaCh3SampleHandlerFactory<SampleHandlerTutorial>(SampleConfig, xsec);
/// ```
template <typename SampleType>
std::vector<SampleType*> MaCh3SampleHandlerFactory(const std::vector<std::string>& SampleConfig,
                                               ParameterHandlerGeneric* xsec) {
// ********************************************
  std::vector<SampleType*> Handlers(SampleConfig.size());
  for (size_t i = 0; i < SampleConfig.size(); ++i)
  {
    // Instantiate the sample using the specified class type
    SampleType* Sample = new SampleType(SampleConfig[i], xsec);
    Sample->Reweight();

    // Obtain sample name and create a TString version for histogram naming
    std::string name = Sample->GetTitle();
    TString NameTString = TString(name.c_str());

    // Clone the 1D histogram with a modified name
    if (Sample->GetNDim() == 1) {
      auto hist = static_cast<TH1D*>(Sample->GetMCHist(1));
      Sample->AddData(hist);
    } else {
      auto hist = static_cast<TH2D*>(Sample->GetMCHist(2));
      Sample->AddData(hist);
    }
    Handlers[i] = Sample;
  }
  return Handlers;
}
