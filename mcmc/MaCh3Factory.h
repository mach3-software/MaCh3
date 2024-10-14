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


/// @brief MaCh3 Factory initiates one of implemented fitting algorithms
/// @param fitMan pointer to Manager class
/// @param Samples vector of Sample PDF which will be used in the fits
/// @param Covariances vector of Systematic objects which will be used in the fits
///
/// @note Example YAML configuration:
/// @code
/// General:
///   FittingAlgorithm: ["MCMC"]
std::unique_ptr<FitterBase> MaCh3FitterFactory(manager *fitMan, std::vector<samplePDFBase>& Samples, std::vector<covarianceBase>& Covariances);


/// @brief Factory function for creating a covariance matrix for systematic handling.
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
covarianceXsec* MaCh3CovarianceFactory(manager *fitMan, const std::string& PreFix);


