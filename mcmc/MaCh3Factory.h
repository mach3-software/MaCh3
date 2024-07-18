#pragma once

// MaCh3 includes
#include "mcmc/FitterBase.h"
#include "mcmc/mcmc.h"
#include "mcmc/PSO.h"
#include "mcmc/LikelihoodFit.h"
#ifdef MaCh3_MINUIT2
#include "mcmc/MinuitFit.h"
#endif

/// @brief MaCh3 Factory initiates one of implemented fitting algorithms
/// @param fitMan pointer to Manager class
/// @param Samples vector of Sample PDF which will be used in the fits
/// @param Covariances vector of Systematic objects which will be used in the fits
FitterBase* MaCh3FitterFactory(manager *fitMan, std::vector<samplePDFBase>& Samples, std::vector<covarianceBase>& Covariances);
