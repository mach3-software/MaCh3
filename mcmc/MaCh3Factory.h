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
covarianceXsec* MaCh3CovarianceFactory(manager *fitMan, const std::string& PreFix);

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
                                               covarianceOsc* osc = nullptr)
{
  std::vector<SampleType*> PDFs(SampleConfig.size());
  for (size_t i = 0; i < SampleConfig.size(); ++i)
  {
    // Instantiate the sample using the specified class type
    SampleType* Sample = new SampleType(SampleConfig[i], xsec);
    Sample->SetXsecCov(xsec);
    if (osc != nullptr) Sample->SetOscCov(osc);
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
