#pragma once

#include "mcmc/FitterBase.h"

/// @brief Implementation of base Likelihood Fit class, it is mostly responsible for likelihood calculation while actual use depend on classes which inherits.
class LikelihoodFit : public FitterBase {
 public:
    /// @brief Constructor
    LikelihoodFit(manager * const fitMan);
    /// @brief Destructor
    virtual ~LikelihoodFit();

    /// @brief Chi2 calculation over all included samples and syst objects
    virtual double CalcChi2(const double* x);
    /// @brief Get total number of params, this sums over all covariance objects
    inline int GetNPars(){return NPars;};

    /// @brief Implementation of fitting algorithm
    virtual void runMCMC() = 0;

    /// @brief Get name of class
    virtual inline std::string GetName()const {return "LikelihoodFit";};
  protected:
    /// @brief prepare output and perform sanity checks
    void PrepareFit();
    /// Number of all parameters from all covariances
    int NPars;
    /// Number of all parameters from all covariances in PCA base
    int NParsPCA;
    /// Flag telling if mirroring is used or not
    bool fMirroring;
};

