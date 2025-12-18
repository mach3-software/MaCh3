#pragma once

#include "Fitters/FitterBase.h"

/// @brief Implementation of base Likelihood Fit class, it is mostly responsible for likelihood calculation while actual use depend on classes which inherits.
class LikelihoodFit : public FitterBase {
 public:
    /// @brief Constructor
    LikelihoodFit(manager * const fitMan);
    /// @brief Destructor
    virtual ~LikelihoodFit();

    /// @brief Chi2 calculation over all included samples and syst objects
    virtual double CalcChi2(const double* x);
    //Expects parameters to be in proposal/PCA basis, transforms back before calling CalcChi2
    virtual double CalcChi2PC(const double* x);
    /// @brief Get total number of params, this sums over all covariance objects
    inline int GetNPars(){return NPars;};

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

