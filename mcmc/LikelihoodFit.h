#pragma once

#include "FitterBase.h"

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
    void PrepareFit();
    int NPars;
    int NParsPCA;
    bool fMirroring;
};

