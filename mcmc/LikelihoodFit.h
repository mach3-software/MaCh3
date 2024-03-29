#pragma once

#include "FitterBase.h"

class LikelihoodFit : public FitterBase {
 public:
    LikelihoodFit(manager * const fitMan);
    virtual ~LikelihoodFit();

    virtual double CalcChi2(const double* x);
    inline int GetNPars(){return NPars;};

    virtual void runMCMC() = 0;

  protected:
    void PrepareFit();
    int NPars;
    int NParsPCA;
    bool fMirroring;
};

