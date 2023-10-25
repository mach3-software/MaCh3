#ifndef __LIKELIHOODFIT_H__
#define __LIKELIHOODFIT_H__

#include "FitterBase.h"

class LikelihoodFit : public FitterBase {
 public:
    LikelihoodFit(manager * const fitMan);
    ~LikelihoodFit();

    double CalcChi2(const double* x);
    int GetNPars(){return NPars;};

  protected:
    void PrepereFit();
    int NPars;
    int NParsPCA;
    bool fMirroring;
};

#endif
