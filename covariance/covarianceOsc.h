#pragma once

// MaCh3 includes
#include "covariance/covarianceBase.h"

class covarianceOsc : public covarianceBase
{
 public:
  covarianceOsc(const char* name, const char *file);
  virtual ~covarianceOsc();
  double GetLikelihood();
  inline int checkBounds();
  double *getPropPars();
  void proposeStep();
  void setFlipDeltaM23(bool flip){flipdelM = flip;}
  void setFlipBeta(bool flip){flipBeta = flip;}
  void useReactorPrior(bool reactor){reactorPrior = reactor;};
  void setExtraBranches(TTree &tree);

  inline double GetPathLength() { return L;}
  inline double GetDensity() { return density;}
  inline bool GetPerformBetaStudy(){return PerformBetaStudy;}
  //KS: Print all useful information's after initialization
  void Print();
  //KS: Currently prob3++/probgp requires particular order so we need to check this is the case
  void CheckOrderOfParams();

 protected:
    double L, density;
    bool flipdelM;
    bool reactorPrior;
    bool flipBeta;
    double *oscpars1;
    bool PerformBetaStudy;

    int kDeltaCP;
    int kDeltaM23;
    int kSinTheta23;
    int kBeta;
};

