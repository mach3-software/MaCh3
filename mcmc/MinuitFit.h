#ifndef __MINUITFIT_H__
#define __MINUITFIT_H__

#include "LikelihoodFit.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


class MinuitFit : public LikelihoodFit {
 public:
  MinuitFit(manager * const fitMan);
  ~MinuitFit();

  void runFit();

  private:
    ROOT::Math::Minimizer* minuit;
};


#endif
