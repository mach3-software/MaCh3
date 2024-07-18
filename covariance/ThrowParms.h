#pragma once

#include <iostream>
#include <assert.h>
#include <algorithm>
#include <math.h>
#include <vector>

#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TVectorT.h>
#include <TDecompChol.h>

/// @brief No idea really...
/// @warning this is deprecated don't use it
class ThrowParms {
 private:
  TVectorD    *pvals;
  TMatrixDSym *covar; 
  TMatrixD    *chel_dec;
  int         npars;
  TRandom3    rand; 

  void CheloskyDecomp(TMatrixD &chel_mat);
  void StdNormRand(double *z);

 public:
  ThrowParms(TVectorD &parms, TMatrixDSym &covm);
  ~ThrowParms();

  void SetSeed (int seed = 0) {rand.SetSeed(seed);};
  int GetSize() {return npars;};
  void ThrowSet(std::vector<double> &parms);

};
