#include "splines/SplineBase.h"

#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"

// *****************************************
SplineBase::SplineBase() {
// *****************************************

}


// *****************************************
SplineBase::~SplineBase() {
// *****************************************

}

// *****************************************
// Get the spline coefficients from the TSpline3 so that we can load ONLY these onto the GPU, not the whole TSpline3 object
// This loads up coefficients into two arrays: one x array and one yabcd array
// This should maximize our cache hits!
void SplineBase::getTF1Coeff(TF1_red* &spl, int &nPoints, float *& coeffs) {
// *****************************************
  // Initialise all arrays to 1.0
  for (int i = 0; i < _nTF1Coeff_; ++i) {
    coeffs[i] = 0.0;
  }

  // Get number of points in spline
  nPoints = spl->GetSize();

  if(nPoints > _nTF1Coeff_)
  {
    MACH3LOG_ERROR("Too big number of points for TF1");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  // TSpline3 can only take doubles, not floats
  // But our GPU is slow with doubles, so need to cast to float
  for (int i = 0; i < nPoints; i++) {
    coeffs[i] = float(spl->GetParameter(M3::int_t(i)));
  }
  // The structure is now coeffs  = {a,b,c,d,e}
}
