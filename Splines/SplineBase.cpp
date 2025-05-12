#include "Splines/SplineBase.h"

#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"

// *****************************************
SplineBase::SplineBase() {
// *****************************************
  nParams = 0;
  SplineSegments = nullptr;
  ParamValues = nullptr;
}


// *****************************************
SplineBase::~SplineBase() {
// *****************************************
}

// *************************
// CW: Only need to do the binary search once per parameter, not once per event!
// Takes down the number of binary searches from 1.2M to 17, ha!
// For ROOT version see root/hist/hist/src/TSpline3.cxx TSpline3::FindX(double)
void SplineBase::FindSplineSegment() {
// *************************
  // Loop over the splines
  //KS: Tried multithreading here with 48 splines and it is faster with one thread, maybe in future multithreading will be worth revisiting
  for (M3::int_t i = 0; i < nParams; ++i)
  {
    const M3::int_t nPoints = SplineInfoArray[i].nPts;
    const std::vector<M3::float_t>& xArray = SplineInfoArray[i].xPts;

    // Get the variation for this reconfigure for the ith parameter
    const float xvar = float(*SplineInfoArray[i].splineParsPointer);
    ParamValues[i] = xvar;

    // EM: if we have a parameter that has no response for any event (i.e. all splines have just one knot), then skip it and avoid a seg fault here
    //     In principle, such parameters shouldn't really be included in the first place, but with new det syst splines this
    //     could happen if say you were just running on one FHC run, then all RHC parameters would be flat and the code below would break.
    if(xArray.size() == 0) continue;

    // The segment we're interested in (klow in ROOT code)
    M3::int_t segment = 0;
    M3::int_t kHigh = nPoints-1;
    //KS: We expect new segment is very close to previous
    const M3::int_t PreviousSegment = SplineInfoArray[i].CurrSegment;

    // If the variation is below the lowest saved spline point
    if (xvar <= xArray[0]) {
      segment = 0;
      // If the variation is above the highest saved spline point
    } else if (xvar >= xArray[nPoints-1]) {
      //CW: Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
      segment = kHigh;
      //KS: It is quite probable the new segment is same as in previous step so try to avoid binary search, first we have to check if it is in bounds to avoid seg fault
    } else if( xArray[PreviousSegment+1] > xvar && xvar >= xArray[PreviousSegment] ) {
      segment = PreviousSegment;
      // If the variation is between the maximum and minimum, perform a binary search
    } else {
      // The top point we've got
      M3::int_t kHalf = 0;
      // While there is still a difference in the points (we haven't yet found the segment)
      // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment
      while (kHigh - segment > 1) {
        // Increment the half-step
        kHalf = M3::int_t((segment + kHigh)/2);
        // If our variation is above the kHalf, set the segment to kHalf
        if (xvar > xArray[kHalf]) {
          segment = kHalf;
          // Else move kHigh down
        } else {
          kHigh = kHalf;
        }
      } // End the while: we've now done our binary search
    } // End the else: we've now found our point

    if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;

    //CW: This way we avoid doing 1.2M+ binary searches on the GPU
    // and literally just multiply lots of numbers together on the GPU without any algorithm
    // Update the values and which segment it belongs to
    SplineInfoArray[i].CurrSegment = segment;
    SplineSegments[i] = short(SplineInfoArray[i].CurrSegment);

    #ifdef DEBUG
    if (SplineInfoArray[i].xPts[segment] > xvar && segment != 0) {
      MACH3LOG_ERROR("Found a segment which is _ABOVE_ the variation!");
      MACH3LOG_ERROR("IT SHOULD ALWAYS BE BELOW! (except when segment 0)");
      MACH3LOG_ERROR("Spline: {}", i);
      MACH3LOG_ERROR("Found segment = {}", segment);
      MACH3LOG_ERROR("Doing variation = {}", xvar);
      MACH3LOG_ERROR("x in spline = {}", SplineInfoArray[i].xPts[segment]);
      for (M3::int_t j = 0; j < SplineInfoArray[j].nPts; ++j) {
        MACH3LOG_ERROR("    {} = {}", j, SplineInfoArray[i].xPts[j]);
      }
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    #endif
  } //end loop over params
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
