#pragma once

// MaCh3  includes
#include "covariance/covarianceXsec.h"
#include "samplePDF/Structs.h"


// *******************
/// CW: Add a struct to hold info about the splinified xsec parameters and help with FindSplineSegment
struct FastSplineInfo {
// *******************
  // Number of points in spline
  _int_ nPts;

  // Array of the knots positions
  _float_ *xPts;

  // Array of what segment of spline we're currently interested in
  // Gets updated once per MCMC iteration
  _int_ CurrSegment;

  // Array of the knots positions
  const double* splineParsPointer;
};

// ********************************************
/// CW: Generic xsec class
/// Can use TF1 or TSpline3 or TSpline5 here, tjoho
template <class T>
class XSecStruct{
// ********************************************
public:
  /// @brief CW: The light constructor
  XSecStruct(_int_ NumberOfSplines) {
    nParams = NumberOfSplines;
    Func.reserve(nParams);
    for (int i = 0; i < nParams; ++i) {
      Func[i] = NULL;
    }
  }

  /// @brief CW: The empty constructor
  XSecStruct() {
    nParams = 0;
    Func = NULL;
  };

  /// @brief CW: The light destructor
  ~XSecStruct() {
    for (int i = 0; i < nParams; ++i) {
      if (Func[i]) delete Func[i];
    }
  }

  /// @brief CW: Get number of splines
  inline _int_ GetNumberOfParams() { return nParams; }

  /// @brief CW: The Printer
  inline void Print() {
    std::cout << "    Splines: " << std::endl;
    for (int i = 0; i < nParams; ++i) {
      if (!Func[i]) continue;
      std::cout << "    " << std::left << std::setw(25) << Func[i]->GetName() << std::endl;
    }
  }

  /// @brief CW: Set the number of splines for this event
  inline void SetSplineNumber(const _int_ NumberOfSplines) {
    nParams = NumberOfSplines;
    Func = new T[nParams];
  }

  /// @brief CW: Get the function for the nth spline
  inline T GetFunc(const _int_ nSpline) { return Func[nSpline]; }
  /// @brief CW: Set the function for the nth spline
  inline void SetFunc(const _int_ nSpline, T Function) { Func[nSpline] = Function; }
  /// @brief CW: Eval the current variation
  inline double Eval(const _int_ nSpline, const _float_ variation) {
    // Some will be NULL, check this
    if (Func[nSpline]) {
      return Func[nSpline]->Eval(variation);
    } else {
      return 1.0;
    }
  }
private:
  /// Number of parameters
  _int_ nParams;
  /// The function
  T* Func;
};


// ***************************************************************************
//EM: Apply capping to knot weight for specified spline parameter. param graph needs to have been set in xsecgraph array first
inline void ApplyKnotWeightCap(TGraph* xsecgraph, int splineParsIndex, covarianceXsec* XsecCov) {
// ***************************************************************************
  if(xsecgraph == NULL){
    MACH3LOG_ERROR("ERROR: hmmm looks like you're trying to apply capping for spline parameter {} but it hasn't been set in xsecgraph yet", XsecCov->GetParFancyName(splineParsIndex));
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // EM: cap the weights of the knots if specified in the config
  if(
    (XsecCov->GetParSplineKnotUpperBound(splineParsIndex) != 9999)
    || (XsecCov->GetParSplineKnotLowerBound(splineParsIndex) != -9999))
  {
    for(int knotId = 0; knotId < xsecgraph->GetN(); knotId++){
      double x,y;

      // EM: get the x and y of the point. Also double check that the requested point was valid just to be super safe
      if( xsecgraph->GetPoint(knotId, x, y) == -1) {
        MACH3LOG_ERROR("Invalid knot requested: {}", knotId);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }

      y = std::min(y, XsecCov->GetParSplineKnotUpperBound(splineParsIndex));
      y = std::max(y, XsecCov->GetParSplineKnotLowerBound(splineParsIndex));

      xsecgraph->SetPoint(knotId, x, y);
    }

    // EM: check if our cap made the spline flat, if so we set to 1 knot to avoid problems later on
    bool isFlat = true;

    for(int knotId = 0; knotId < xsecgraph->GetN(); knotId++){
      double x,y;

      // EM: get the x and y of the point. Also double check that the requested point was valid just to be super safe
      if( xsecgraph->GetPoint(knotId, x, y) == -1) {
        MACH3LOG_ERROR("Invalid knot requested: {}", knotId);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }

      if(std::abs(y - 1.0) > 1e-5) isFlat = false;
    }

    if(isFlat){
      xsecgraph->Set(1);
      xsecgraph->SetPoint(0, 0.0, 1.0);
    }
  }
}


// TODO KS: Consider moving here FastSplineInfo, and TSpline3_red etc
