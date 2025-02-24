#pragma once

// MaCh3  includes
#include "ParameterHandler/ParameterHandlerGeneric.h"
#include "SampleHandler/Structs.h"

#include <cmath>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wnull-dereference"

/// @file SplineStructs.h
/// @brief Contains structures and helper functions for handling spline representations of systematic parameters in the MaCh3.
/// @author Clarence Wret
/// @author Kamil Skwarczynski

// *******************
/// @brief CW: Add a struct to hold info about the splinified xsec parameters and help with FindSplineSegment
struct FastSplineInfo {
// *******************
  /// @brief Constructor
  FastSplineInfo() {
    nPts = -999;
    CurrSegment = 0;
    splineParsPointer = nullptr;
  }

  /// @brief Destructor
  virtual ~FastSplineInfo() {

  }

  /// Number of points in spline
  M3::int_t nPts;

  /// Array of the knots positions
  std::vector<M3::float_t> xPts;

  /// Array of what segment of spline we're currently interested in. Gets updated once per MCMC iteration
  M3::int_t CurrSegment;

  /// Array of the knots positions
  const double* splineParsPointer;
};


// ***************************************************************************
/// @brief EM: Apply capping to knot weight for specified spline parameter. param graph needs to have been set in xsecgraph array first
inline void ApplyKnotWeightCap(TGraph* xsecgraph, int splineParsIndex, ParameterHandlerGeneric* XsecCov) {
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

// ************************
/// @brief KS: A reduced ResponseFunction Generic function used for evaluating weight
class TResponseFunction_red {
// ************************
public:
  /// @brief Empty constructor
  TResponseFunction_red() { }
  /// @brief Empty destructor
  virtual ~TResponseFunction_red() { }
  /// @brief Evaluate a variation
  virtual double Eval(const double var)=0;
  /// @brief KS: Printer
  virtual void Print()=0;
  /// @brief DL: Get number of points
  virtual M3::int_t GetNp()=0;
};

// ************************
/// @brief CW: A reduced TF1 class only. Only saves parameters for each TF1 and how many parameters each parameter set has
class TF1_red: public TResponseFunction_red {
// ************************
public:
  /// @brief Empty constructor
  TF1_red() : TResponseFunction_red() {
    length = 0;
    Par = NULL;
  }

  /// @brief Empty destructor
  virtual ~TF1_red() {
    if (Par != NULL) {
      delete[] Par;
      Par = NULL;
    }
  }

  /// @brief The useful constructor with deep copy
  TF1_red(M3::int_t nSize, M3::float_t* Array) : TResponseFunction_red() {
    length = nSize;
    for (int i = 0; i < length; ++i) {
      Par[i] = Array[i];
    }
  }

  /// @brief The TF1 constructor with deep copy
  TF1_red(TF1* &Function) : TResponseFunction_red() {
    Par = NULL;
    SetFunc(Function);
  }

  /// @brief Set the function
  inline void SetFunc(TF1* &Func) {
    length = M3::int_t(Func->GetNpar());
    if (Par != NULL) delete[] Par;
    Par = new M3::float_t[length];
    for (int i = 0; i < length; ++i) {
      Par[i] = M3::float_t(Func->GetParameter(i));
    }
    delete Func;
    Func = NULL;
  }

  /// @brief Evaluate a variation
  inline double Eval(const double var) override {
    return Par[1]+Par[0]*var;

    /* FIXME in future we might introduce more TF1
    //If we have 5 parameters we're using a fifth order polynomial
    if (Type == kFifthOrderPolynomial) {
      return 1+Par[0]*var+Par[1]*var*var+Par[2]*var*var*var+Par[3]*var*var*var*var+Par[4]*var*var*var*var*var;
    } else if (Type == kTwoLinears) {
      return (var <= 0)*(Par[2]+Par[0]*var)+(var > 0)*(Par[2]+Par[1]*var);
    } else if (Type == kLinear) {
      return (Par[1]+Par[0]*var);
    } else if (Type == kPseudoHeaviside) {
      return (var <= 0)*(1+Par[0]*var) + (1 >= var)*(var > 0)*(1+Par[1]*var) + (var > 1)*(Par[3]+Par[2]*var);
    }else {
      MACH3LOG_ERROR("    Class only knows about 5th order polynomial, two superposed linear functions, linear function, or pseudo Heaviside.");
      MACH3LOG_ERROR("    You have tried something else than this, which remains unimplemented.");
      MACH3LOG_ERROR("{}: {}", __FILE__, __LINE__);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    */
  }

  /// @brief Set a parameter to a value
  inline void SetParameter(M3::int_t Parameter, M3::float_t Value) {
    Par[Parameter] = Value;
  }

  /// @brief Get a parameter value
  double GetParameter(M3::int_t Parameter) {
    if (Parameter > length) {
      MACH3LOG_ERROR("You requested parameter number {} but length is {} parameters", Parameter, length);
      throw MaCh3Exception(__FILE__ , __LINE__ );
      return -999.999;
    }
    return Par[Parameter];
  }

  /// @brief Set the size
  inline void SetSize(M3::int_t nSpline) {
    length = nSpline;
    Par = new M3::float_t[length];
  }
  /// @brief Get the size
  inline int GetSize() { return length; }
  /// @brief Print detailed info
  inline void Print() override {
    MACH3LOG_INFO("Printing TF1_red:");
    MACH3LOG_INFO("  Length  = {}", length);
    for (int i = 0; i < length; i++) {
      MACH3LOG_INFO("  Coeff {}  = {}", i, Par[i]);
    }
  }

  /// @brief KS: Make a TF1 from the reduced TF1
  inline TF1* ConstructTF1(const std::string& function, const int xmin, const int xmax) {
    TF1 *func = new TF1("TF1", function.c_str(), xmin, xmax);
    for(int i = 0; i < length; ++i) {
      func->SetParameter(i, Par[i]);
    }
    return func;
  }

  /// @brief DL: Get number of points
  inline M3::int_t GetNp() override { return length; }

private:
  /// The parameters
  M3::float_t* Par;
  M3::int_t length;
};

// ************************
/// CW: Reduced TSpline3 class
class TSpline3_red: public TResponseFunction_red {
// ************************
public:
  /// @brief Empty constructor
  TSpline3_red() : TResponseFunction_red() {
    nPoints = 0;
    Par = NULL;
    XPos = NULL;
    YResp = NULL;
  }

  /// @brief The constructor that takes a TSpline3 pointer and copies in to memory
  TSpline3_red(TSpline3* &spline, SplineInterpolation InterPolation = kTSpline3) : TResponseFunction_red() {
    Par = NULL;
    XPos = NULL;
    YResp = NULL;
    SetFunc(spline, InterPolation);
  }

  /// @brief constructor taking parameters
  TSpline3_red(M3::float_t *X, M3::float_t *Y, M3::int_t N, M3::float_t **P) : TResponseFunction_red() {
    nPoints = N;
    // Save the parameters for each knot
    Par = new M3::float_t*[nPoints];
    // Save the positions of the knots
    XPos = new M3::float_t[nPoints];
    // Save the y response at each knot
    YResp = new M3::float_t[nPoints];
    for(int j = 0; j < N; ++j){
      Par[j] = new M3::float_t[3];
      Par[j][0] = P[j][0];
      Par[j][1] = P[j][1];
      Par[j][2] = P[j][2];
      XPos[j] = X[j];
      YResp[j] = Y[j];

      if((Par[j][0] == -999) | (Par[j][1] ==-999) | (Par[j][2] ==-999) | (XPos[j] ==-999) | (YResp[j] ==-999)){
        MACH3LOG_ERROR("******************* Bad parameter values when constructing TSpline3_red *********************");
        MACH3LOG_ERROR("Passed values (i, x, y, b, c, d): {}, {}, {}, {}, {}, {}", j, X[j], Y[j], P[j][0], P[j][1], P[j][2]);
        MACH3LOG_ERROR("Set values (i, x, y, b, c, d): {}, {}, {}, {}, {}, {}", j, XPos[j], YResp[j], Par[j][0], Par[j][1], Par[j][2]);
        MACH3LOG_ERROR("*********************************************************************************************");
      }
    }
  }
  /// @brief Set the function
  inline void SetFunc(TSpline3* &spline, SplineInterpolation InterPolation = kTSpline3) {
    nPoints = M3::int_t(spline->GetNp());
    if (Par != NULL) {
      for (int i = 0; i < nPoints; ++i) {
        delete[] Par[i];
        Par[i] = NULL;
      }
      delete[] Par;
      Par = NULL;
    }
    if (XPos != NULL) delete[] XPos;
    if (YResp != NULL) delete[] YResp;
    // Save the parameters for each knot
    Par = new M3::float_t*[nPoints];
    // Save the positions of the knots
    XPos = new M3::float_t[nPoints];
    // Save the y response at each knot
    YResp = new M3::float_t[nPoints];

    //KS: Default TSpline3 ROOT implementation
    if(InterPolation == kTSpline3)
    {
      for (int i = 0; i < nPoints; ++i) {
        // 3 is the size of the TSpline3 coefficients
        Par[i] = new M3::float_t[3];
        double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
        spline->GetCoeff(i, x, y, b, c, d);
        XPos[i]   = M3::float_t(x);
        YResp[i]  = M3::float_t(y);
        Par[i][0] = M3::float_t(b);
        Par[i][1] = M3::float_t(c);
        Par[i][2] = M3::float_t(d);
      }
    }
    //CW: Reduce to use linear spline interpolation for certain parameters
    // Not the most elegant way: use TSpline3 object but set coefficients to zero and recalculate spline points; the smart way (but more human intensive) would be to save memory here and simply not store the zeros at all
    // Get which parameters should be linear from the fit manager
    // Convert the spline number to global xsec parameter
    // Loop over the splines points
    // KS: kLinearFunc should be used with TF1, this is just as safety
    else if(InterPolation == kLinear || InterPolation == kLinearFunc)
    {
      for (int k = 0; k < nPoints; ++k) {
        // 3 is the size of the TSpline3 coefficients
        Par[k] = new M3::float_t[3];
        Double_t x1, y1, b1, c1, d1, x2, y2, b2, c2, d2 = 0;
        spline->GetCoeff(k, x1, y1, b1, c1, d1);
        spline->GetCoeff(k+1, x2, y2, b2, c2, d2);
        double tempb = (y2-y1)/(x2-x1);

        XPos[k]   = M3::float_t(x1);
        YResp[k]  = M3::float_t(y1);
        Par[k][0] = M3::float_t(tempb);
        Par[k][1] = M3::float_t(0);
        Par[k][2] = M3::float_t(0);
      }
    }
    //EM: Akima spline is similar to regular cubic spline but is allowed to be discontinuous in 2nd derivative and coefficients in any segment
    // only depend on th 2 nearest points on either side
    else if(InterPolation == kAkima)
    {
      // get the knot values for the spline
      for (int i = 0; i < nPoints; ++i) {
        // 3 is the size of the TSpline3 coefficients
        Par[i] = new M3::float_t[3];

        double x = -999.99, y = -999.99;
        spline->GetKnot(i, x, y);

        XPos[i]   = M3::float_t(x);
        YResp[i]  = M3::float_t(y);
      }

      M3::float_t* mvals = new M3::float_t[nPoints + 3];
      M3::float_t* svals = new M3::float_t[nPoints + 1];

      for (int i = -2; i <= nPoints; ++i) {
        // if segment is first or last or 2nd to first or last, needs to be dealt with slightly differently;
        // need to estimate the values for additinal points which would lie outside of the spline
        if(i ==-2){
          mvals[i+2] = M3::float_t(3.0 * (YResp[1] - YResp[0]) / (XPos[1] - XPos[0]) - 2.0*(YResp[2] - YResp[1]) / (XPos[2] - XPos[1]));
        }
        else if(i==-1){
          mvals[i+2] = M3::float_t(2.0 * (YResp[1] - YResp[0]) / (XPos[1] - XPos[0]) - (YResp[2] - YResp[1]) / (XPos[2] - XPos[1]));
        }
        else if(i==nPoints){
          mvals[i+2] = M3::float_t(3.0 * (YResp[nPoints-1] - YResp[nPoints-2]) / (XPos[nPoints-1] - XPos[nPoints-2]) - 2.0*(YResp[nPoints-2] - YResp[nPoints-3]) / (XPos[nPoints-2] - XPos[nPoints-3]));
        }
        else if(i == nPoints - 1){
          mvals[i+2] = M3::float_t(2.0 * (YResp[nPoints-1] - YResp[nPoints-2]) / (XPos[nPoints-1] - XPos[nPoints-2]) - (YResp[nPoints-2] - YResp[nPoints-3]) / (XPos[nPoints-2] - XPos[nPoints-3]));
        }
        //standard internal segment
        else{
          mvals[i+2] = (YResp[i+1] - YResp[i])/ (XPos[i+1] - XPos[i]);
        }
      }

      for(int i =2; i<=nPoints+2; i++){
        if (std::abs(mvals[i+1] - mvals[i]) + std::abs(mvals[i-1] - mvals[i-2]) != 0.0){
          svals[i-2] = (std::abs(mvals[i+1] - mvals[i]) * mvals[i-1] + std::abs(mvals[i-1] - mvals[i-2]) *mvals[i]) / (std::abs(mvals[i+1] - mvals[i]) + std::abs(mvals[i-1] - mvals[i-2]));
        }
        else{svals[i-2] = mvals[i];}
      }

      // calculate the coefficients for the spline
      for(int i = 0; i <nPoints; i++){
        M3::float_t b, c, d = M3::float_t(-999.999);

        b = svals[i];
        c = M3::float_t(3.0* (YResp[i+1] - YResp[i]) / (XPos[i+1] - XPos[i]) -2.0 *svals[i] - svals[i +1]) /(XPos[i+1] - XPos[i]);
        d = M3::float_t((svals[i + 1] +svals[i]) - 2.0*(YResp[i+1] - YResp[i]) / (XPos[i+1] - XPos[i])) / ((XPos[i+1] - XPos[i]) * (XPos[i+1] - XPos[i]));

        Par[i][0] = b;
        Par[i][1] = c;
        Par[i][2] = d;
      }

      // check the input spline for linear segments, if there are any then overwrite the calculated coefficients
      // this will pretty much only ever be the case if they are set to be linear in samplePDFND i.e. the user wants it to be linear
      for(int i = 0; i <nPoints-1; i++){
        double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
        spline->GetCoeff(i, x, y, b, c, d);

        if((c == 0.0 && d == 0.0)){
          Par[i][0] = M3::float_t(b);
          Par[i][1] = M3::float_t(0.0);
          Par[i][2] = M3::float_t(0.0);
        }
      }
      delete[] mvals;
      delete[] svals;
    }
    //EM: Monotone spline is similar to regular cubic spline but enforce the condition that the interpolated value at any point
    // must be between its two nearest knots, DOES NOT make the entire spline monotonic, only the segments
    else if(InterPolation == kMonotonic)
    {
      // values of the secants at each point (for calculating monotone spline)
      M3::float_t * Secants = new M3::float_t[nPoints -1];
      // values of the tangens at each point (for calculating monotone spline)
      M3::float_t *  Tangents = new M3::float_t[nPoints];

      // get the knot values for the spline
      for (int i = 0; i < nPoints; ++i) {
        // 3 is the size of the TSpline3 coefficients
        Par[i] = new M3::float_t[3];

        double x = -999.99, y = -999.99;
        spline->GetKnot(i, x, y);

        XPos[i]   = M3::float_t(x);
        YResp[i]  = M3::float_t(y);

        Tangents[i] = 0.0;
      }

      // deal with the case of two points (just do linear interpolation between them)
      if (nPoints == 2){
        Par[0][0] = (YResp[1] - YResp[0]) / (XPos[1] - XPos[0]);
        Par[0][1] = 0.0;
        Par[0][2] = 0.0;
        // extra "virtual" segment at end to make Par array shape fit with knot arrays shapes
        Par[1][1] = 0.0;
        Par[1][2] = 0.0;

        return;
      } // if nPoints !=2 do full monotonic spline treatment:
      else
      {
        // first pass over knots to calculate the secants
        for (int i = 0; i < nPoints-1; ++i) {
          Secants[i] = (YResp[i+1] - YResp[i]) / (XPos[i+1] - XPos[i]);
          MACH3LOG_TRACE("Secant {}: {}", i, Secants[i]);
        }

        Tangents[0] = Secants[0];
        Tangents[nPoints-1] = Secants[nPoints -2];

        M3::float_t alpha;
        M3::float_t beta;

        // second pass over knots to calculate tangents
        for (int i = 1; i < nPoints-1; ++i) {
          if ((Secants[i-1] >= 0.0 && Secants[i] >= 0.0) | (Secants[i-1] < 0.0 && Secants[i] < 0.0)){ //check for same sign
            Tangents[i] = M3::float_t((Secants[i-1] + Secants[i]) /2.0);
          }
        }

        // third pass over knots to rescale tangents
        for (int i = 0; i < nPoints-1; ++i) {
          if (Secants[i] == 0.0){
            Tangents[i] = 0.0;
            Tangents[i+1] = 0.0;
          }

          else{
            alpha = Tangents[i]  / Secants[i];
            beta = Tangents[i+1] / Secants[i];

            if (alpha <0.0){
              Tangents[i] = 0.0;
            }
            if (beta < 0.0){
              Tangents[i+1] = 0.0;
            }

            if (alpha * alpha + beta * beta >9.0){
              M3::float_t tau = M3::float_t(3.0 / std::sqrt(alpha * alpha + beta * beta));
              Tangents[i]   = tau * alpha * Secants[i];
              Tangents[i+1] = tau * beta  * Secants[i];
            }
          }
        } // finished rescaling tangents
        // fourth pass over knots to calculate the coefficients for the spline
        M3::float_t dx;
        for(int i = 0; i <nPoints-1; i++){
          M3::float_t b, c, d = M3::float_t(-999.999);
          dx = XPos[i+1] - XPos[i];

          b = Tangents[i] * dx;
          c = M3::float_t(3.0* (YResp[i+1] - YResp[i]) -2.0 *dx * Tangents[i] - dx * Tangents[i +1]);
          d = M3::float_t(2.0* (YResp[i] - YResp[i+1]) + dx * (Tangents[i] + Tangents[i+1]));

          Par[i][0] = b /  dx;
          Par[i][1] = c / (dx * dx);
          Par[i][2] = d / (dx * dx * dx);

          if((Par[i][0] == -999) | (Par[i][1] == -999) | (Par[i][2] ==-999) | (Par[i][0] == -999.999) | (Par[i][1] == -999.999) | (Par[i][2] ==-999.999)){
            MACH3LOG_INFO("Bad spline parameters for segment {}: (b, c, d) = {}, {}, {}. This will cause problems with GPU.",
                          i, Par[i][0], Par[i][1], Par[i][2]);
          }
          MACH3LOG_TRACE("b: {}", b);
          MACH3LOG_TRACE("dx: {}, x_0: {}, x_1: {}", dx, XPos[i], XPos[i+1]);
          MACH3LOG_TRACE("    y_0: {}, y_1: {}", YResp[i], YResp[i+1]);
        }

        // include params for final "segment" outside of the spline so that par array fits with x and y arrays,
        // should never actually get used but if not set then the GPU code gets very angry
        Par[nPoints-1][0] = 0.0;
        Par[nPoints-1][1] = 0.0;
        Par[nPoints-1][2] = 0.0;

        // check the input spline for linear segments, if there are any then overwrite the calculated coefficients
        // this will pretty much only ever be the case if they are set to be linear in samplePDFND i.e. the user wants it to be linear
        for(int i = 0; i <nPoints-1; i++){
          double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
          spline->GetCoeff(i, x, y, b, c, d);

          if((c == 0.0 && d == 0.0)){
            Par[i][0] = M3::float_t(b);
            Par[i][1] = 0.0;
            Par[i][2] = 0.0;
          }
        }
        delete[] Secants;
        delete[] Tangents;
      } // end of if(nPoints !=2)
    }
    else
    {
      MACH3LOG_ERROR("Unsupported interpolation type {}", static_cast<int>(InterPolation));
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    delete spline;
    spline = NULL;
  }

  /// @brief Empty destructor
  virtual ~TSpline3_red() {
    if(Par != NULL) {
      for (int i = 0; i < nPoints; ++i) {
        if (Par[i] != NULL) {
          delete[] Par[i];
        }
      }
      delete[] Par;
    }
    if(XPos != NULL) delete[] XPos;
    if(YResp != NULL) delete[] YResp;
    Par = NULL;
    XPos = YResp = NULL;
  }

  /// @brief Find the segment relevant to this variation in x
  /// See root/hist/hist/src/TSpline3::FindX(double) or samplePDFND....::FindSplineSegment
  inline int FindX(double x) {
    // The segment we're interested in (klow in ROOT code)
    int segment = 0;
    int kHigh = nPoints-1;
    // If the variation is below the lowest saved spline point
    if (x <= XPos[0]){
      segment = 0;
      // If the variation is above the highest saved spline point
    } else if (x >= XPos[nPoints-1]) {
      // Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
      segment = kHigh;
      // If the variation is between the maximum and minimum, perform a binary search
    } else {
      // The top point we've got
      int kHalf = 0;
      // While there is still a difference in the points (we haven't yet found the segment)
      // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment
      while (kHigh - segment > 1) {
        // Increment the half-step
        kHalf = (segment + kHigh)/2;
        // If our variation is above the kHalf, set the segment to kHalf
        if (x > XPos[kHalf]) {
          segment = kHalf;
          // Else move kHigh down
        } else {
          kHigh = kHalf;
        }
      } // End the while: we've now done our binary search
    } // End the else: we've now found our point
    if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;
    return segment;
  }

  /// @brief CW: Evaluate the weight from a variation
  inline double Eval(double var) override {
    // Get the segment for this variation
    int segment = FindX(var);
    // The get the coefficients for this variation
    M3::float_t x = M3::float_t(-999.99), y = M3::float_t(-999.99), b = M3::float_t(-999.99), c = M3::float_t(-999.99), d = M3::float_t(-999.99);
    GetCoeff(segment, x, y, b, c, d);
    double dx = var - x;
    // Evaluate the third order polynomial
    double weight = y+dx*(b+dx*(c+d*dx));
    return weight;
  }

  /// @brief CW: Get the number of points
  inline M3::int_t GetNp() override { return nPoints; }
  // Get the ith knot's x and y position
  inline void GetKnot(int i, M3::float_t &xtmp, M3::float_t &ytmp) {
    xtmp = XPos[i];
    ytmp = YResp[i];
  }

  /// @brief CW: Get the coefficient of a given segment
  inline void GetCoeff(int segment, M3::float_t &x, M3::float_t &y, M3::float_t &b, M3::float_t &c, M3::float_t &d) {
    b = Par[segment][0];
    c = Par[segment][1];
    d = Par[segment][2];
    x = XPos[segment];
    y = YResp[segment];
  }

  /// @brief CW: Make a TSpline3 from the reduced splines
  inline TSpline3* ConstructTSpline3() {
    // KS: Sadly ROOT only accepts double...
    #ifdef _LOW_MEMORY_STRUCTS_
    std::vector<Double_t> xPosDoubles(nPoints);
    std::vector<Double_t> yPosDoubles(nPoints);
    for (Int_t i = 0; i < nPoints; ++i) {
      xPosDoubles[i] = static_cast<Double_t>(XPos[i]); // Convert float to double
      yPosDoubles[i] = static_cast<Double_t>(YResp[i]); // Convert float to double
    }
    TSpline3 *spline = new TSpline3("Spline", xPosDoubles.data(), yPosDoubles.data(), static_cast<int>(nPoints));
    #else
    TSpline3 *spline = new TSpline3("Spline", XPos, YResp, nPoints);
    #endif
    return spline;
  }

  /// @brief Print detailed info
  inline void Print() override {
    MACH3LOG_INFO("Printing TSpline_red:");
    MACH3LOG_INFO(" Nknots = {}", nPoints);
    for (int i = 0; i < nPoints; ++i) {
      MACH3LOG_INFO("  i = {} x = {} y = {} b = {} c = {} d = {}",
                    i, XPos[i], YResp[i], Par[i][0], Par[i][1], Par[i][2]);
    }
  }

  protected: //changed to protected from private so can be accessed by derived classes
    /// Number of points/knot in TSpline3
    M3::int_t nPoints;
    /// Always uses a third order polynomial, so hard-code the number of coefficients in implementation
    M3::float_t **Par;
    /// Positions of each x for each knot
    M3::float_t *XPos;
    /// y-value for each knot
    M3::float_t *YResp;
};

// *****************************************
/// @brief CW: Helper function used in the constructor, tests to see if the spline is flat
/// @param spl pointer to TSpline3_red that will be checked
inline bool isFlat(TSpline3_red* &spl) {
// *****************************************
  int Np = spl->GetNp();
  M3::float_t x, y, b, c, d;
  // Go through spline segment parameters,
  // Get y values for each spline knot,
  // Every knot must evaluate to 1.0 to create a flat spline
  for(int i = 0; i < Np; i++) {
    spl->GetCoeff(i, x, y, b, c, d);
    if (y != 1) {
      return false;
    }
  }
  return true;
}

// *********************************
/// @brief CW: Reduced the TSpline3 to TSpline3_red
/// @param MasterSpline Vector of TSpline3_red pointers which we strip back
inline std::vector<std::vector<TSpline3_red*> > ReduceTSpline3(std::vector<std::vector<TSpline3*> > &MasterSpline) {
// *********************************
  std::vector<std::vector<TSpline3*> >::iterator OuterIt;
  std::vector<TSpline3*>::iterator InnerIt;

  // The return vector
  std::vector<std::vector<TSpline3_red*> > ReducedVector;
  ReducedVector.reserve(MasterSpline.size());

  // Loop over each parameter
  int OuterCounter = 0;
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt, ++OuterCounter) {
    // Make the temp vector
    std::vector<TSpline3_red*> TempVector;
    TempVector.reserve(OuterIt->size());
    int InnerCounter = 0;
    // Loop over each TSpline3 pointer
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt, ++InnerCounter) {
      // Here's our delicious TSpline3 object
      TSpline3 *spline = (*InnerIt);
      // Now make the reduced TSpline3 pointer
      TSpline3_red *red = NULL;
      if (spline != NULL) {
        red = new TSpline3_red(spline);
        (*InnerIt) = spline;
      }
      // Push back onto new vector
      TempVector.push_back(red);
    } // End inner for loop
    ReducedVector.push_back(TempVector);
  } // End outer for loop
  // Now have the reduced vector
  return ReducedVector;
}

// *********************************
/// @brief CW: Reduced the TF1 to TF1_red
/// @param MasterSpline Vector of TF1_red pointers which we strip back
inline std::vector<std::vector<TF1_red*> > ReduceTF1(std::vector<std::vector<TF1*> > &MasterSpline) {
// *********************************
  std::vector<std::vector<TF1*> >::iterator OuterIt;
  std::vector<TF1*>::iterator InnerIt;

  // The return vector
  std::vector<std::vector<TF1_red*> > ReducedVector;
  ReducedVector.reserve(MasterSpline.size());

  // Loop over each parameter
  int OuterCounter = 0;
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt, ++OuterCounter) {
    // Make the temp vector
    std::vector<TF1_red*> TempVector;
    TempVector.reserve(OuterIt->size());
    int InnerCounter = 0;
    // Loop over each TSpline3 pointer
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt, ++InnerCounter) {
      // Here's our delicious TSpline3 object
      TF1* spline = (*InnerIt);
      // Now make the reduced TSpline3 pointer (which deleted TSpline3)
      TF1_red* red = NULL;
      if (spline != NULL) {
        red = new TF1_red(spline);
        (*InnerIt) = spline;
      }
      // Push back onto new vector
      TempVector.push_back(red);
    } // End inner for loop
    ReducedVector.push_back(TempVector);
  } // End outer for loop
  // Now have the reduced vector
  return ReducedVector;
}

#pragma GCC diagnostic pop

