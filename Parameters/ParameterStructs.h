#pragma once

// C++ includes
#include <set>
#include <list>
#include <unordered_map>

// MaCh3 includes
#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"
#include "Manager/Core.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT include
#include "TSpline.h"
#include "TObjString.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TH1.h"
_MaCh3_Safe_Include_End_ //}

/// @file ParameterStructs.h
/// @author Asher Kaboth
/// @author Clarence Wret
/// @author Patrick Dunne
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski

// *******************
/// @brief Template to make vector out of an array of any length
template< typename T, size_t N >
std::vector<T> MakeVector( const T (&data)[N] ) {
// *******************
  return std::vector<T>(data, data+N);
}

// *******************
/// @brief Generic cleanup function
template <typename T>
void CleanVector(std::vector<T>& vec) {
// *******************
  vec.clear();
  vec.shrink_to_fit();
}

// *******************
/// @brief Generic cleanup function
template <typename T>
void CleanVector(std::vector<std::vector<T>>& vec) {
// *******************
  for (auto& innerVec : vec) {
    innerVec.clear();
    innerVec.shrink_to_fit();
  }
  vec.clear();
  vec.shrink_to_fit();
}

// *******************
/// @brief Generic cleanup function
template <typename T>
void CleanVector(std::vector<std::vector<std::vector<T>>>& vec) {
// *******************
  for (auto& inner2DVec : vec) {
    for (auto& innerVec : inner2DVec) {
      innerVec.clear();
      innerVec.shrink_to_fit();
    }
    inner2DVec.clear();
    inner2DVec.shrink_to_fit();
  }
  vec.clear();
  vec.shrink_to_fit();
}

// *******************
/// @brief Generic cleanup function
/// @todo Use recursive to make it more scalable in future
template <typename T>
void CleanVector(std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>>> &vec) {
// *******************
  for (auto& v6 : vec) {
    for (auto& v5 : v6) {
      for (auto& v4 : v5) {
        for (auto& v3 : v4) {
          for (auto& v2 : v3) {
            for (auto& v1 : v2) {
              v1.clear();
              v1.shrink_to_fit();
            }
            v2.clear();
            v2.shrink_to_fit();
          }
          v3.clear();
          v3.shrink_to_fit();
        }
        v4.clear();
        v4.shrink_to_fit();
      }
      v5.clear();
      v5.shrink_to_fit();
    }
    v6.clear();
    v6.shrink_to_fit();
  }
  vec.clear();
  vec.shrink_to_fit();
}

// *******************
/// @brief Generic cleanup function
template <typename T>
void CleanContainer(std::vector<T*>& container) {
// *******************
  for (T* ptr : container) {
    if(ptr) delete ptr;
    ptr = nullptr;
  }
  container.clear();
  container.shrink_to_fit();
}

// *******************
/// @brief Generic cleanup function
template <typename T>
void CleanContainer(std::vector< std::vector< std::vector<T*> > >& container) {
// *******************
  for (auto& container2D : container) {
    for (auto& container1D : container2D) {
      for (T* ptr : container1D) {
        if(ptr) delete ptr;
      }
      container1D.clear();
      container1D.shrink_to_fit();
    }
    container2D.clear();
    container2D.shrink_to_fit();
  }
  container.clear();
  container.shrink_to_fit();
}

// *******************
/// @brief Base class storing info for parameters types, helping unify codebase
struct TypeParameterBase {
// *******************
  /// Name of parameters
  std::string name;

  /// Parameter number of this normalisation in current systematic model
  int index = M3::_BAD_INT_;
};

// *******************
/// @brief ETA - Normalisations for cross-section parameters
/// Carrier for whether you want to apply a systematic to an event or not
struct NormParameter : public TypeParameterBase {
  /// Mode which parameter applies to
  std::vector<int> modes;
  /// Horn currents which parameter applies to
  std::vector<int> horncurrents;
  /// PDG which parameter applies to
  std::vector<int> pdgs;
  /// Preosc PDG which parameter applies to
  std::vector<int> preoscpdgs;
  /// Targets which parameter applies to
  std::vector<int> targets;
  /// Does this parameter have kinematic bounds
  bool hasKinBounds;
  /// Generic vector contain enum relating to a kinematic variable
  /// and lower and upper bounds. This can then be passed to IsEventSelected
  std::vector< std::vector< std::vector<double> > > Selection;

  /// Generic vector containing the string of kinematic type
  /// This then needs to be converted to a kinematic type enum
  /// within a SampleHandler daughter class
  /// The bounds for each kinematic variable are given in Selection
  std::vector< std::string > KinematicVarStr;
};

// HH - a shorthand type for funcpar functions
using FuncParFuncType = std::function<void (const double*, std::size_t)>;
// *******************
/// @brief HH - Functional parameters
/// Carrier for whether you want to apply a systematic to an event or not
struct FunctionalParameter : public TypeParameterBase {
// *******************
  /// Mode which parameter applies to
  std::vector<int> modes;
  /// Horn currents which parameter applies to
  std::vector<int> horncurrents;
  /// PDG which parameter applies to
  std::vector<int> pdgs;
  /// Preosc PDG which parameter applies to
  std::vector<int> preoscpdgs;
  /// Targets which parameter applies to
  std::vector<int> targets;
  /// Does this parameter have kinematic bounds
  bool hasKinBounds;
  /// Generic vector contain enum relating to a kinematic variable
  /// and lower and upper bounds. This can then be passed to IsEventSelected
  std::vector< std::vector< std::vector<double> > > Selection;

  /// Generic vector containing the string of kinematic type
  /// This then needs to be converted to a kinematic type enum
  /// within a SampleHandler daughter class
  /// The bounds for each kinematic variable are given in Selection
  std::vector< std::string > KinematicVarStr;

  /// Parameter value pointer
  const double* valuePtr;

  /// Function pointer
  FuncParFuncType* funcPtr;
};

/// Make an enum of the spline interpolation type
enum RespFuncType {
  kTSpline3_red,  //!< Uses TSpline3_red for interpolation
  kTF1_red,       //!< Uses TF1_red for interpolation
  kRespFuncTypes  //!< This only enumerates
};

/// Make an enum of the spline interpolation type
enum SplineInterpolation {
  kTSpline3,             //!< Default TSpline3 interpolation
  kLinear,               //!< Linear interpolation between knots
  kMonotonic,            //!< EM: DOES NOT make the entire spline monotonic, only the segments
  kAkima,                //!< EM: Akima spline iis allowed to be discontinuous in 2nd derivative and coefficients in any segment
  kLinearFunc,           //!< Liner interpolation using TF1 not spline
  kSplineInterpolations  //!< This only enumerates
};

// **************************************************
/// @brief Get function for TF1_red
/// @param i Interpolation type
inline std::string GetTF1(const SplineInterpolation i) {
  // **************************************************
  std::string Func = "";
  switch(i) {
    case SplineInterpolation::kLinearFunc:
      Func = "([1]+[0]*x)";
      break;
    case SplineInterpolation::kTSpline3:
    case SplineInterpolation::kLinear:
    case SplineInterpolation::kMonotonic:
    case SplineInterpolation::kAkima:
    case SplineInterpolation::kSplineInterpolations:
      MACH3LOG_ERROR("Interpolation type {} not supported for TF1!", static_cast<int>(i));
      throw MaCh3Exception(__FILE__, __LINE__);
    default:
      MACH3LOG_ERROR("UNKNOWN SPLINE INTERPOLATION SPECIFIED!");
      MACH3LOG_ERROR("You gave {}", static_cast<int>(i));
      throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return Func;
}

// **************************************************
/// @brief Convert a RespFuncType type to a SplineInterpolation
/// @param i Interpolation type
inline RespFuncType SplineInterpolation_ToRespFuncType(const SplineInterpolation i) {
// **************************************************
  RespFuncType Type = kRespFuncTypes;
  switch(i) {
    case SplineInterpolation::kTSpline3:
    case SplineInterpolation::kLinear:
    case SplineInterpolation::kMonotonic:
    case SplineInterpolation::kAkima:
      Type = RespFuncType::kTSpline3_red;
      break;
    case SplineInterpolation::kLinearFunc:
      Type = RespFuncType::kTF1_red;
      break;
    case SplineInterpolation::kSplineInterpolations:
      MACH3LOG_ERROR("kSplineInterpolations is not a valid interpolation type!");
      throw MaCh3Exception(__FILE__, __LINE__);
    default:
      MACH3LOG_ERROR("UNKNOWN SPLINE INTERPOLATION SPECIFIED!");
      MACH3LOG_ERROR("You gave {}", static_cast<int>(i));
      throw MaCh3Exception(__FILE__, __LINE__);
  }
  return Type;
}

// **************************************************
/// @brief Convert a LLH type to a string
inline std::string SplineInterpolation_ToString(const SplineInterpolation i) {
// **************************************************
  std::string name = "";
  switch(i) {
    //  TSpline3 (third order spline in ROOT)
    case SplineInterpolation::kTSpline3:
      name = "TSpline3";
      break;
    case SplineInterpolation::kLinear:
      name = "Linear";
      break;
    case SplineInterpolation::kMonotonic:
      name = "Monotonic";
      break;
    //  (Experimental) Akima_Spline (crd order spline which is allowed to be discontinuous in 2nd deriv)
    case SplineInterpolation::kAkima:
      name = "Akima";
      break;
    case SplineInterpolation::kLinearFunc:
      name = "LinearFunc";
      break;
    case SplineInterpolation::kSplineInterpolations:
      MACH3LOG_ERROR("kSplineInterpolations is not a valid interpolation type!");
      throw MaCh3Exception(__FILE__, __LINE__);
    default:
      MACH3LOG_ERROR("UNKNOWN SPLINE INTERPOLATION SPECIFIED!");
      MACH3LOG_ERROR("You gave {}", static_cast<int>(i));
      throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return name;
}

/// Make an enum of systematic type recognised by covariance class
/// @todo KS: Consider using enum class, it is generally recommended as safer. It will require many static_cast
enum SystType {
  kNorm,      //!< For normalisation parameters
  kSpline,    //!< For splined parameters (1D)
  kFunc,      //!< For functional parameters
  kOsc,       //!< For oscillation parameters
  kSystTypes  //!< This only enumerates
};

// *******************
/// @brief KS: Struct holding info about Spline Systematics
struct SplineParameter : public TypeParameterBase {
// *******************
  /// Spline interpolation vector
  SplineInterpolation _SplineInterpolationType;

  /// Modes to which spline applies (valid only for binned splines)
  std::vector<int> _fSplineModes;

  /// EM: Cap spline knot lower value
  double _SplineKnotLowBound;
  /// EM: Cap spline knot higher value
  double _SplineKnotUpBound;
};

// **************************************************
/// @brief Convert a Syst type type to a string
inline std::string SystType_ToString(const SystType i) {
// **************************************************
  std::string name = "";
  switch(i) {
    case SystType::kNorm:
      name = "Norm";
      break;
    case SystType::kSpline:
      name = "Spline";
      break;
    case SystType::kFunc:
      name = "Functional";
      break;
    case SystType::kOsc:
      name = "Oscillation";
      break;
    case SystType::kSystTypes:
      MACH3LOG_ERROR("kSystTypes is not a valid SystType!");
      throw MaCh3Exception(__FILE__, __LINE__);
    default:
      MACH3LOG_ERROR("UNKNOWN SYST TYPE SPECIFIED!");
      MACH3LOG_ERROR("You gave {}", static_cast<int>(i));
      throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return name;
}

// *******************
/// @brief KS: Struct holding info about oscillation Systematics
/// @note right now it is empty
struct OscillationParameter : public TypeParameterBase {
// *******************

};
