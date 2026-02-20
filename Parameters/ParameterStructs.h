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
/// @brief Definitions of generic parameter structs and utility templates for MaCh3.
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

/// @brief Base case: do nothing for non-vector types.
template <typename T>
void CleanVector(T&) {}

// *******************
/// @brief Recursively releases memory held by a std::vector (including nested).
///
/// Recursively clears a vector and requests release of its unused capacity.
/// `clear()` or `vec = {}` only destroy elements but keep the allocated buffer.
/// For deeply nested vectors we must clean children first, then call
/// `shrink_to_fit()` to ask the implementation to return the capacity.
template <typename T>
void CleanVector(std::vector<T>& v) {
// *******************
  for (auto& x : v)
    CleanVector(x);

  v.clear();
  v.shrink_to_fit();
}

// *******************
/// @brief Base case: do nothing for non-pointer types.
template <typename T>
void CleanContainer(T&) {}
// *******************

// *******************
/// @brief Deletes an owned raw pointer and nulls it.
///
/// Only use when ownership is explicit (allocated with `new`).
template <typename T>
void CleanContainer(T*& ptr) {
  // *******************
  delete ptr;
  ptr = nullptr;
}

// *******************
/// @brief Recursively deletes elements stored as pointers inside containers
/// and releases the container capacity.
///
/// This is required for structures like `std::vector<T*>` or deeply nested
/// pointer containers. Without deleting first, memory is leaked; without
/// `shrink_to_fit()`, capacity is retained.
template <typename T>
void CleanContainer(std::vector<T>& container) {
  // *******************
  for (auto& element : container)
    CleanContainer(element);

  container.clear();
  container.shrink_to_fit();
}

// *******************
/// @brief Base class storing info for parameters types, helping unify codebase
/// @ingroup SamplesAndParameters
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
/// @author Ed Atkin
struct NormParameter : public TypeParameterBase {
  /// Mode which parameter applies to
  std::vector<int> modes;
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

/// HH - a shorthand type for funcpar functions
using FuncParFuncType = std::function<void (const double*, std::size_t)>;
// *******************
/// @brief HH - Functional parameters
/// Carrier for whether you want to apply a systematic to an event or not
/// @author Hank Hua
struct FunctionalParameter : public TypeParameterBase {
// *******************
  /// Mode which parameter applies to
  std::vector<int> modes;
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
  kKochanekBartels,      //!< KS: Kochanek-Bartels spline: allows local control of tension, continuity, and bias
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
    case SplineInterpolation::kKochanekBartels:
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
    case SplineInterpolation::kKochanekBartels:
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
    case SplineInterpolation::kKochanekBartels:
      name = "KochanekBartels";
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
