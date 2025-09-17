#pragma once

/// @file Core.h
/// @author Clarence Wret
/// @author Kamil Skwarczynski
/// @author Luke Pickering

// C++ includes
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>

/// Run low or high memory versions of structs
/// N.B. for 64 bit systems sizeof(float) == sizeof(double) so not a huge effect
namespace M3 {
  #ifdef _LOW_MEMORY_STRUCTS_
  /// Custom floating point (float or double)
  using float_t = float;
  /// Custom integer (int or short int)
  using int_t = short;
  /// Custom unsigned integer (unsigned short int or unsigned int)
  using uint_t = unsigned short;
  /// Custom variable allowing us to save in either double or float
  constexpr static const char* float_t_str_repr = "F";
  #else
  using float_t = double;
  using int_t = int;
  using uint_t = unsigned;
  constexpr static const char* float_t_str_repr = "D";
  #endif

  /// Function template for fused multiply-add
  template <typename T>
  constexpr T fmaf_t(T x, T y, T z) {
    #ifdef _LOW_MEMORY_STRUCTS_
    return std::fmaf(x, y, z);
    #else
    return std::fma(x, y, z);
    #endif
  }
  /// Default value used for double initialisation
  constexpr static const double _BAD_DOUBLE_ = -999.99;
  /// Default value used for int initialisation
  constexpr static const int _BAD_INT_ = -999;
  constexpr static const double _DEFAULT_RETURN_VAL_ = -999999.123456;

  /// Some commonly used variables to which we set pointers to
  constexpr static const double Unity_D = 1.;
  constexpr static const float Unity_F = 1.;
  #ifdef _LOW_MEMORY_STRUCTS_
  constexpr static const float_t Unity = Unity_F;
  #else
  constexpr static const float_t Unity = Unity_D;
  #endif
  constexpr static const double Zero_D = 0.;
  constexpr static const float Zero_F = 0.;
  #ifdef _LOW_MEMORY_STRUCTS_
  constexpr static const float_t Zero = Zero_F;
  #else
  constexpr static const float_t Zero = Zero_D;
  #endif

  /// When parameter has no bound this serves as it. Lowest possible value the system
  constexpr static const double KinematicLowBound = std::numeric_limits<double>::lowest();
  /// When parameter has no bound this serves as it. Highest possible value the system
  constexpr static const double KinematicUpBound = std::numeric_limits<double>::max();

  /// Large Likelihood is used it parameter go out of physical boundary, this indicates in MCMC that such step should be removed
  constexpr static const double _LARGE_LOGL_ = 1234567890.0;

  /// MC prediction lower bound in bin to identify problematic binning definitions and handle LogL calculation
  constexpr static const double _LOW_MC_BOUND_ = .00001;

  /// Default value for spline knot capping, default mean not capping is being applied
  constexpr static const double DefSplineKnotUpBound = 9999;
  /// Default value for spline knot capping, default mean not capping is being applied
  constexpr static const double DefSplineKnotLowBound = -9999;
}

/// KS: noexcept can help with performance but is terrible for debugging, this is meant to help easy way of of turning it on or off. In near future move this to struct or other central class.
#ifndef DEBUG
#define _noexcept_ noexcept
#else
#define _noexcept_
#endif

/// KS: Using restrict limits the effects of pointer aliasing, aiding optimizations. While reading I found that there might be some compilers which don't have __restrict__. As always we use _restrict_ to more easily turn off restrict in the code
#ifndef DEBUG
#define _restrict_ __restrict__
#else
#define _restrict_
#endif

/// Number of overflow bins in TH2Poly,
#define _TH2PolyOverflowBins_ 9

#ifdef MULTITHREAD
#include "omp.h"
#endif

/// @brief KS: Avoiding warning checking for headers
/// @details Many external files don't strictly adhere to rigorous C++ standards.
/// Inline functions in such headers may cause errors when compiled with strict MaCh3 compiler flags.
/// @warning Use this for any external header file to avoid warnings.
#define _MaCh3_Safe_Include_Start_ \
_Pragma("GCC diagnostic push") \
_Pragma("GCC diagnostic ignored \"-Wuseless-cast\"") \
_Pragma("GCC diagnostic ignored \"-Wfloat-conversion\"") \
_Pragma("GCC diagnostic ignored \"-Wold-style-cast\"") \
_Pragma("GCC diagnostic ignored \"-Wformat-nonliteral\"") \
_Pragma("GCC diagnostic ignored \"-Wswitch-enum\"") \
_Pragma("GCC diagnostic ignored \"-Wconversion\"") \
_Pragma("GCC diagnostic ignored \"-Wshadow\"") \
_Pragma("GCC diagnostic ignored \"-Wswitch-enum\"")
/// @brief KS: Restore warning checking after including external headers
#define _MaCh3_Safe_Include_End_ \
_Pragma("GCC diagnostic pop")

// clang and IntelLLVM need slightly different diagnostics
#if defined(__clang__) || defined(__INTEL_LLVM_COMPILER)
  #undef _MaCh3_Safe_Include_Start_
  #define _MaCh3_Safe_Include_Start_ \
  _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wfloat-conversion\"") \
  _Pragma("clang diagnostic ignored \"-Wold-style-cast\"") \
  _Pragma("clang diagnostic ignored \"-Wformat-nonliteral\"") \
  _Pragma("clang diagnostic ignored \"-Wswitch-enum\"") \
  _Pragma("clang diagnostic ignored \"-Wconversion\"") \
  _Pragma("clang diagnostic ignored \"-Wshadow\"") \
  _Pragma("clang diagnostic ignored \"-Wdeprecated-literal-operator\"") \
  _Pragma("clang diagnostic ignored \"-Wswitch-enum\"")
  #undef _MaCh3_Safe_Include_End_
  #define _MaCh3_Safe_Include_End_ \
  _Pragma("clang diagnostic pop")
#endif
