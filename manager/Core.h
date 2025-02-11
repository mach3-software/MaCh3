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
  #else
  using float_t = double;
  using int_t = int;
  using uint_t = unsigned;
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
  constexpr static const double _BAD_DOUBLE_ = -999.99;
  constexpr static const int _BAD_INT_ = -999;
  constexpr static const double _DEFAULT_RETURN_VAL_ = -999999.123456;

  /// Some commonly used variables to which we set pointers to
  constexpr static const double Unity = 1.;
  constexpr static const double Zero = 0.;
  constexpr static const float Unity_F = 1.;
  constexpr static const float Zero_F = 0.;
  constexpr static const int Unity_Int = 1;

  /// Large Likelihood is used it parameter go out of physical boundary, this indicates in MCMC that such step should be removed
  constexpr static const double _LARGE_LOGL_ = 1234567890.0;
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
_Pragma("GCC diagnostic ignored \"-Wshadow\"")

/// @brief KS: Restore warning checking after including external headers
#define _MaCh3_Safe_Include_End_ \
_Pragma("GCC diagnostic pop")

// clang need slightly different diagnostics
#if defined(__clang__)
  #undef _MaCh3_Safe_Include_Start_
  #define _MaCh3_Safe_Include_Start_ \
  _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wfloat-conversion\"") \
  _Pragma("clang diagnostic ignored \"-Wold-style-cast\"") \
  _Pragma("clang diagnostic ignored \"-Wformat-nonliteral\"") \
  _Pragma("clang diagnostic ignored \"-Wswitch-enum\"") \
  _Pragma("clang diagnostic ignored \"-Wconversion\"") \
  _Pragma("clang diagnostic ignored \"-Wshadow\"")

  #undef _MaCh3_Safe_Include_End_
  #define _MaCh3_Safe_Include_End_ \
  _Pragma("clang diagnostic pop")
#endif
