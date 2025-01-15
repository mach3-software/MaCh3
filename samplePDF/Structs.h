#pragma once

/// Run low or high memory versions of structs
/// N.B. for 64 bit systems sizeof(float) == sizeof(double) so not a huge effect
/// KS: Need more testing on FD
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
}

/// KS: noexcept can help with performance but is terrible for debugging, this is meant to help easy way of of turning it on or off. In near future move this to struct or other central class.
//#define SafeException
#ifndef SafeException
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

constexpr static const double _BAD_DOUBLE_ = -999.99;
constexpr static const int _BAD_INT_ = -999;

constexpr static const double _DEFAULT_RETURN_VAL_ = -999999.123456;

// Some commonly used variables to which we set pointers to
constexpr static const double Unity = 1.;
constexpr static const double Zero = 0.;
constexpr static const float Unity_F = 1.;
constexpr static const float Zero_F = 0.;
constexpr static const int Unity_Int = 1;

// C++ includes
#include <sstream>
#include <fstream> 
#include <iostream>
#include <vector>
#include <iomanip>
#include <set>
#include <list>
#include <unordered_map>

#ifdef MULTITHREAD
#include "omp.h"
#endif

#include "manager/MaCh3Exception.h"
#include "manager/MaCh3Logger.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
// ROOT include
#include "TSpline.h"
#include "TObjString.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TH1.h"
// NuOscillator includes
#include "Constants/OscillatorConstants.h"
#pragma GCC diagnostic pop


/// @file Structs.h
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
/// @brief KS: This is mad way of converting string to int. Why? To be able to use string with switch
constexpr unsigned int str2int(const char* str, int h = 0) {
// *******************
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

// *******************
/// @brief ETA - Normalisations for cross-section parameters
/// Carrier for whether you want to apply a systematic to an event or not
struct XsecNorms4 {
// *******************
  /// Name of parameters
  std::string name;
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
  std::vector< std::vector<double> > Selection;

  /// Generic vector containing the string of kinematic type
  /// This then needs to be converted to a kinematic type enum
  /// within a samplePDF daughter class
  /// The bounds for each kinematic variable are given in Selection
  std::vector< std::string > KinematicVarStr;

  /// Parameter number of this normalisation in current systematic model
  int index;
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
    //  TSpline3 (third order spline in ROOT)
    case SplineInterpolation::kTSpline3:
      Type = RespFuncType::kTSpline3_red;
      break;
    case SplineInterpolation::kLinear:
      Type = RespFuncType::kTSpline3_red;
      break;
    case SplineInterpolation::kMonotonic:
      Type = RespFuncType::kTSpline3_red;
      break;
    //  (Experimental) Akima_Spline (crd order spline which is allowed to be discontinuous in 2nd deriv)
    case SplineInterpolation::kAkima:
      Type = RespFuncType::kTSpline3_red;
      break;
    case SplineInterpolation::kLinearFunc:
      Type = RespFuncType::kTF1_red;
      break;
    default:
      MACH3LOG_ERROR("UNKNOWN SPLINE INTERPOLATION SPECIFIED!");
      MACH3LOG_ERROR("You gave {}", static_cast<int>(i));
      throw MaCh3Exception(__FILE__ , __LINE__ );
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
  kSystTypes  //!< This only enumerates
};

// *******************
/// @brief KS: Struct holding info about Spline Systematics
struct XsecSplines1 {
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
    default:
      MACH3LOG_ERROR("UNKNOWN SYST TYPE SPECIFIED!");
      MACH3LOG_ERROR("You gave {}", static_cast<int>(i));
      throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return name;
}

// *****************
/// Enum to track the target material
enum TargetMat {
// *****************
  kTarget_H  = 1,   //!< Hydrogen
  kTarget_C  = 12,  //!< Carbon 12
  kTarget_N  = 14,
  kTarget_O  = 16,  //!< Oxygen 16
  kTarget_Al = 27,
  kTarget_Ar = 40,
  kTarget_Ti = 48,
  kTarget_Fe = 56,
  kTarget_Pb = 207
};

// *****************
/// @brief Converted the Target Mat to a string
inline std::string TargetMat_ToString(const TargetMat i) {
// *****************
  std::string name;

  switch(i) {
    case kTarget_H:
      name = "Hydrogen";
      break;
    case kTarget_C:
      name = "Carbon";
      break;
    case kTarget_N:
      name = "Nitrogen";
      break;
    case kTarget_O:
      name = "Oxygen";
      break;
    case kTarget_Al:
      name = "Aluminium";
      break;
    case kTarget_Ti:
      name = "Titanium";
      break;
    case kTarget_Fe:
      name = "Iron";
      break;
    case kTarget_Pb:
      name = "Lead";
      break;
    default:
      name = "TargetMat_Undefined";
      break;
  }

  return name;
}

// *****************
/// Enum to track the incoming neutrino species
enum NuPDG {
// *****************
  kNue = 12,
  kNumu = 14,
  kNutau = 16,
  kNueBar = -12,
  kNumuBar = -14,
  kNutauBar = -16
};

// *****************
/// Enum to track neutrino species for Prob3
enum ProbNu {
// *****************
  kProbNue = 1,
  kProbNumu = 2,
  kProbNutau = 3,
  kProbNueBar = -1,
  kProbNumuBar = -2,
  kProbNutauBar = -3
};

/// @brief ETA - Probs3++ doesn't use the PDG codes for the neutrino type so add in a small converter
inline int PDGToProbs(NuPDG pdg){

  int ReturnProbNu = -999;

  switch (pdg){
    case kNue:
      ReturnProbNu = kProbNue;
      break;
    case kNumu:
      ReturnProbNu = kProbNumu;
      break;
    case kNutau:
      ReturnProbNu = kProbNutau;
      break;
    case kNueBar:
      ReturnProbNu = kProbNueBar;
      break;
    case kNumuBar:
      ReturnProbNu = kProbNumuBar;
      break;
    case kNutauBar:
      ReturnProbNu = kProbNutauBar;
      break;
    default:
      MACH3LOG_WARN("Unrecognised pdg for the neutrino so can't map this to an int for Prob3++");
      break;
  }

  return ReturnProbNu;
}

inline int ProbsToPDG(ProbNu NuType){
  int ReturnNuPDG = -999;

  switch (NuType){
    case kProbNue:
      ReturnNuPDG = static_cast<int>(kNue);
      break;
    case kProbNumu:
      ReturnNuPDG = static_cast<int>(kNumu);
      break;
    case kProbNutau:
      ReturnNuPDG = static_cast<int>(kNutau);
      break;
    case kProbNueBar:
      ReturnNuPDG = static_cast<int>(kNueBar);
      break;
    case kProbNumuBar:
      ReturnNuPDG = static_cast<int>(kNumuBar);
      break;
    case kProbNutauBar:
      ReturnNuPDG = static_cast<int>(kNutauBar);
      break;
    default:
      MACH3LOG_WARN("Unrecognised NuType for the neutrino so can't map this to a PDG code");
      break;
  }
  return ReturnNuPDG;
}

/// Make an enum of the test statistic that we're using
enum TestStatistic {
  kPoisson,                //!< Standard Poisson likelihood
  kBarlowBeeston,          //!< Barlow-Beeston following Conway \cite Conway:2011in
  kIceCube,                //!< Based on \cite Arguelles:2019izp
  kPearson,                //!< Standard Pearson likelihood
  kDembinskiAbdelmottele,  //!< Based on \cite Dembinski:2022ios
  kNTestStatistics         //!< Number of test statistics
};

// **************************************************
/// @brief Convert a LLH type to a string
inline std::string TestStatistic_ToString(TestStatistic i) {
// **************************************************
  std::string name = "";

  switch(i) {
    case kPoisson:
    name = "Poisson";
    break;
    case kBarlowBeeston:
    name = "BarlowBeeston";
    break;
    case kIceCube:
    name = "IceCube";
    break;
    case kPearson:
    name = "Pearson";
    break;
    case kDembinskiAbdelmottele:
    name = "DembinskiAbdelmottele";
    break;
    default:
      MACH3LOG_ERROR("UNKNOWN LIKELIHOOD SPECIFIED!");
      MACH3LOG_ERROR("You gave test-statistic {}", static_cast<int>(i));
      throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return name;
}


// ***************************
// A handy namespace for variables extraction
namespace MaCh3Utils {
  // ***************************

  // ***************************
  /// @brief Return mass for given PDG
  // *****************************
  // Get the mass of a particle from the PDG
  // In GeV, not MeV!
  inline double GetMassFromPDG(const int PDG) {
    // *****************************
    switch (abs(PDG)) {
      case 11:
        return 0.511E-3;
        break;
      case 13:
        return 105.658E-3;
        break;
      case 15:
        return 1.77682;
        break;
      case 22:
        return 0.;
        break;
      case 211:
        return 139.57E-3;
        break;
      case 111:
        return 134.98E-3;
        break;
      case 2112:
        return 939.565E-3;
        break;
      case 2212:
        return 938.27E-3;
        break;
        //Oxygen nucleus
      case 1000080160:
        return 14.89926;
        break;
        //eta
      case 221:
        return 547.862E-3;
        break;
        //K^0 (s or l)
      case 311:
      case 130:
      case 310:
        return 497.611E-3;
        break;
      case 321:
        return 493.677E-3;
        break;
        // Lambda baryon
      case 3122:
        return 1115.683E-3;
        break;
      case 12:
      case 14:
      case 16:
        return 0.0;
        break;
      default:
        MACH3LOG_ERROR("Haven't got a saved mass for PDG: {}", PDG);
        MACH3LOG_ERROR("Please implement me!");
        throw MaCh3Exception(__FILE__, __LINE__);
    } // End switch
    MACH3LOG_ERROR("Warning, didn't catch a saved mass");
    return 0;
  }
  // ***************************

  // ***************************
  /// @brief Convert from PDG flavour to NuOscillator type
  /// beware that in the case of anti-neutrinos the NuOscillator
  /// type simply gets multiplied by -1
  inline int PDGToNuOscillatorFlavour(int NuPdg){
    int NuOscillatorFlavour = _BAD_INT_;
    switch(std::abs(NuPdg)){
      case NuPDG::kNue:
        NuOscillatorFlavour = NuOscillator::kElectron;
        break;
      case NuPDG::kNumu:
        NuOscillatorFlavour = NuOscillator::kMuon;
        break;
      case NuPDG::kNutau:
        NuOscillatorFlavour = NuOscillator::kTau;
        break;
      default:
        MACH3LOG_ERROR("Unknown Nuetrino PDG {}, cannot convert to NuOscillator type", NuPdg);
        break;
    }

    //This is very cheeky but if the PDG is negative then multiply the PDG by -1
    // This is consistent with the treatment that NuOscillator expects as enums only
    // exist for the generic matter flavour and not the anti-matter version
    if(NuPdg < 0){NuOscillatorFlavour *= -1;}

    return NuOscillatorFlavour;
  }
  // ***************************
} // end MaCh3Utils namespace
