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
// NuOscillator includes
#include "Constants/OscillatorConstants.h"
_MaCh3_Safe_Include_End_ //}

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
/// @brief KS: This is mad way of converting string to int. Why? To be able to use string with switch
constexpr unsigned int str2int(const char* str, int h = 0) {
// *******************
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

// *******************
/// @brief ETA - Normalisations for cross-section parameters
/// Carrier for whether you want to apply a systematic to an event or not
struct NormParameter {
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
  /// within a SampleHandler daughter class
  /// The bounds for each kinematic variable are given in Selection
  std::vector< std::string > KinematicVarStr;

  /// Parameter number of this normalisation in current systematic model
  int index;
};

// HH - a shorthand type for funcpar functions
using FuncParFuncType = std::function<void (const double*, std::size_t, std::size_t)>;
// *******************
/// @brief HH - Functional parameters
/// Carrier for whether you want to apply a systematic to an event or not
struct FunctionalParameter {
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
  /// within a SampleHandler daughter class
  /// The bounds for each kinematic variable are given in Selection
  std::vector< std::string > KinematicVarStr;

  /// Parameter number of this functional in current systematic model
  int index;

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
  kSystTypes  //!< This only enumerates
};

// *******************
/// @brief KS: Struct holding info about Spline Systematics
struct SplineParameter {
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

// *****************
/// Enum to track the target material
enum TargetMat {
// *****************
  kTarget_H  = 1,    //!< Hydrogen (Atomic number 1)
  kTarget_C  = 12,   //!< Carbon 12 (Atomic number 6)
  kTarget_N  = 14,   //!< Nitrogen (Atomic number 7)
  kTarget_O  = 16,   //!< Oxygen 16 (Atomic number 8)
  kTarget_Al = 27,   //!< Aluminum (Atomic number 13)
  kTarget_Ar = 40,   //!< Argon (Atomic number 18)
  kTarget_Ti = 48,   //!< Titanium (Atomic number 22)
  kTarget_Fe = 56,   //!< Iron (Atomic number 26)
  kTarget_Pb = 207   //!< Lead (Atomic number 82)
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
    case kTarget_Ar:
      name = "Argon";
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

/// Make an enum of the test statistic that we're using
/// @todo KS: Consider adding BakerCousins based on Baker & Cousins, Nucl.Instrum.Meth.A 221 (1984) 437-442
enum TestStatistic {
  kPoisson,                 //!< Standard Poisson likelihood
  kBarlowBeeston,           //!< Barlow-Beeston following Conway \cite Conway:2011in
  kIceCube,                 //!< Based on \cite Arguelles:2019izp
  kPearson,                 //!< Standard Pearson likelihood
  kDembinskiAbdelmotteleb,  //!< Based on \cite Dembinski:2022ios
  kNTestStatistics          //!< Number of test statistics
};

// **************************************************
/// @brief Convert a LLH type to a string
inline std::string TestStatistic_ToString(TestStatistic i) {
// **************************************************
  std::string name = "";

  switch(i) {
    case TestStatistic::kPoisson:
    name = "Poisson";
    break;
    case TestStatistic::kBarlowBeeston:
    name = "Barlow-Beeston";
    break;
    case TestStatistic::kIceCube:
    name = "IceCube";
    break;
    case TestStatistic::kPearson:
    name = "Pearson";
    break;
    case TestStatistic::kDembinskiAbdelmotteleb:
    name = "Dembinski-Abdelmotteleb";
    break;
    case TestStatistic::kNTestStatistics:
      MACH3LOG_ERROR("kNTestStatistics is not a valid TestStatistic!");
      throw MaCh3Exception(__FILE__, __LINE__);
    default:
      MACH3LOG_ERROR("UNKNOWN LIKELIHOOD SPECIFIED!");
      MACH3LOG_ERROR("You gave test-statistic {}", static_cast<int>(i));
      throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return name;
}

// ***************************
/// @brief KS: Small struct used for applying kinematic cuts
struct KinematicCut {
// ***************************
  /// Index or enum value identifying the kinematic variable to cut on.
  int ParamToCutOnIt = M3::_BAD_INT_;
  /// Lower bound on which we apply cut
  double LowerBound = M3::_BAD_DOUBLE_;
  /// Upper bound on which we apply cut
  double UpperBound = M3::_BAD_DOUBLE_;
};

// ***************************
// A handy namespace for variables extraction
namespace MaCh3Utils {
// ***************************

  /// @brief Return mass for given PDG
  /// @note Get the mass of a particle from the PDG In GeV, not MeV!
  /// @cite pdg2024 (particle masses)
  /// @cite ame2020 (nuclear masses)
  inline double GetMassFromPDG(const int PDG) {
    // *****************************
    switch (abs(PDG)) {
      // Leptons
      case 11: return 0.00051099895; // e
      case 13: return 0.1056583755;  // mu
      case 15: return 1.77693;       // tau
      // Neutrinos
      case 12:
      case 14:
      case 16:
        return 0.; 
      // Photon
      case 22: return 0.; 
      // Mesons
      case 211: return 0.13957039; // pi_+/-
      case 111: return 0.1349768;  // pi_0
      case 221: return 0.547862;   // eta
      case 311:                    // K_0
      case 130:                    // K_0_L
      case 310:                    // K_0_S
        return 0.497611; 
      case 321: return 0.493677;   // K_+/-
      // Baryons
      case 2112: return 0.939565; // n
      case 2212: return 0.938272; // p
      case 3122: return 1.115683; // lambda
      case 3222: return 1.118937; // sig_+
      case 3112: return 1.197449; // sig_-
      case 3212: return 1.192642; // sig_0
      // Nuclei
      case 1000050110: return 10.255103; // Boron-11
      case 1000060120: return 11.177929; // Carbon-12
      case 1000070140: return 13.043781; // Nitrogen-14
      case 1000080160: return 14.899169; // Oxygen-16
      case 1000090190: return 17.696901; // Fluorine-19
      case 1000110230: return 21.414835; // Sodium-23
      case 1000130270: return 25.133144; // Aluminum-27
      case 1000140280: return 26.060342; // Silicon-28
      case 1000190390: return 36.294463; // Potassium-39
      case 1000180400: return 37.224724; // Argon-40
      case 1000220480: return 44.663224; // Titanium-48
      case 1000300640: return 59.549619; // Zinc-64
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
    int NuOscillatorFlavour = M3::_BAD_INT_;
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
        MACH3LOG_ERROR("Unknown Neutrino PDG {}, cannot convert to NuOscillator type", NuPdg);
        break;
    }

    //This is very cheeky but if the PDG is negative then multiply the PDG by -1
    // This is consistent with the treatment that NuOscillator expects as enums only
    // exist for the generic matter flavour and not the anti-matter version
    if(NuPdg < 0){NuOscillatorFlavour *= -1;}

    return NuOscillatorFlavour;
  }
  // ***************************

  // ***************************
  /// @brief Convert double into string for precision, useful for playing with yaml if you don't want to have in config floating point precision...
  inline std::string FormatDouble(const double value, const int precision) {
  // ***************************
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
  }

} // end MaCh3Utils namespace
