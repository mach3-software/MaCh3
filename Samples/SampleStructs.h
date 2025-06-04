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

/// @file SampleStructs.h
/// @author Asher Kaboth
/// @author Clarence Wret
/// @author Patrick Dunne
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski


// *******************
/// @brief KS: This is mad way of converting string to int. Why? To be able to use string with switch
constexpr unsigned int str2int(const char* str, int h = 0) {
// *******************
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
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
  /// JM Tag of whether this is a event or sub-event -level variable (default to event-level)
  bool isSubEventVar = false;
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
