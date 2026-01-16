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
constexpr unsigned int str2int(const char* str, const int h = 0) {
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
enum TestStatistic {
  kPoisson,                 //!< Standard Poisson likelihood @cite BakerCousins1984
  kBarlowBeeston,           //!< Barlow-Beeston (@cite Barlow:1993dm) following Conway approximation (@cite Conway:2011in)
  kIceCube,                 //!< Based on @cite Arguelles:2019izp
  kPearson,                 //!< Standard Pearson likelihood @cite Pearson1900
  kDembinskiAbdelmotteleb,  //!< Based on @cite Dembinski:2022ios
  kNTestStatistics          //!< Number of test statistics
};

// **************************************************
/// @brief Convert a LLH type to a string
inline std::string TestStatistic_ToString(const TestStatistic TestStat) {
// **************************************************
  std::string name = "";

  switch(TestStat) {
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
      MACH3LOG_ERROR("You gave test-statistic {}", static_cast<int>(TestStat));
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
/// @brief KS: Store bin lookups allowing to quickly find bin after migration
struct BinShiftLookup {
// ***************************
  /// lower to check if Eb has moved the erec bin
  double lower_binedge;
  /// lower to check if Eb has moved the erec bin
  double lower_lower_binedge;
  /// upper to check if Eb has moved the erec bin
  double upper_binedge;
  /// upper to check if Eb has moved the erec bin
  double upper_upper_binedge;
};

// ***************************
/// @brief KS: Small struct storying info about used binning
struct SampleBinningInfo {
// ***************************
  /// Vector to hold N-axis bin-edges
  std::vector<std::vector<double>> BinEdges;
  /// Number of N-axis bins in the histogram used for likelihood calculation
  std::vector<size_t> AxisNBins;

  /// Number of total bins
  size_t nBins = M3::_BAD_INT_;
  /// If you have binning for multiple samples and trying to define 1D vector let's
  size_t GlobalOffset = M3::_BAD_INT_;
  /// Bin lookups for all dimensions
  std::vector <std::vector<BinShiftLookup> > BinLookup;
  /// Stride factors for converting N-dimensional bin indices to a linear index.
  std::vector<int> Strides;

  /// @brief Get linear bin index from 2D bin indices
  /// @param xBin The bin index along the X axis (0-based)
  /// @param yBin The bin index along the Y axis (0-based)
  /// @return The linear bin index corresponding to (xBin, yBin)
  int GetBin(const int xBin, const int yBin) const {
    return static_cast<int>(yBin * AxisNBins[0] + xBin);
  }

  /// @brief Convert N-dimensional bin indices to a linear bin index.
  /// @param Bins Vector of bin indices along each dimension
  /// @details
  /// Mapping follows row-major order:
  ///  - 1D: xBin
  ///  - 2D: yBin * nXBins + xBin
  ///  - 3D: zBin * nXBins * nYBins + yBin * nXBins + xBin
  int GetBin(const std::vector<int>& Bins) const {
    int BinNumber = 0;
    for(size_t i = 0; i < Bins.size(); ++i) {
      BinNumber += Bins[i]*Strides[i];
    }
    return BinNumber;
  }

  /// @brief Get linear bin index from 2D bin indices with additional checks
  /// @param xBin The bin index along the X axis (0-based)
  /// @param yBin The bin index along the Y axis (0-based)
  /// @return The linear bin index corresponding to (xBin, yBin)
  /// @warning this performs additional checks so do not use in parts of code used during fit
  int GetBinSafe(const int xBin, const int yBin) const {
    if (xBin < 0 || yBin < 0 || static_cast<size_t>(xBin) >= AxisNBins[0] || static_cast<size_t>(yBin) >= AxisNBins[1]) {
      MACH3LOG_ERROR("GetBinSafe: Bin indices out of range: xBin={}, yBin={}, max xBin={}, max yBin={}",
                     xBin, yBin, AxisNBins[0] - 1, AxisNBins[1] - 1);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    return GetBin(xBin, yBin);
  }

  /// @brief Calculates the global bin number for a given 2D bin, accounting for multiple binning samples.
  /// @param xBin The bin index along the X axis (0-based)
  /// @param yBin The bin index along the Y axis (0-based)
  int GetBinGlobal(const int xBin, const int yBin) const {
    return static_cast<int>(GlobalOffset + GetBin(xBin, yBin));
  }

  /// @brief DB Find the relevant bin in the PDF for each event
  int FindBin(const int Dimension, const double Var, const int NomBin) const {
    return FindBin(Var, NomBin, AxisNBins[Dimension], BinEdges[Dimension], BinLookup[Dimension]);
  }

  /// @brief DB Find the relevant bin in the PDF for each event
  /// @param KinVar       The value of the kinematic variable for the event.
  /// @param NomBin       The nominal bin index where the event would fall without any shifts.
  /// @param N_Bins       The total number of bins in this dimension.
  /// @param Bin_Edges    Vector of bin edge values (size = N_Bins + 1).
  /// @param Bin_Lookup   Vector of `BinShiftLookup` structs providing precomputed lower and upper
  ///                     edges for the nominal bin and its neighbors to efficiently handle
  ///                     shifted events
  int FindBin(const double KinVar,
              const int NomBin,
              const size_t N_Bins,
              const std::vector<double>& Bin_Edges,
              const std::vector<BinShiftLookup>& Bin_Lookup) const {
    //DB Check to see if momentum shift has moved bins
    //DB - First , check to see if the event is outside of the binning range and skip event if it is
    if (KinVar < Bin_Edges[0] || KinVar >= Bin_Edges[N_Bins]) {
      return M3::UnderOverFlowBin;
    }

    // KS: Get reference to avoid repeated indexing and help with performance
    const BinShiftLookup& _restrict_ Bin = Bin_Lookup[NomBin];
    const double lower = Bin.lower_binedge;
    const double upper = Bin.upper_binedge;
    const double lower_lower = Bin.lower_lower_binedge;
    const double upper_upper = Bin.upper_upper_binedge;

    //DB - Second, check to see if the event is still in the nominal bin
    if (KinVar < upper && KinVar >= lower) {
      return NomBin;
    }
    //DB - Thirdly, check the adjacent bins first as Eb+CC+EScale shifts aren't likely to move an Erec more than 1bin width
    //Shifted down one bin from the event bin at nominal
    if (KinVar < lower && KinVar >= lower_lower) {
      return NomBin-1;
    }
    //Shifted up one bin from the event bin at nominal
    if (KinVar < upper_upper && KinVar >= upper) {
      return NomBin+1;
    }
    //DB - If we end up in this loop, the event has been shifted outside of its nominal bin, but is still within the allowed binning range
    // KS: Perform binary search to find correct bin. We already checked if isn't outside of bounds
    return static_cast<int>(std::distance(Bin_Edges.begin(), std::upper_bound(Bin_Edges.begin(), Bin_Edges.end(), KinVar)) - 1);
  }

  /// @brief Initializes lookup arrays for efficient bin migration in a single dimension.
  /// @param Bin_Lookup Reference to the BinShiftLookup struct to be initialized.
  /// @param Bin_Edges Vector of bin edges defining the bin boundaries.
  /// @param TotBins Number of bins in the dimension.
  void InitialiseLookUpSingleDimension(std::vector<BinShiftLookup>& Bin_Lookup, const std::vector<double>& Bin_Edges, const size_t TotBins) {
    Bin_Lookup.resize(TotBins);
    //Set rw_pdf_bin and upper_binedge and lower_binedge for each skmc_base
    for(size_t bin_i = 0; bin_i < TotBins; bin_i++){
      double low_lower_edge = M3::_DEFAULT_RETURN_VAL_;
      double low_edge = Bin_Edges[bin_i];
      double upper_edge = Bin_Edges[bin_i+1];
      double upper_upper_edge = M3::_DEFAULT_RETURN_VAL_;

      if (bin_i == 0) {
        low_lower_edge = Bin_Edges[0];
      } else {
        low_lower_edge = Bin_Edges[bin_i-1];
      }

      if (bin_i + 2 < TotBins) {
        upper_upper_edge = Bin_Edges[bin_i + 2];
      } else if (bin_i + 1 < TotBins) {
        upper_upper_edge = Bin_Edges[bin_i + 1];
      }

      Bin_Lookup[bin_i].lower_binedge = low_edge;
      Bin_Lookup[bin_i].upper_binedge = upper_edge;
      Bin_Lookup[bin_i].lower_lower_binedge = low_lower_edge;
      Bin_Lookup[bin_i].upper_upper_binedge = upper_upper_edge;
    }
  }

  /// @brief Initialise stride factors for linear bin index calculation.
  /// @details
  /// Strides define how N-dimensional bin indices are converted into a single
  /// linear bin index using row-major ordering.
  ///
  /// For example:
  ///  - 1D: stride[0] = 1
  ///  - 2D: stride[1] = nXBins
  ///  - 3D: stride[2] = nXBins * nYBins
  void InitialiseStrides(const int Dimension) {
    Strides.resize(Dimension);
    int stride = 1;
    for (int i = 0; i < Dimension; ++i) {
      Strides[i] = stride;
      // Multiply stride by the number of bins in this axis
      stride *= static_cast<int>(AxisNBins[i]);
    }
  }

  /// @brief Initialise special lookup arrays allowing to more efficiently perform bin-migration
  ///        These arrays store the lower and upper edges of each bin and their neighboring bins.
  void InitialiseBinMigrationLookUp(const int Dimension) {
    BinLookup.resize(Dimension);
    for (int i = 0; i < Dimension; ++i) {
      InitialiseLookUpSingleDimension(BinLookup[i], BinEdges[i], AxisNBins[i]);
    }

    InitialiseStrides(Dimension);
  }
};

/// @brief Get the sample index corresponding to a global bin number.
/// @param BinningInfo Vector of SampleBinningInfo structs.
/// @param GlobalBin The global bin number.
/// @return The index of the sample, or garbage if not found.
inline int GetSampleFromGlobalBin(const std::vector<SampleBinningInfo>& BinningInfo, const size_t GlobalBin) {
  for (size_t iSample = 0; iSample < BinningInfo.size(); ++iSample) {
    const SampleBinningInfo& info = BinningInfo[iSample];
    if (GlobalBin >= info.GlobalOffset && GlobalBin < info.GlobalOffset + info.nBins) {
      return static_cast<int>(iSample);
    }
  }

  MACH3LOG_ERROR("Couldn't find sample corresponding to bin {}", GlobalBin);
  throw MaCh3Exception(__FILE__, __LINE__);

  // GlobalBin is out of range for all samples
  return M3::_BAD_INT_;
}

/// @brief Get the local (sample) bin index from a global bin number.
/// @param BinningInfo Vector of SampleBinningInfo structs.
/// @param GlobalBin The global bin number.
/// @return The bin index within the sample.
inline int GetLocalBinFromGlobalBin(const std::vector<SampleBinningInfo>& BinningInfo,
                                    const size_t GlobalBin) {
  for (size_t iSample = 0; iSample < BinningInfo.size(); ++iSample) {
    const SampleBinningInfo& info = BinningInfo[iSample];

    if (GlobalBin >= info.GlobalOffset &&
      GlobalBin <  info.GlobalOffset + info.nBins)
    {
      return static_cast<int>(GlobalBin - info.GlobalOffset);
    }
  }

  MACH3LOG_ERROR("Couldn't find local bin corresponding to bin {}", GlobalBin);
  throw MaCh3Exception(__FILE__, __LINE__);

  return M3::_BAD_INT_;
}

// ***************************
// A handy namespace for variables extraction
namespace MaCh3Utils {
// ***************************
  // *****************************
  /// @brief Return mass for given PDG
  /// @note Get the mass of a particle from the PDG In GeV, not MeV!
  /// @todo this could be constexpr in c++17
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
  inline int PDGToNuOscillatorFlavour(const int NuPdg){
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
