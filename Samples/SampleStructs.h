#pragma once

// C++ includes
#include <set>
#include <list>
#include <unordered_map>

// MaCh3 includes
#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"
#include "Manager/Core.h"
#include "Parameters/ParameterStructs.h"

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
/// @brief Small struct used for applying shifts due to functional params
/// @author Hank Hua
struct FunctionalShifter {
// ***************************
  /// Pointer to parameter value
  const double* valuePtr = nullptr;
  /// Pointer to shifting function
  FuncParFuncType* funcPtr = nullptr;
};

// ***************************
/// @brief KS: Store bin lookups allowing to quickly find bin after migration
struct BinShiftLookup {
// ***************************
  /// lower to check if shift has moved the event to different bin
  double lower_binedge = M3::_BAD_DOUBLE_;
  /// lower to check if shift has moved the event to different bin
  double lower_lower_binedge = M3::_BAD_DOUBLE_;
  /// upper to check if shift has moved the event to different bin
  double upper_binedge = M3::_BAD_DOUBLE_;
  /// upper to check if shift has moved the event to different bin
  double upper_upper_binedge = M3::_BAD_DOUBLE_;
};

// ***************************
/// @brief KS: This hold bin extents in N-Dimensions allowing to check if Bin falls into
struct BinInfo {
// ***************************
  /// The extent of the bin, stored as {lower, upper} bounds for each dimension.
  /// Extent[d][0] = lower edge
  /// Extent[d][1] = upper edge
  std::vector<std::array<double, 2>> Extent;
  /// @brief Checks if a given event (point) falls inside the bin.
  bool IsEventInside(const std::vector<double>& KinVars) const _noexcept_ {
    bool inside = true;
    const size_t N = KinVars.size();
    #ifdef MULTITHREAD
    #pragma omp simd reduction(&:inside)
    #endif
    for(size_t i = 0; i < N; ++i) {
      const double Var = KinVars[i];
      const bool in_bin = (Var > Extent[i][0]) & (Var <= Extent[i][1]);
      inside &= in_bin;
    }
    return inside;
  }
  /// @brief Checks if a given event (point) falls inside the bin using pointer array
  bool IsEventInside(const std::vector<const double*>& KinVars) const _noexcept_ {
    bool inside = true;
    const size_t N = KinVars.size();
    #ifdef MULTITHREAD
    #pragma omp simd reduction(&:inside)
    #endif
    for (size_t i = 0; i < N; ++i) {
      const double Var = *KinVars[i];
      const bool in_bin = (Var > Extent[i][0]) & (Var <= Extent[i][1]);
      inside &= in_bin;
    }
    return inside;
  }
};

// ***************************
/// @brief KS: Struct storing all information required for sample binning
///
/// @details
/// This struct encapsulates the full binning definition for a single analysis
/// sample. It stores the bin edges, number of bins per dimension, stride factors,
/// and lookup tables required to efficiently map multi-dimensional kinematic
/// variables to linear bin indices.
///
/// @author Kamil Skwarczynski
struct SampleBinningInfo {
// ***************************
  /// Vector to hold N-axis bin-edges
  std::vector<std::vector<double>> BinEdges;
  /// Number of N-axis bins in the histogram used for likelihood calculation
  std::vector<int> AxisNBins;

  /// Number of total bins
  int nBins = M3::_BAD_INT_;
  /// If you have binning for multiple samples and trying to define 1D vector let's
  int GlobalOffset = M3::_BAD_INT_;
  /// Bin lookups for all dimensions
  std::vector <std::vector<BinShiftLookup> > BinLookup;
  /// Stride factors for converting N-dimensional bin indices to a linear index.
  std::vector<int> Strides;
  /// Tells whether to use inform binning grid or non-uniform
  bool Uniform = true;
  /// Bins used only for non-uniform
  std::vector<BinInfo> Bins;
  /// This grid tells what bins are associated with with what BinEdges of Grid Binnins
  std::vector<std::vector<int>> BinGridMapping;

  /// @brief Initialise Uniform Binning
  void InitUniform(const std::vector<std::vector<double>>& InputEdges) {
    BinEdges = InputEdges;
    Uniform = true;
    AxisNBins.resize(BinEdges.size());

    nBins = 1;
    for(size_t iDim = 0; iDim < BinEdges.size(); iDim++)
    {
      const auto& Edges = BinEdges[iDim];
      if (!std::is_sorted(Edges.begin(), Edges.end())) {
        MACH3LOG_ERROR("VarBins for Dim {} must be in increasing order in sample config, VarBins: [{}]",
                       iDim, fmt::join(Edges, ", "));
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      //Sanity check that some binning has been specified
      if(BinEdges[iDim].size() == 0){
        MACH3LOG_ERROR("No binning specified for Dim {} of sample binning, please add some binning to the sample config", iDim);
        MACH3LOG_ERROR("Please ensure BinEdges are correctly configured for all dimensions");
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      // Set the number of bins for this dimension
      AxisNBins[iDim] = static_cast<int>(BinEdges[iDim].size()) - 1;
      // Update total number of bins
      nBins *= AxisNBins[iDim];

      MACH3LOG_INFO("{}-Dim Binning: [{:.2f}]", iDim, fmt::join(BinEdges[iDim], ", "));
    }
    GlobalOffset = 0;
    /// Lastly prepare special histograms used for event migration
    InitialiseBinMigrationLookUp(static_cast<int>(BinEdges.size()));
  }

  /// @brief Check that non-uniform bin extents do not overlap
  void CheckBinsDoNotOverlap(const std::vector<BinInfo>& TestedBins) const {
    if (TestedBins.empty()) return;

    const size_t ExtentDim = TestedBins[0].Extent.size();

    for (size_t i = 0; i < TestedBins.size(); ++i) {
      for (size_t j = i + 1; j < TestedBins.size(); ++j) {
        bool OverlapsInAllDims = true;

        for (size_t iDim = 0; iDim < ExtentDim; ++iDim) {
          const double a_lo = TestedBins[i].Extent[iDim][0];
          const double a_hi = TestedBins[i].Extent[iDim][1];
          const double b_lo = TestedBins[j].Extent[iDim][0];
          const double b_hi = TestedBins[j].Extent[iDim][1];

          // [low, high) overlap check
          if (!(a_lo < b_hi && b_lo < a_hi)) {
            OverlapsInAllDims = false;
            break;
          }
        }

        if (OverlapsInAllDims) {
          MACH3LOG_ERROR("Overlapping non-uniform bins detected: Bin {} and Bin {}", i, j);
          for (size_t iDim = 0; iDim < ExtentDim; ++iDim) {
            MACH3LOG_ERROR("  Dim {}: Bin {} [{}, {}), Bin {} [{}, {})",
                           iDim,
                           i, TestedBins[i].Extent[iDim][0], TestedBins[i].Extent[iDim][1],
                           j, TestedBins[j].Extent[iDim][0], TestedBins[j].Extent[iDim][1]);
          }
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
    }
  }

  /// @brief Check that non-uniform bins fully cover the bounding box (no gaps)
  ///
  /// The idea is:
  ///  1. Build a Cartesian test grid from all bin edges (TestGridEdges)
  ///  2. For every N-dimensional grid cell, test its midpoint
  ///  3. If any midpoint is not inside at least one bin, a gap exists
  void CheckBinsHaveNoGaps(const std::vector<BinInfo>& TestedBins,
                           const std::vector<double>& MinVal,
                           const std::vector<double>& MaxVal,
                           size_t ValidationBinsPerDim = 100) const {
    bool gap_found = false;
    if (TestedBins.empty()) return;
    const size_t Dim = TestedBins[0].Extent.size();
    if (MinVal.size() != Dim || MaxVal.size() != Dim) {
      MACH3LOG_ERROR("MinVal/MaxVal size does not match dimension of bins");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Build a fine validation grid from the provided min/max
    std::vector<std::vector<double>> TestGridEdges(Dim);
    for (size_t d = 0; d < Dim; ++d) {
      TestGridEdges[d].resize(ValidationBinsPerDim + 1);
      const double width = (MaxVal[d] - MinVal[d]) / static_cast<double>(ValidationBinsPerDim);
      for (size_t i = 0; i <= ValidationBinsPerDim; ++i)
        TestGridEdges[d][i] = MinVal[d] + static_cast<double>(i) * width;
    }
    // Precompute midpoints of each test cell
    std::vector<size_t> indices(Dim, 0);

    std::function<void(size_t)> scan = [&](size_t d) {
      if (gap_found) return; // stop recursion

      // Base case: we have selected a cell in every dimension
      if (d == Dim) {
        std::vector<double> point(Dim);

        // Compute the midpoint of the current N-D grid cell
        for (size_t i = 0; i < Dim; ++i) {
          const double lo = TestGridEdges[i][indices[i]];
          const double hi = TestGridEdges[i][indices[i] + 1];
          point[i] = 0.5 * (lo + hi);
        }

        // Check coverage
        for (const auto& bin : TestedBins) {
          if (bin.IsEventInside(point)) {
            return; // covered
          }
        }

        // Not covered by any bin → gap
        MACH3LOG_WARN("Gap detected in non-uniform binning at point [{:.2f}]", fmt::join(point, ", "));
        gap_found = true;
        return;
      }

      for (size_t i = 0; i + 1 < TestGridEdges[d].size(); ++i) {
        indices[d] = i;
        scan(d + 1);
      }
    };
    // Start recursion at dimension 0
    scan(0);
  }

  /// @brief Initialise Non-Uniform Binning
  void InitialiseGridMapping() {
    const size_t Dim = BinEdges.size();

    // Compute total number of "mega bins"
    int NGridBins = 1;
    for (size_t d = 0; d < Dim; ++d) {
      NGridBins *= AxisNBins[d];
    }
    BinGridMapping.resize(NGridBins);

    // Helper: convert linear index to multi-dimensional index
    std::vector<int> MultiIndex(Dim, 0);

    std::function<void(size_t, int)> scan;
    scan = [&](size_t d, int LinearIndex) {
      if (d == Dim) {
        // Compute the edges of the current mega-bin
        std::vector<std::array<double, 2>> CellEdges(Dim);
        for (size_t i = 0; i < Dim; ++i) {
          CellEdges[i][0] = BinEdges[i][MultiIndex[i]];
          CellEdges[i][1] = BinEdges[i][MultiIndex[i] + 1];
        }

        // Check overlap with all non-uniform bins
        for (size_t iBin = 0; iBin < Bins.size(); ++iBin) {
          auto& bin = Bins[iBin];
          bool overlap = true;
          for (size_t i = 0; i < Dim; ++i) {
            const double a_lo = bin.Extent[i][0];
            const double a_hi = bin.Extent[i][1];
            const double b_lo = CellEdges[i][0];
            const double b_hi = CellEdges[i][1];
            if (!(a_hi > b_lo && a_lo < b_hi)) { // overlap condition
              overlap = false;
              break;
            }
          }
          if (overlap) {
            BinGridMapping[LinearIndex].push_back(static_cast<int>(iBin));

            // Debug statement: show both small bin extent AND mega-bin edges
            std::vector<std::string> bin_extent_str(Dim);
            std::vector<std::string> mega_edges_str(Dim);

            for (size_t i = 0; i < Dim; ++i) {
              bin_extent_str[i] = fmt::format("[{:.3f}, {:.3f}]", bin.Extent[i][0], bin.Extent[i][1]);
              mega_edges_str[i]  = fmt::format("[{:.3f}, {:.3f}]", CellEdges[i][0], CellEdges[i][1]);
            }

            MACH3LOG_DEBUG("MegaBin {} (multi-index [{}], edges {}) assigned Bin {} with extents {}",
                           LinearIndex, fmt::join(MultiIndex, ","), fmt::join(mega_edges_str, ", "),
                           iBin, fmt::join(bin_extent_str, ", "));
          }
        }
        return;
      }

      // Loop over all bins along this dimension
      for (int i = 0; i < AxisNBins[d]; ++i) {
        MultiIndex[d] = i;
        int NewLinearIndex = LinearIndex;
        if (d > 0) {
          int stride = 1;
          for (size_t s = 0; s < d; ++s) stride *= AxisNBins[s];
          NewLinearIndex += i * stride;
        } else {
          NewLinearIndex = i;
        }
        scan(d + 1, NewLinearIndex);
      }
    };
    // Start the recursive scan over all dimensions, beginning at dimension 0 with linear index 0
    scan(0, 0);
  }

  /// @brief Initialise Non-Uniform Binning
  void InitNonUniform(const std::vector<std::vector<std::vector<double>>>& InputBins) {
    Uniform = false;
    nBins = 0;
    Bins.resize(InputBins.size());

    size_t ExtentDim = InputBins[0].size();
    if (ExtentDim == 1) {
      MACH3LOG_ERROR("Trying to initialise Non-Uniform binning for single dimension, this is silly...");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    for(size_t iBin = 0; iBin < InputBins.size(); iBin++) {
      const auto& NewExtent = InputBins[iBin];
      if (NewExtent.size() != ExtentDim) {
        MACH3LOG_ERROR("Dimension of Bin {} is {}, while others have {}", iBin, NewExtent.size(), ExtentDim);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      BinInfo NewBin;
      for (const auto& extent : NewExtent) {
        if (extent.size() != 2) {
          MACH3LOG_ERROR("Extent size is not 2 for Bin {}", iBin);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        NewBin.Extent.push_back({extent[0], extent[1]});
        MACH3LOG_DEBUG("Adding extent for Bin {} Dim {}: [{:.2f}, {:.2f}]",
                       iBin, NewBin.Extent.size()-1, NewBin.Extent.back()[0], NewBin.Extent.back()[1]);
      }
      Bins[iBin] = std::move(NewBin);
      nBins++;
    }
    // Ensure we do not have weird overlaps
    CheckBinsDoNotOverlap(Bins);
    /// @todo KS: Now we create "Large Bins" automatically, in future we can expand to add more user control
    constexpr int BinsPerDimension = 10;
    // Now we create huge map, which will allow to easily find non uniform bins
    BinEdges.resize(ExtentDim);
    AxisNBins.resize(ExtentDim);
    std::vector<double> MinVal(ExtentDim), MaxVal(ExtentDim);
    for (size_t iDim = 0; iDim < ExtentDim; iDim++) {
      MinVal[iDim] = std::numeric_limits<double>::max();
      MaxVal[iDim] = std::numeric_limits<double>::lowest();

      // Find min and max for this dimension
      for (const auto& bin : Bins) {
        MinVal[iDim] = std::min(MinVal[iDim], bin.Extent[iDim][0]);
        MaxVal[iDim] = std::max(MaxVal[iDim], bin.Extent[iDim][1]);
      }

      MACH3LOG_DEBUG("Mapping binning: Dim {} Min = {:.2f}, Max = {:.2f}", iDim, MinVal[iDim], MaxVal[iDim]);
      BinEdges[iDim].resize(BinsPerDimension + 1);
      double BinWidth = (MaxVal[iDim] - MinVal[iDim]) / static_cast<double>(BinsPerDimension);
      for (size_t iEdge = 0; iEdge <= BinsPerDimension; iEdge++) {
        BinEdges[iDim][iEdge] = MinVal[iDim] + static_cast<double>(iEdge) * BinWidth;
      }
      AxisNBins[iDim] = BinsPerDimension;
      MACH3LOG_DEBUG("Mapping binning: Dim {} BinEdges = [{:.2f}]", iDim, fmt::join(BinEdges[iDim], ", "));
    }
    CheckBinsHaveNoGaps(Bins, MinVal, MaxVal, 200);
    GlobalOffset = 0;
    /// prepare special histograms used for event migration
    InitialiseBinMigrationLookUp(static_cast<int>(BinEdges.size()));
    /// Lastly GridMap
    InitialiseGridMapping();
  }

  /// @brief Get linear bin index from ND bin indices with additional checks
  /// @param BinIndices Vector of bin indices along each dimension
  /// @warning this performs additional checks so do not use in parts of code used during fit
  int GetBinSafe(const std::vector<int>& BinIndices) const {
    for(int iDim = 0; iDim < static_cast<int>(BinIndices.size()); iDim++){
      if (BinIndices[iDim] < 0 || BinIndices[iDim] >= AxisNBins[iDim]) {
        MACH3LOG_ERROR("{}: Bin indices out of range: Dim = {}, Bin={}, max Ndim Bin={}",
                       __func__, iDim, BinIndices[iDim], AxisNBins[iDim]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    return GetBin(BinIndices);
  }

  /// @brief Convert N-dimensional bin indices to a linear bin index.
  /// @param BinIndices Vector of bin indices along each dimension
  /// @details
  /// Mapping follows row-major order:
  ///  - 1D: xBin
  ///  - 2D: yBin * nXBins + xBin
  ///  - 3D: zBin * nXBins * nYBins + yBin * nXBins + xBin
  int GetBin(const std::vector<int>& BinIndices) const {
    int BinNumber = 0;
    for(size_t i = 0; i < BinIndices.size(); ++i) {
      BinNumber += BinIndices[i]*Strides[i];
    }
    return BinNumber;
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
              const int N_Bins,
              const std::vector<double>& Bin_Edges,
              const std::vector<BinShiftLookup>& Bin_Lookup) const {
    //DB Check to see if momentum shift has moved bins
    //DB - First , check to see if the event is outside of the binning range and skip event if it is
    if (KinVar < Bin_Edges[0] || KinVar >= Bin_Edges[N_Bins]) {
      return M3::UnderOverFlowBin;
    }
    // KS: If NomBin is UnderOverFlowBin we must do binary search :(
    if(NomBin > M3::UnderOverFlowBin) {
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
    }
    //DB - If we end up in this loop, the event has been shifted outside of its nominal bin, but is still within the allowed binning range
    // KS: Perform binary search to find correct bin. We already checked if isn't outside of bounds
    return static_cast<int>(std::distance(Bin_Edges.begin(), std::upper_bound(Bin_Edges.begin(), Bin_Edges.end(), KinVar)) - 1);
  }

  /// @brief Initializes lookup arrays for efficient bin migration in a single dimension.
  /// @param Bin_Lookup Reference to the BinShiftLookup struct to be initialized.
  /// @param Bin_Edges Vector of bin edges defining the bin boundaries.
  /// @param TotBins Number of bins in the dimension.
  void InitialiseLookUpSingleDimension(std::vector<BinShiftLookup>& Bin_Lookup, const std::vector<double>& Bin_Edges, const int TotBins) {
    Bin_Lookup.resize(TotBins);
    //Set rw_pdf_bin and upper_binedge and lower_binedge for each skmc_base
    for(int bin_i = 0; bin_i < TotBins; bin_i++){
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
      stride *= AxisNBins[i];
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
inline int GetSampleFromGlobalBin(const std::vector<SampleBinningInfo>& BinningInfo, const int GlobalBin) {
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
                                    const int GlobalBin) {
  for (size_t iSample = 0; iSample < BinningInfo.size(); ++iSample) {
    const SampleBinningInfo& info = BinningInfo[iSample];

    if (GlobalBin >= info.GlobalOffset &&
      GlobalBin <  info.GlobalOffset + info.nBins)
    {
      return GlobalBin - info.GlobalOffset;
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
