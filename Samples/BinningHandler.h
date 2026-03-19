#pragma once
#include <vector>
#include <string>

#include "Samples/SampleStructs.h"
#include "Samples/FarDetectorCoreInfoStruct.h"
#include "Manager/Manager.h"

// ***************************
/// @brief KS: Class handling binning for multiple samples
/// @details
/// ## Introduction
/// Each sample can define its own binning in an arbitrary number of dimensions.
/// Internally, every sample's multi-dimensional binning is linearised into a
/// single 1D array. All samples are then concatenated into one global bin index
/// space, allowing the entire analysis to be treated as a single large vector
/// of bins.
///
/// The concept of a "global bin" refers to the position of a bin in this
/// linearised, analysis-wide bin index space. Local (sample) bins are always
/// enumerated starting from zero, while global bins span all samples
/// consecutively.
///
/// Example layout of global bins with offsets:
/// @code
///   Sample 0 (GlobalOffset = 0,  nBins = 4):
///     Local bins:   [0] [1] [2] [3]
///     Global bins:  [0] [1] [2] [3]
///
///   Sample 1 (GlobalOffset = 4,  nBins = 3):
///     Local bins:   [0] [1] [2]
///     Global bins:  [4] [5] [6]
///
///   Sample 2 (GlobalOffset = 7,  nBins = 2):
///     Local bins:   [0] [1]
///     Global bins:  [7] [8]
///
///   Global bin index space:
///   ------------------------------------------------
///   | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
///   ------------------------------------------------
/// @endcode
///
/// ## Uniform and Non-Uniform Binning Scheme
/// MaCh3 supports Uniform and Non-Uniform binning scheme
///
/// In the non-uniform scheme, bin sizes may vary along each dimension,
/// but all bins are required to be axis-aligned hyper-rectangles.
/// Arbitrary or irregular bin shapes are not supported like banana-shape.
///
/// Example of Uniform
/// @code
///
/// +--------+--------+--------+--------+
/// | Bin 0  | Bin 1  | Bin 2  | Bin 3  |
/// | (0,0)  | (1,0)  | (2,0)  | (3,0)  |
/// +--------+--------+--------+--------+
/// | Bin 4  | Bin 5  | Bin 6  | Bin 7  |
/// | (0,1)  | (1,1)  | (2,1)  | (3,1)  |
/// +--------+--------+--------+--------+
/// | Bin 8  | Bin 9  | Bin 10 | Bin 11 |
/// | (0,2)  | (1,2)  | (2,2)  | (3,2)  |
/// +--------+--------+--------+--------+
///
/// @endcode
///
/// Example of Non-Uniform
///
/// @code
///
/// +--------+------------+-------+---------------------+
/// | Bin 0  |    Bin 1   | Bin 2 |         Bin 3       |
/// |        |            |       |                     |
/// +--------+------------+-------+---------------------+
/// | Bin 4  |            |       |                     |
/// |        |            | Bin 6 |                     |
/// +--------+   Bin 5    +-------+                     |
/// | Bin 7  |            |       |         Bin 9       |
/// |        |            | Bin 8 |                     |
/// +--------+------------+-------+---------------------+
///
/// @endcode
///
///
/// ## Bin Finding Algorithm
/// Since MaCh3 supports event migration bin finding algorithm must be fast
/// to efficiently be able find bin during running fit.
/// MaCh3 is caching nominal bin with idea that during fit migration should be around this nominal bin.
/// Thus MaCh3 first checks if after shift event falls into Nom-bin and later adjacent. If not backs to binary search.
///
/// ### Uniform
/// In case of uniform binning above algorithm is easy to test as one performs it for every dimension independently
/// i.e. find X-bin, then Y etc. After which can find bin in flattened 1D space.
///
///
/// ### Non-Uniform
/// Internally, non-uniform binning is implemented using two levels:
///
/// 1. **MegaBins (mapping bins)**
///    These form a coarse, *uniform* grid that spans the full phase space.
///    Each MegaBin acts as a container for one or more non-uniform bins.
///
/// 2. **Non-Uniform bins**
///    The actual analysis bins, defined as hyper-rectangles with arbitrary
///    extents inside a MegaBin.
///
/// #### Bin finding procedure
/// For a given event, the bin-finding algorithm proceeds as follows:
///
/// 1. Locate the MegaBin using the same fast per-dimension logic as in
///    uniform binning.
/// 2. Once the MegaBin is identified, loop over all non-uniform bins
///    associated with that MegaBin.
/// 3. The first bin whose extents fully contain the event is selected.
/// 4. If no bin matches, the event is assigned to the under/overflow bin.
///
///
///
/// @author Kamil Skwarczynski
/// @author Dan Barrow
class BinningHandler {
// ***************************
 public:
  /// @brief Constructor
  BinningHandler();

  /// @brief destructor
  virtual ~BinningHandler() {};

  /// @brief Find Global bin including
  /// @param iSample index of a given sample
  /// @param KinVar Vector of pointers to kinematic variable like Erec
  /// @param NomBin Vector of nominal bin indices for this event, one per dimension.
  int FindGlobalBin(const int iSample,
                    const std::vector<const double*>& KinVar,
                    const std::vector<int>& NomBin) const;

  /// @brief Find the nominal bin for a given variable in a given sample and dimension
  /// @param iSample Sample index
  /// @param iDim Dimension index (0 = X, 1 = Y, ...)
  /// @param Var Kinematic variable value
  int FindNominalBin(const int iSample, const int iDim, const double Var) const;

  /// @brief Get gloabl bin based on sample, and dimension of each sample with additional checks
  /// @param iSample index of a given sample
  /// @param Bins Vector of bin indices along each dimension
  int GetGlobalBinSafe(const int iSample, const std::vector<int>& Bins) const;
  /// @brief Get gloabl bin based on sample, and dimension of each sample without any safety checks
  /// @param iSample index of a given sample
  /// @param Bins Vector of bin indices along each dimension
  int GetBinSafe(const int iSample, const std::vector<int>& Bins) const;
  /// @brief Get total number of bins over all samples/kinematic bins etc
  int GetNBins() const {return TotalNumberOfBins;};
  /// @brief Get total number of bins over for a given sample
  int GetNBins(const int iSample) const {return static_cast<int>(SampleBinning[iSample].nBins);};
  /// @brief Get fancy name for a given bin, to help match it with global properties
  /// @param GlobalBin Global Bin integrated over all samples
  std::string GetBinName(const int GlobalBin) const;
  /// @brief Get fancy name for a given bin, to help match it with global properties
  /// @param iSample index of a given sample
  /// @param SampleBin Global Bin for a given sample
  std::string GetBinName(const int iSample, const int SampleBin) const;
    /// @brief Get fancy name for a given bin, to help match it with global properties
  /// @param iSample index of a given sample
  /// @param Bins Vector of bin indices along each dimension
  std::string GetBinName(const int iSample, const std::vector<int>& Bins) const;

  /// @brief Get N-dim bin edges for a given sample
  /// @param iSample index of a given sample
  /// @param iDim dimension for which we extract bin edges
  std::vector<double> GetBinEdges(const int iSample, const int iDim) const {return SampleBinning[iSample].BinEdges.at(iDim);};

  /// @brief Get Number of N-axis bins for a given sample
  /// @param iSample index of a given sample
  /// @param iDim dimension for which we extract number of bins
  int GetNAxisBins(const int iSample, const int iDim) const;
  /// @brief Tells whether given sample is using unform binning
  /// @param iSample index of a given sample
  bool IsUniform(const int iSample) const;
  /// @brief Get bin number corresponding to where given sample starts
  /// @param iSample index of a given sample
  int GetSampleStartBin(const int iSample) const;
  /// @brief Get bin number corresponding to where given sample ends
  /// @param iSample index of a given sample
  int GetSampleEndBin(const int iSample) const;
  /// @brief Return NonUnifomr bins to for example check extent etc
  const std::vector<BinInfo> GetNonUniformBins(const int iSample) const;
  /// @brief Sets the GlobalOffset for each SampleBinningInfo to enable linearization of multiple 2D binning samples.
  void SetGlobalBinNumbers();

  /// @brief Function to setup the binning of your sample histograms and the underlying
  /// arrays that get handled in fillArray() and fillArray_MP().
  void SetupSampleBinning(const YAML::Node& Settings, SampleInfo& SingleSample);

 private:
  /// Total number of bins
  int TotalNumberOfBins;
  /// Binning info for individual sample
  std::vector<SampleBinningInfo> SampleBinning;
};
