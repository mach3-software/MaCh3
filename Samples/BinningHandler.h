#pragma once
#include <vector>
#include <string>

#include "Samples/SampleStructs.h"

// ***************************
/// @brief KS: Class handling binning for multiple samples
/// @details Sample ix treated as additional dimension so everything is treated as one huge 1D array
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
  /// @todo add Y binning migration
  /// @param iSample index of a given sample
  int FindGlobalBin(const int iSample, const double XVar, const int NomXBin, const int NomYBin) const;

  /// @brief Get gloabl bin based on sample, and dimension of each sample with additional checks
  /// @param iSample index of a given sample
  int GetGlobalBinSafe(const int iSample, const int xBin, const int yBin) const;
  /// @brief Get gloabl bin based on sample, and dimension of each sample without any safety checks
  /// @param iSample index of a given sample
  int GetBinSafe(const int iSample, const int xBin, const int yBin) const;
  /// @brief Get total number of bins over all samples/kinematic bins etc
  int GetNBins() const {return TotalNumberOfBins;};

  /// @brief Get X bin edges for a given sample
  /// @param iSample index of a given sample
  std::vector<double> GetXBinEdges(const int iSample) const {return SampleBinning[iSample].XBinEdges;};
  /// @brief Get Y bin edges for a given sample
  /// @param iSample index of a given sample
  std::vector<double> GetYBinEdges(const int iSample) const {return SampleBinning[iSample].YBinEdges;};

  /// @brief Get Number of X bins for a given sample
  /// @param iSample index of a given sample
  int GetNXBins(const int iSample) const {return static_cast<int>(SampleBinning[iSample].nXBins);};
  /// @brief Get Number of Y bins for a given sample
  /// @param iSample index of a given sample
  int GetNYBins(const int iSample) const {return static_cast<int>(SampleBinning[iSample].nYBins);};

  /// @brief Get bin number corresponding to where given sample starts
  /// @param iSample index of a given sample
  int GetSampleStartBin(const int iSample) const;
  /// @brief Get bin number corresponding to where given sample ends
  /// @param iSample index of a given sample
  int GetSampleEndBin(const int iSample) const;

  /// @brief Sets the GlobalOffset for each SampleBinningInfo to enable linearization of multiple 2D binning samples.
  void SetGlobalBinNumbers();

  /// @brief Function to setup the binning of your sample histograms and the underlying
  /// arrays that get handled in fillArray() and fillArray_MP().
  /// The Binning.XBinEdges are filled in the daughter class from the sample config file.
  /// This "passing" can be removed.
  void SetupSampleBinning(const std::vector<double>& X_BinEdges, const std::vector<double>& Y_BinEdges);

 private:
  /// Total number of bins
  int TotalNumberOfBins;
  /// Binning info for individual sample
  std::vector<SampleBinningInfo> SampleBinning;
};
