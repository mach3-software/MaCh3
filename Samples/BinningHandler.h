#pragma once
#include <vector>
#include <string>

#include "Samples/SampleStructs.h"

// ***************************
/// @brief KS: Class handling binning
/// @author Kamil Skwarczynski
/// @author Dan Barrow
class BinningHandler {
// ***************************
 public:
  /// @brief Constructor
  BinningHandler();

  /// @brief destructor
  virtual ~BinningHandler() {};

  /// Find Bin for
  /// @todo add Y binning migration
  int FindGlobalBin(const int NomSample, const double XVar, const int NomXBin, const int NomYBin) const;


  int GetGlobalBinSafe(const int Sample, const int xBin, const int yBin) const;
  int GetBinSafe(const int Sample, const int xBin, const int yBin) const;

  int GetNBins() const {return TotalNumberOfBins;};

  std::vector<double> GetXBinEdges(const int Samples) const {return SampleBinning[Samples].XBinEdges;};
  std::vector<double> GetYBinEdges(const int Samples) const {return SampleBinning[Samples].YBinEdges;};

  int GetNXBins(const int Samples) const {return static_cast<int>(SampleBinning[Samples].nXBins);};
  int GetNYBins(const int Samples) const {return static_cast<int>(SampleBinning[Samples].nYBins);};

  int GetSampleStartBin(const int Samples) const;
  int GetSampleEndBin(const int Samples) const;

  /// @brief Sets the GlobalOffset for each SampleBinningInfo to enable linearization of multiple 2D binning samples.
  void SetGlobalBinNumbers();


  /// @brief Function to setup the binning of your sample histograms and the underlying
  /// arrays that get handled in fillArray() and fillArray_MP().
  /// The Binning.XBinEdges are filled in the daughter class from the sample config file.
  /// This "passing" can be removed.
  void SetupSampleBinning(const std::vector<double>& X_BinEdges, const std::vector<double>& Y_BinEdges);

 private:
  int TotalNumberOfBins;
  std::vector<SampleBinningInfo> SampleBinning;
};
