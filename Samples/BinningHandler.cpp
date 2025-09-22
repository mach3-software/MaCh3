#include "Samples/BinningHandler.h"


// ************************************************
BinningHandler::BinningHandler() {
// ************************************************

}

// ************************************************
// Function to setup the binning of your sample histograms and the underlying
// arrays that get handled in fillArray() and fillArray_MP().
// The Binning.XBinEdges are filled in the daughter class from the sample config file.
// This "passing" can be removed.
void BinningHandler::SetupSampleBinning(const std::vector<double>& X_BinEdges, const std::vector<double>& Y_BinEdges) {
// ************************************************
  MACH3LOG_INFO("Setting up Sample Binning");

  SampleBinningInfo SingleBinning;

  SingleBinning.XBinEdges = X_BinEdges;
  SingleBinning.YBinEdges = Y_BinEdges;

  //A string to store the binning for a nice print out
  std::string XBinEdgesStr = "";
  std::string YBinEdgesStr = "";

  for(auto XBinEdge : SingleBinning.XBinEdges){
    XBinEdgesStr += std::to_string(XBinEdge);
    XBinEdgesStr += ", ";
  }
  MACH3LOG_INFO("XBinning:");
  MACH3LOG_INFO("{}", XBinEdgesStr);

  //And now the YBin Edges
  for(auto YBinEdge : SingleBinning.YBinEdges){
    YBinEdgesStr += std::to_string(YBinEdge);
    YBinEdgesStr += ", ";
  }
  MACH3LOG_INFO("YBinning:");
  MACH3LOG_INFO("{}", YBinEdgesStr);


  //Sanity check that some binning has been specified
  if(SingleBinning.XBinEdges.size() == 0 && SingleBinning.YBinEdges.size() == 0){
    MACH3LOG_ERROR("No binning specified for either X or Y of sample binning, please add some binning to the sample config");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //Set the number of X and Y bins now
  SingleBinning.nXBins = SingleBinning.XBinEdges.size() - 1;
  SingleBinning.nYBins = SingleBinning.YBinEdges.size() - 1;

  // Set total number of bins
  SingleBinning.nBins = SingleBinning.nXBins * SingleBinning.nYBins;
  SingleBinning.GlobalOffset = 0;
  /// Lastly prepare special histograms used for event migration
  SingleBinning.InitialiseBinMigrationLookUp();

  SampleBinning.emplace_back(SingleBinning);

  // now setup global numbering
  SetGlobalBinNumbers();
}

// ************************************************
int BinningHandler::FindGlobalBin(const int NomSample, const double XVar, const int NomXBin, const int NomYBin) const {
// ************************************************
  //DB Find the relevant bin in the PDF for each event
  const int XBinToFill = SampleBinning[NomSample].FindXBin(XVar, NomXBin);
  const int YBinToFill = NomYBin;

  //DB Fill relevant part of thread array
  if (XBinToFill != -1 && YBinToFill != -1) {
    int GlobalBin = SampleBinning[NomSample].GetBin(XBinToFill, YBinToFill);
    GlobalBin += static_cast<int>(SampleBinning[NomSample].GlobalOffset);
    return GlobalBin;
  } else {
    return -1;
  }
}

// ************************************************
int BinningHandler::GetBinSafe(const int Sample, const int xBin, const int yBin) const {
  // ************************************************
  const int GlobalbBin = SampleBinning[Sample].GetBinSafe(xBin, yBin);
  return GlobalbBin;
}

// ************************************************
int BinningHandler::GetGlobalBinSafe(const int Sample, const int xBin, const int yBin) const {
// ************************************************
  const int GlobalbBin = SampleBinning[Sample].GetBinSafe(xBin, yBin) + static_cast<int>(SampleBinning[Sample].GlobalOffset);
  return GlobalbBin;
}

// ************************************************
int BinningHandler::GetSampleStartBin(const int Sample) const {
// ************************************************
 return static_cast<int>(SampleBinning[Sample].GlobalOffset);
}

// ************************************************
int BinningHandler::GetSampleEndBin(const int Sample) const {
// ************************************************
  if (Sample == static_cast<int>(SampleBinning.size()) - 1) {
    return GetNBins();
  } else {
    return static_cast<int>(SampleBinning[Sample+1].GlobalOffset);
  }
}

// ************************************************
/// @brief Sets the GlobalOffset for each SampleBinningInfo to enable linearization of multiple 2D binning samples.
void BinningHandler::SetGlobalBinNumbers() {
// ************************************************
  if (SampleBinning.empty()) {
    MACH3LOG_ERROR("No binning samples provided.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  size_t GlobalOffsetCounter = 0;
  for(size_t iSample = 0; iSample < SampleBinning.size(); iSample++){
    SampleBinning[iSample].GlobalOffset = GlobalOffsetCounter;
    GlobalOffsetCounter += SampleBinning[iSample].nBins;
  }
  // lastly modify total number of bins
  TotalNumberOfBins = static_cast<int>(GlobalOffsetCounter);
}
