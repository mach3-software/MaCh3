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
void BinningHandler::SetupSampleBinning(const YAML::Node& Settings, SampleInfo& SingleSample) {
// ************************************************
  MACH3LOG_INFO("Setting up Sample Binning");

  //Binning
  SingleSample.nDimensions = 0;
  SingleSample.XVarStr = GetFromManager(Settings["XVarStr"], std::string(""));
  auto X_BinEdges = GetFromManager(Settings["XVarBins"], std::vector<double>());
  const auto& edgesx = X_BinEdges;
  if (!std::is_sorted(edgesx.begin(), edgesx.end())) {
    MACH3LOG_ERROR("XVarBins must be in increasing order in sample config {}\n  XVarBins: [{}]",
                   SingleSample.SampleTitle, fmt::join(edgesx, ", "));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if(SingleSample.XVarStr.length() > 0){
    SingleSample.nDimensions++;
  } else{
    MACH3LOG_ERROR("Please specify an X-variable string in sample config");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  SingleSample.YVarStr = GetFromManager(Settings["YVarStr"], std::string(""));
  auto Y_BinEdges = GetFromManager(Settings["YVarBins"], std::vector<double>());
  const auto& edgesy = Y_BinEdges;
  if (!std::is_sorted(edgesy.begin(), edgesy.end())) {
    MACH3LOG_ERROR("Y_BinEdges must be in increasing order in sample config {}\n  Y_BinEdges: [{}]",
                   SingleSample.SampleTitle, fmt::join(edgesy, ", "));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if(SingleSample.YVarStr.length() > 0){
    if(SingleSample.XVarStr.length() == 0){
      MACH3LOG_ERROR("Please specify an X-variable string in sample config. I won't work only with a Y-variable");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    SingleSample.nDimensions++;
  }

  if(SingleSample.nDimensions == 0){
    MACH3LOG_ERROR("Error setting up the sample binning");
    MACH3LOG_ERROR("Number of dimensions is {}", SingleSample.nDimensions);
    MACH3LOG_ERROR("Check that an XVarStr has been given in the sample config");
    throw MaCh3Exception(__FILE__, __LINE__);
  } else{
    MACH3LOG_INFO("Found {} dimensions for sample binning", SingleSample.nDimensions);
  }

  //Check whether you are setting up 1D or 2D binning
  if(SingleSample.nDimensions == 1){
    MACH3LOG_INFO("Setting up {}D binning with {}", SingleSample.nDimensions, SingleSample.XVarStr);
    Y_BinEdges = {-1e8, 1e8};
  } else if(SingleSample.nDimensions == 2){
    MACH3LOG_INFO("Setting up {}D binning with {} and {}", SingleSample.nDimensions, SingleSample.XVarStr, SingleSample.YVarStr);
  } else{
    MACH3LOG_ERROR("Number of dimensions is not 1 or 2, this is unsupported at the moment");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

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
    MACH3LOG_ERROR("Please ensure XVarBins and/or YVarStr are correctly configured");
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
  const int GlobalBin = SampleBinning[Sample].GetBinSafe(xBin, yBin);
  return GlobalBin;
}

// ************************************************
int BinningHandler::GetGlobalBinSafe(const int Sample, const int xBin, const int yBin) const {
// ************************************************
  const int GlobalBin = SampleBinning[Sample].GetBinSafe(xBin, yBin) + static_cast<int>(SampleBinning[Sample].GlobalOffset);
  return GlobalBin;
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
