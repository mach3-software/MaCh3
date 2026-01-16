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
  /// @warning for now we hardcode to 2D...
  SingleSample.VarStr.resize(2);
  SingleSample.VarStr[0] = GetFromManager(Settings["XVarStr"], std::string(""));
  auto X_BinEdges = GetFromManager(Settings["XVarBins"], std::vector<double>());
  const auto& edgesx = X_BinEdges;
  if (!std::is_sorted(edgesx.begin(), edgesx.end())) {
    MACH3LOG_ERROR("XVarBins must be in increasing order in sample config {}\n  XVarBins: [{}]",
                   SingleSample.SampleTitle, fmt::join(edgesx, ", "));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if(SingleSample.VarStr[0].length() > 0){
    SingleSample.nDimensions++;
  } else{
    MACH3LOG_ERROR("Please specify an X-variable string in sample config");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  SingleSample.VarStr[1] = GetFromManager(Settings["YVarStr"], std::string(""));
  auto Y_BinEdges = GetFromManager(Settings["YVarBins"], std::vector<double>());
  const auto& edgesy = Y_BinEdges;
  if (!std::is_sorted(edgesy.begin(), edgesy.end())) {
    MACH3LOG_ERROR("Y_BinEdges must be in increasing order in sample config {}\n  Y_BinEdges: [{}]",
                   SingleSample.SampleTitle, fmt::join(edgesy, ", "));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if(SingleSample.VarStr[1].length() > 0){
    if(SingleSample.VarStr[0].length() == 0){
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
    MACH3LOG_INFO("Setting up {}D binning with {}", SingleSample.nDimensions, SingleSample.VarStr[0]);
    Y_BinEdges = {-1e8, 1e8};
  } else if(SingleSample.nDimensions == 2){
    MACH3LOG_INFO("Setting up {}D binning with {} and {}", SingleSample.nDimensions, SingleSample.VarStr[0], SingleSample.VarStr[1]);
  } else{
    MACH3LOG_ERROR("Number of dimensions is not 1 or 2, this is unsupported at the moment");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  SampleBinningInfo SingleBinning;
  /// @warning for now we hardcode to 2D...
  SingleBinning.BinEdges.resize(2);
  SingleBinning.AxisNBins.resize(2);
  SingleBinning.BinEdges[0] = X_BinEdges;
  SingleBinning.BinEdges[1] = Y_BinEdges;

  //A string to store the binning for a nice print out
  std::string XBinEdgesStr = "";
  std::string YBinEdgesStr = "";

  for(auto XBinEdge : SingleBinning.BinEdges[0]){
    XBinEdgesStr += std::to_string(XBinEdge);
    XBinEdgesStr += ", ";
  }
  MACH3LOG_INFO("XBinning:");
  MACH3LOG_INFO("{}", XBinEdgesStr);

  //And now the YBin Edges
  for(auto YBinEdge : SingleBinning.BinEdges[1]){
    YBinEdgesStr += std::to_string(YBinEdge);
    YBinEdgesStr += ", ";
  }
  MACH3LOG_INFO("YBinning:");
  MACH3LOG_INFO("{}", YBinEdgesStr);


  //Sanity check that some binning has been specified
  if(SingleBinning.BinEdges[0].size() == 0 && SingleBinning.BinEdges[1].size() == 0){
    MACH3LOG_ERROR("No binning specified for either X or Y of sample binning, please add some binning to the sample config");
    MACH3LOG_ERROR("Please ensure XVarBins and/or YVarStr are correctly configured");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //Set the number of X and Y bins now
  SingleBinning.AxisNBins[0] = SingleBinning.BinEdges[0].size() - 1;
  SingleBinning.AxisNBins[1] = SingleBinning.BinEdges[1].size() - 1;

  // Set total number of bins
  SingleBinning.nBins = SingleBinning.AxisNBins[0] * SingleBinning.AxisNBins[1];
  SingleBinning.GlobalOffset = 0;
  /// Lastly prepare special histograms used for event migration
  SingleBinning.InitialiseBinMigrationLookUp(SingleSample.nDimensions);

  SampleBinning.emplace_back(SingleBinning);

  // now setup global numbering
  SetGlobalBinNumbers();
}

// ************************************************
int BinningHandler::FindGlobalBin(const int NomSample,
                                  const std::vector<const double*>& KinVar,
                                  const std::vector<int>& NomBin) const {
// ************************************************
  //DB Find the relevant bin in the PDF for each event
  const int Dim = static_cast<int>(KinVar.size());
  const SampleBinningInfo& _restrict_ SB = SampleBinning[NomSample];
  int GlobalBin = 0;

  for(int i = 0; i < Dim; ++i) {
    const double Var = *KinVar[i];
    const int Bin = SB.FindBin(i, Var, NomBin[i]);
    // KS: If we are outside of range in only one dimension this mean out of bounds, we can simply quickly finish
    if(Bin < 0) return M3::UnderOverFlowBin;
    // KS: inline GetBin computation to avoid any allocation
    GlobalBin += Bin * SB.Strides[i];
  }

  GlobalBin += static_cast<int>(SB.GlobalOffset);
  return GlobalBin;
}

// ************************************************
int BinningHandler::FindNominalBin(const int iSample,
                   const int iDim,
                   const double Var) const {
// ************************************************
  const SampleBinningInfo& info = SampleBinning[iSample];

  const auto& edges = info.BinEdges[iDim];

  // Outside binning range
  if (Var < edges.front() || Var >= edges.back()) {
    return M3::UnderOverFlowBin;
  }
  return static_cast<int>(std::distance(edges.begin(), std::upper_bound(edges.begin(), edges.end(), Var)) - 1);
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

// ************************************************
// Get fancy name for a given bin, to help match it with global properties
std::string BinningHandler::GetBinName(const int iSample, const int xBin, const int yBin) const {
// ************************************************
  return GetBinName(iSample, GetBinSafe(iSample, xBin, yBin));
}

// ************************************************
// Get fancy name for a given bin, to help match it with global properties
std::string BinningHandler::GetBinName(const int iSample, const int GlobSampleBin) const {
// ************************************************
  const auto& Binning = SampleBinning[iSample];

  // Safety checks
  if (GlobSampleBin < 0 || GlobSampleBin >= static_cast<int>(Binning.nBins)) {
    MACH3LOG_ERROR("Requested bin {} is out of range for sample {}", GlobSampleBin, iSample);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  int Dim = static_cast<int>(Binning.Strides.size());
  std::vector<int> Bins(Dim, 0);
  int Remaining = GlobSampleBin;

  // Convert the flat/global bin index into per-dimension indices
  // Dim0 is the fastest-changing axis, Dim1 the next, etc.
  //
  // For example (2D):
  //   x = bin % Nx
  //   y = bin / Nx
  //
  // For 3D:
  //   x = bin % Nx
  //   y = (bin / Nx) % Ny
  //   z = bin / (Nx * Ny)
  for (int i = 0; i < Dim; ++i) {
    const int nBinsDim = static_cast<int>(Binning.BinEdges[i].size()) - 1;
    Bins[i] = Remaining % nBinsDim;
    Remaining /= nBinsDim;
  }

  std::string BinName;
  for (int i = 0; i < Dim; ++i) {
    if (i > 0) BinName += ", ";
    const double min = Binning.BinEdges[i].at(Bins[i]);
    const double max = Binning.BinEdges[i].at(Bins[i] + 1);
    BinName += fmt::format("Dim{} ({:g}, {:g})", i, min, max);
  }

  return BinName;
}

// ************************************************
std::string BinningHandler::GetBinName(const int GlobalBin) const {
// ************************************************
  int SampleBin = GetSampleFromGlobalBin(SampleBinning, GlobalBin);
  int LocalBin  = GetLocalBinFromGlobalBin(SampleBinning, GlobalBin);
  return GetBinName(SampleBin, LocalBin);
}
