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
  SingleSample.VarStr = Get<std::vector<std::string>>(Settings["VarStr"], __FILE__ , __LINE__);
  SingleSample.nDimensions = static_cast<int>(SingleSample.VarStr.size());

  SampleBinningInfo SingleBinning;
  bool Uniform = Get<bool>(Settings["Uniform"], __FILE__ , __LINE__);
  if(Uniform == false) {
    auto Bins = Get<std::vector<std::vector<std::vector<double>>>>(Settings["Bins"], __FILE__, __LINE__);
    SingleBinning.InitNonUniform(Bins);
  } else {
    auto Bin_Edges = Get<std::vector<std::vector<double>>>(Settings["VarBins"], __FILE__ , __LINE__);
    SingleBinning.InitUniform(Bin_Edges);
  }
  if(SingleSample.VarStr.size() != SingleBinning.BinEdges.size()) {
    MACH3LOG_ERROR("Number of variables ({}) does not match number of bin edge sets ({}) in sample config '{}'",
                   SingleSample.VarStr.size(), SingleBinning.BinEdges.size(),SingleSample.SampleTitle);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

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
  const SampleBinningInfo& _restrict_ Binning = SampleBinning[NomSample];
  int GlobalBin = 0;

  for(int i = 0; i < Dim; ++i) {
    const double Var = *KinVar[i];
    const int Bin = Binning.FindBin(i, Var, NomBin[i]);
    // KS: If we are outside of range in only one dimension this mean out of bounds, we can simply quickly finish
    if(Bin < 0) return M3::UnderOverFlowBin;
    // KS: inline GetBin computation to avoid any memory allocation, which in reweight loop is very costly
    GlobalBin += Bin * Binning.Strides[i];
  }

  if(Binning.Uniform) {
    GlobalBin += static_cast<int>(Binning.GlobalOffset);
    return GlobalBin;
  } else{
    const auto& _restrict_ BinMapping = Binning.BinGridMapping[GlobalBin];
    for(size_t iBin = 0; iBin < BinMapping.size(); iBin++) {
      const int BinNumber = BinMapping[iBin];
      const auto& _restrict_ NonUniBin = Binning.Bins[BinNumber];
      if(NonUniBin.IsEventInside(KinVar)){
        return BinNumber + Binning.GlobalOffset;
      }
    }
    MACH3LOG_DEBUG("Didn't find any bin so returning UnderOverFlowBin");
    return M3::UnderOverFlowBin;
  }
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
int BinningHandler::GetBinSafe(const int Sample, const std::vector<int>& Bins) const {
// ************************************************
  const int GlobalBin = SampleBinning[Sample].GetBinSafe(Bins);
  return GlobalBin;
}

// ************************************************
int BinningHandler::GetGlobalBinSafe(const int Sample, const std::vector<int>& Bins) const {
// ************************************************
  const int GlobalBin = SampleBinning[Sample].GetBinSafe(Bins) + static_cast<int>(SampleBinning[Sample].GlobalOffset);
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

  int GlobalOffsetCounter = 0;
  for(size_t iSample = 0; iSample < SampleBinning.size(); iSample++){
    SampleBinning[iSample].GlobalOffset = GlobalOffsetCounter;
    GlobalOffsetCounter += SampleBinning[iSample].nBins;
  }
  // lastly modify total number of bins
  TotalNumberOfBins = GlobalOffsetCounter;
}

// ************************************************
// Get fancy name for a given bin, to help match it with global properties
std::string BinningHandler::GetBinName(const int iSample, const std::vector<int>& Bins) const {
// ************************************************
  const auto& Binning = SampleBinning[iSample];
  if(!Binning.Uniform) {
    MACH3LOG_ERROR("When using Non-Uniform binning for sample {} please use One bin instead of Axis bins", iSample);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return GetBinName(iSample, GetBinSafe(iSample, Bins));
}

// ************************************************
// Get fancy name for a given bin, to help match it with global properties
std::string BinningHandler::GetBinName(const int iSample, const int SampleBin) const {
// ************************************************
  const auto& Binning = SampleBinning[iSample];

  // Safety checks
  if (SampleBin < 0 || SampleBin >= static_cast<int>(Binning.nBins)) {
    MACH3LOG_ERROR("Requested bin {} is out of range for sample {}", SampleBin, iSample);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  std::string BinName;

  if(Binning.Uniform) {
    int Dim = static_cast<int>(Binning.Strides.size());
    std::vector<int> Bins(Dim, 0);
    int Remaining = SampleBin;

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

    for (int i = 0; i < Dim; ++i) {
      if (i > 0) BinName += ", ";
      const double min = Binning.BinEdges[i].at(Bins[i]);
      const double max = Binning.BinEdges[i].at(Bins[i] + 1);
      BinName += fmt::format("Dim{} ({:g}, {:g})", i, min, max);
    }
  } else{
    const BinInfo& bin = Binning.Bins[SampleBin];
    const int Dim = static_cast<int>(bin.Extent.size());

    for (int i = 0; i < Dim; ++i) {
      if (i > 0) BinName += ", ";
      const double min = bin.Extent[i][0];
      const double max = bin.Extent[i][1];
      BinName += fmt::format("Dim{} ({:g}, {:g})", i, min, max);
    }
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

// ************************************************
int BinningHandler::GetNAxisBins(const int iSample, const int iDim) const {
// ************************************************
  const auto& Binning = SampleBinning[iSample];
  if(!Binning.Uniform) {
    MACH3LOG_ERROR("When using Non-Uniform binning for sample {} please use global bin instead of {}", iSample, __func__);
    throw MaCh3Exception(__FILE__, __LINE__);
  } else{
    return static_cast<int>(Binning.AxisNBins.at(iDim));
  }
}

// ************************************************
bool BinningHandler::IsUniform(const int iSample) const {
// ************************************************
  const auto& Binning = SampleBinning[iSample];
  return Binning.Uniform;
}

// ************************************************
const std::vector<BinInfo> BinningHandler::GetNonUniformBins(const int iSample) const {
// ************************************************
  const auto& Binning = SampleBinning[iSample];
  if(!Binning.Uniform) {
    return Binning.Bins;
  } else{
    MACH3LOG_ERROR("{} for sample {} will not work becasue binnin is unfiorm", __func__, iSample);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

