#include "Samples/BinningHandler.h"

#include "spdlog/fmt/ranges.h"

// ************************************************
BinningHandler::BinningHandler() {
// ************************************************
}

auto BinRangeToBinEdges(YAML::Node const &bin_range) {
  bool is_lin = true;
  YAML::Node bin_range_specifier;
  if (bin_range["linspace"]) {
    bin_range_specifier = bin_range["linspace"];
  } else if (bin_range["logspace"]) {
    is_lin = false;
    bin_range_specifier = bin_range["logspace"];
  } else {
    std::stringstream ss;
    ss << bin_range;
    MACH3LOG_ERROR("When parsing binning, expected bin range specifier with "
                   "key linspace or logspace, but found,\n{}",
                   ss.str());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  auto nb = Get<int>(bin_range_specifier["nb"], __FILE__, __LINE__);
  auto low = Get<double>(bin_range_specifier["low"], __FILE__, __LINE__);
  auto up = Get<double>(bin_range_specifier["up"], __FILE__, __LINE__);

  std::vector<double> edges(nb + 1, low);
  // force the last bin to be exactly as parsed to avoid numerical instabilities
  // not quite lining back up with the end of the range, which could
  // cause spurious errors or infinitesimally small bins for specifications
  // like: [ { logspace: { nb: 10, 1E-1, 10}, 10, 11} ]
  edges.back() = up;

  if (is_lin) {
    double bw = (up - low) / nb;
    for (int i = 0; i < (nb - 1); ++i) {
      edges[i + 1] = edges[i] + bw;
    }
  } else {
    double llow = std::log10(low);
    double lup = std::log10(up);
    double lbw = (lup - llow) / nb;
    for (int i = 0; i < (nb - 1); ++i) {
      edges[i + 1] = std::pow(10, llow + (i + 1) * lbw);
    }
  }

  return edges;
}

/// @brief Builds a single dimension's bin edges from YAML::Node
/// @details
/// BinEdges:  [ <dim0bin0lowedge>, <dim0bin1upedge>, <dim0bin2upedge>, ...
/// <dim0binNupedge> ] BinEdges:  { linspace: { nb: 100, low: 0, up: 10} }
/// BinEdges:  { logspace: { nb: 100, low: 1E-1, up: 10} }
/// BinEdges:  [ { linspace: { nb: 100, low: 0, up: 10} }, 10, 15, { logspace: {
/// nb: 5, low: 15, up: 100} } ]
auto BuildBinEdgesFromNode(YAML::Node const &bin_edges_node,
                           bool &found_range_specifier) {
  if (bin_edges_node.IsMap()) {
    found_range_specifier = true;
    return BinRangeToBinEdges(bin_edges_node);
  }
  std::vector<double> edges_builder;
  if (!bin_edges_node.IsSequence()) {
    std::stringstream ss;
    ss << bin_edges_node;
    MACH3LOG_ERROR(
        "When parsing binning, expected to find a YAML map or sequence, "
        "but found:\n{}",
        ss.str());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for (auto const &it : bin_edges_node) {
    if (it.IsScalar()) {
      edges_builder.push_back(it.as<double>());
    } else if (it.IsMap()) {
      found_range_specifier = true;
      auto range_edges = BinRangeToBinEdges(it);
      std::copy(range_edges.begin(), range_edges.end(),
                std::back_inserter(edges_builder));
    } else {
      std::stringstream ss;
      ss << bin_edges_node;
      MACH3LOG_ERROR(
          "When parsing binning, expected elements in outer sequence to all be "
          "either scalars or maps, but found:\n{}",
          ss.str());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  // Check for duplicates or out-of-order bins
  std::vector<double> edges;
  for (size_t eb_it = 0; eb_it < edges_builder.size(); ++eb_it) {
    if (edges.size()) {
      if (edges_builder[eb_it] == edges.back()) { // remove duplicate edges
        continue;
      } else if (edges_builder[eb_it] < edges.back()) {
        MACH3LOG_ERROR(
            "When parsing binning, found edges that were not monotonically "
            "increasing, problem bin at index: {}:\n{}",
            eb_it, edges_builder);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    edges.push_back(edges_builder[eb_it]);
  }

  return edges;
}

/// @brief Parses YAML node describing multidim uniform binning
/// @details
/// # dimensional list implicit for 1D binnings
/// BinEdges:  [<dim0bin0lowedge>, <dim0bin1upedge>, <dim0bin2upedge>, ...
/// <dim0binNupedge>]
///
/// BinEdges: [ [<dim0bin0lowedge>, <dim0bin1upedge>, <dim0bin2upedge>, ...
/// <dim0binNupedge>],
///            ...
///            [<dimNbin0lowedge>, <dimNbin1upedge>, <dimNbin2upedge>, ...
///            <dimNbinNupedge>] ]
///
/// # bin edge list implicit for 1D range-only binnings
/// BinEdges: { linspace: { nb: 100, low: 0, up: 10} }
/// # mixed syntax binnings allowed
/// BinEdges: [ { linspace: { nb: 100, low: 0, up: 10} }, 10, 15, { logspace: {
/// nb: 5, low: 15, up: 100} }  ] # for ND range-only binnings, lists are
/// required to disambiguate 1D mixed specifier binnings from ND binnings
/// BinEdges: [ [ { linspace: { nb: 100, low: 0, up: 10} } ],
///            ...
///            [<dimNbin0lowedge>, <dimNbin1upedge>, <dimNbin2upedge>, ...
///            <dimNbinNupedge>] ]
/// BinEdges: [ [ { linspace: { nb: 100, low: 0, up: 10} }, 10, 15, { logspace:
/// { nb: 5, low: 15, up: 100} }  ],
///            ...
///            [<dimNbin0lowedge>, <dimNbin1upedge>, <dimNbin2upedge>, ...
///            <dimNbinNupedge>] ]
auto UniformBinEdgeConfigParser(YAML::Node const &bin_edges_node,
                                bool &found_range_specifier) {
  if (bin_edges_node.IsMap()) {
    found_range_specifier = true;
    return std::vector<std::vector<double>>{
        BinRangeToBinEdges(bin_edges_node),
    };
  } else if (bin_edges_node.IsSequence()) {
    if (!bin_edges_node.size()) {
      std::stringstream ss;
      ss << bin_edges_node;
      MACH3LOG_ERROR("When parsing binning, found an empty sequence:\n{}",
                     ss.str());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    auto const &first_el = bin_edges_node[0];
    if (first_el.IsScalar() || first_el.IsMap()) { // 1D binning
      return std::vector<std::vector<double>>{
          BuildBinEdgesFromNode(bin_edges_node, found_range_specifier),
      };
    }
    // ND binning
    std::vector<std::vector<double>> dims;
    for (auto const &dim_node : bin_edges_node) {
      dims.push_back(BuildBinEdgesFromNode(dim_node, found_range_specifier));
    }
    return dims;
  } else {
    std::stringstream ss;
    ss << bin_edges_node;
    MACH3LOG_ERROR(
        "When parsing binning, expected to find a YAML map or sequence, "
        "but found:\n{}",
        ss.str());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
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
  bool found_range_specifier = false;
  if(Uniform == false) {
    auto Bins = Get<std::vector<std::vector<std::vector<double>>>>(Settings["Bins"], __FILE__, __LINE__);
    SingleBinning.InitNonUniform(Bins);
  } else {
    YAML::Node const & bin_edges_node = Settings["BinEdges"] ? Settings["BinEdges"] : Settings["VarBins"];
    if(!bin_edges_node){
      MACH3LOG_ERROR("When setting up Uniform sample binning, didn't find expected key: BinEdges (or VarBins for backward compatibility).");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    SingleBinning.InitUniform(UniformBinEdgeConfigParser(bin_edges_node, found_range_specifier));
  }
  if(SingleSample.VarStr.size() != SingleBinning.BinEdges.size()) {
    MACH3LOG_ERROR("Number of variables ({}) does not match number of bin edge sets ({}) in sample config '{}'",
                   SingleSample.VarStr.size(), SingleBinning.BinEdges.size(),SingleSample.SampleTitle);
    if(found_range_specifier){
      std::stringstream ss;
      ss << (Settings["BinEdges"] ? Settings["BinEdges"] : Settings["VarBins"]);
      MACH3LOG_ERROR(R"(A bin range specifier was found in node:

  {}

Please carefully check the number of square brackets used, a 2D binning
  comprised of just 2 range specifiers must explicitly include the axis
  list specifier like:

  BinEdges: [ [ {{linspace: {{nb:10, low:0, up: 10}}}} ], [ {{linspace: {{nb:5, low:10, up: 100}}}} ] ]
)", ss.str());
    }
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
  } else {
    const auto& _restrict_ BinMapping = Binning.BinGridMapping[GlobalBin];
    const size_t nNonUniBins = BinMapping.size();
    for(size_t iBin = 0; iBin < nNonUniBins; iBin++) {
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

