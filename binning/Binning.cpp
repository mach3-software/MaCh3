#include "binning/Binning.h"

#include "manager/MaCh3Exception.h"
#include "manager/MaCh3Logger.h"

#include "spdlog/fmt/bundled/ranges.h"

std::vector<double> get_bin_edges_from_node(YAML::Node const &ax_node) {
  // builds a uniform binning from an nbins,start,stop type
  // specification
  if (ax_node["linspace"]) {

    int nbins = ax_node["linspace"]["n"].as<int>();
    double start = ax_node["linspace"]["low"].as<double>();
    double stop = ax_node["linspace"]["high"].as<double>();

    double width = (stop - start) / double(nbins);
    std::vector<double> bin_edges;
    bin_edges.push_back(start);
    for (int i = 0; i < nbins; ++i) {
      bin_edges.push_back(bin_edges.back() + width);
    }
    return bin_edges;
  } else if (ax_node["variable"]) {
    return ax_node["variable"].as<std::vector<double>>();
  } else {
    MACH3LOG_ERROR(
        "No valid axis binning definition found, valid values: linspace: "
        "{n: 10, low: 0, high: 10}, variable: [edge0, edge1, ....]");
    std::stringstream ss;
    ss << ax_node;
    MACH3LOG_ERROR("YAML::Dump follows:\n{}\n", ss.str());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

std::unique_ptr<Binning> Binning::MakeBinning(YAML::Node const &config) {
  if (config["axes"]) {
    return std::make_unique<Binning_Product1D>(config);
  } else if (config["bins"]) {
    return std::make_unique<Binning_ArbitraryHyperRect>(config);
  }
  MACH3LOG_ERROR("No valid binning definition found, expected a top level node "
                 "keyed \"axes\" or \"bins\".");
  std::stringstream ss;
  ss << config;
  MACH3LOG_ERROR("YAML::Dump follows:\n{}\n", ss.str());
  throw MaCh3Exception(__FILE__, __LINE__);
}

Binning_Product1D::Binning_Product1D(YAML::Node const &config) {
  nbins_per_slice = {
      1,
  };
  for (auto const &ax : config["axes"]) {
    auto edges = get_bin_edges_from_node(ax);
    Axes.emplace_back(edges.size() - 1, edges.data());
    nbins_per_slice.push_back(nbins_per_slice.back() * Axes.back().GetNbins());
  }
  if (!Axes.size()) {
    MACH3LOG_ERROR("After constructing Binning_Product1D, had no axes.");
    std::stringstream ss;
    ss << config;
    MACH3LOG_ERROR("YAML::Dump follows:\n{}\n", ss.str());
  }
}

int Binning_Product1D::GetBinNumber(
    std::vector<int> const &axis_bin_numbers) const {

  if (axis_bin_numbers.size() != Axes.size()) {
    MACH3LOG_ERROR("Binning_Product1D::GetBinNumber called with "
                   "axis_bin_numbers: {}, but the binning has {} axes.",
                   axis_bin_numbers, Axes.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  int gbin = 0;
  for (size_t ax_i = 0; ax_i < Axes.size(); ++ax_i) {

    int ax_bin = axis_bin_numbers[ax_i];

    // flow bin on this axis. Since we don't track per-axis flow bins, just
    // check if the bin number is valid and if not return the  flow bin at
    // gbin = -1
    if ((ax_bin < 0) || (ax_bin >= Axes[ax_i].GetNbins())) {
      return npos;
    }

    gbin += ax_bin * nbins_per_slice[ax_i];
  }
  return gbin;
}

int Binning_Product1D::GetBinNumber(std::vector<double> const &values,
                                    int gbi_hint) const {

  (void)gbi_hint;

  if (values.size() != Axes.size()) {
    MACH3LOG_ERROR("Binning_Product1D::GetBinNumber called with values: {}, "
                   "but the binning has {} axes.",
                   values, Axes.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::vector<int> bin_numbers(values.size());

  for (size_t ax_i = 0; ax_i < Axes.size(); ++ax_i) {

    int ax_bin = Axes[ax_i].FindFixBin(values[ax_i]);
    // flow bin on this axis. The flow bins here are in ROOT convention
    if ((ax_bin == 0) || (ax_bin == (Axes[ax_i].GetNbins() + 1))) {
      return npos;
    }
    // get rid of the ROOT underflow along each axis as we have a
    // underflow at gbin == -1
    bin_numbers[ax_i] = ax_bin - 1;
  }
  return GetBinNumber(bin_numbers);
}

bool Binning_Product1D::IsInBin(std::vector<double> const &values,
                                int gbi) const {
  if ((gbi < 0) || (gbi == npos)) {
    return false;
  }
  auto ax_bins = DecomposeBinNumber(gbi);

  if (!ax_bins.size()) {
    return false;
  }

  for (int axi = 0; axi < GetNDimensions(); ++axi) {
    if ((values[axi] < Axes[axi].GetBinLowEdge(ax_bins[axi] + 1)) ||
        values[axi] >= Axes[axi].GetBinUpEdge(ax_bins[axi] +
                                              1)) { // not in bin in this dim
      return false;
    }
  }
  return true;
}

std::vector<int> Binning_Product1D::DecomposeBinNumber(int gbi) const {
  if ((gbi < 0) || (gbi == npos)) {
    return {};
  }

  std::vector<int> axis_binning;
  for (int ax_i = int(Axes.size() - 1); ax_i >= 0; --ax_i) {
    int ax_bin = gbi / nbins_per_slice[ax_i];
    axis_binning.insert(axis_binning.begin(), ax_bin);
    gbi = gbi % nbins_per_slice[ax_i];
  }
  return axis_binning;
}

double Binning_Product1D::GetBinHyperVolume(int gbi) const {

  if ((gbi < 0) || (gbi == npos)) {
    return 0xdeadbeef;
  }

  double hv = 1;
  auto ax_bins = DecomposeBinNumber(gbi);

  if (!ax_bins.size()) {
    return 0xdeadbeef;
  }

  for (int axi = 0; axi < GetNDimensions(); ++axi) {
    if (ax_bins[axi] >= Axes[axi].GetNbins()) {
      return 0xdeadbeef;
    }

    hv *= Axes[axi].GetBinWidth(ax_bins[axi] +
                                1); // back to ROOT underflow convention
  }
  return hv;
}

std::string Binning_Product1D::to_string() const {
  std::stringstream ss;

  int ndims = GetNDimensions();
  ss << "axes: [\n";
  for (int iax = 0; iax < ndims; ++iax) {
    ss << "  [";
    ss << fmt::format("{:4.3g}", Axes[iax].GetBinLowEdge(1)) << ", ";
    for (int ibin = 0; ibin < Axes[iax].GetNbins(); ++ibin) {
      ss << fmt::format("{:4.3g}", Axes[iax].GetBinUpEdge(ibin + 1))
         << (((ibin + 1) == Axes[iax].GetNbins()) ? " ]" : ", ");
    }
    ss << (((iax + 1) == ndims) ? "\n" : ",\n");
  }
  ss << "]";

  return ss.str();
}

YAML::Node Binning_Product1D::to_YAML() const {
  YAML::Node doc = YAML::Load("");
  int ndims = GetNDimensions();
  for (int iax = 0; iax < ndims; ++iax) {
    std::vector<double> bin_edges = {
        Axes[iax].GetBinLowEdge(1),
    };
    for (int ibin = 0; ibin < Axes[iax].GetNbins(); ++ibin) {
      bin_edges.push_back(Axes[iax].GetBinUpEdge(ibin + 1));
    }
    doc["axes"].push_back(bin_edges);
  }
  return doc;
}

Binning_ArbitraryHyperRect::Binning_ArbitraryHyperRect(
    YAML::Node const &config) {
  for (auto const &bin : config["bins"]) {
    BinExtents.emplace_back(bin.as<std::vector<std::array<double, 2>>>());
    if (BinExtents.size() > 1) {
      if (int(BinExtents.back().size()) != GetNDimensions()) {
        MACH3LOG_ERROR("When constructing Binning_ArbitraryHyperRect, bin {} "
                       "had {} extents, but the dimensionality of the binning "
                       "(determined from the first bin) was {}.",
                       (BinExtents.size() - 1), BinExtents.back().size(),
                       GetNDimensions());
        std::stringstream ss;
        ss << config;
        MACH3LOG_ERROR("YAML::Dump follows:\n{}\n", ss.str());
      }
    }
  }
  if (!BinExtents.size()) {
    MACH3LOG_ERROR(
        "After constructing Binning_ArbitraryHyperRect, had no bins.");
    std::stringstream ss;
    ss << config;
    MACH3LOG_ERROR("YAML::Dump follows:\n{}\n", ss.str());
  }
}

int Binning_ArbitraryHyperRect::GetBinNumber(std::vector<double> const &values,
                                             int gbi_hint) const {
  (void)gbi_hint;

  int ndims = GetNDimensions();
  for (int ibin = 0; ibin < int(BinExtents.size()); ++ibin) {
    for (int iax = 0; iax < ndims; ++iax) {
      if ((values[iax] < BinExtents[ibin][iax][0]) ||
          values[iax] >= BinExtents[ibin][iax][1]) { // not in bin in this dim
        break;
      }
      if ((iax + 1) == ndims) { // last dimension in bin, return
        return ibin;
      }
    }
  }

  return Binning::npos;
}

bool Binning_ArbitraryHyperRect::IsInBin(std::vector<double> const &values,
                                         int gbi) const {

  if ((gbi < 0) || (gbi >= int(BinExtents.size()))) {
    return false;
  }

  int ndims = GetNDimensions();
  for (int iax = 0; iax < ndims; ++iax) {
    if ((values[iax] < BinExtents[gbi][iax][0]) ||
        values[iax] >= BinExtents[gbi][iax][1]) { // not in bin in this dim
      return false;
    }
  }
  return true;
}

std::string Binning_ArbitraryHyperRect::to_string() const {
  std::stringstream ss;

  int ndims = GetNDimensions();
  ss << "bins: [\n";
  for (int ibin = 0; ibin < int(BinExtents.size()); ++ibin) {
    ss << "  " << ibin << ": [";
    for (int iax = 0; iax < ndims; ++iax) {
      ss << fmt::format("[{:4.3g},{:4.3g})", BinExtents[ibin][iax][0],
                        BinExtents[ibin][iax][1])
         << (((iax + 1) == ndims) ? " ]" : ", ");
    }
    ss << (((ibin + 1) == int(BinExtents.size())) ? "\n" : ",\n");
  }
  ss << "]";
  return ss.str();
}

YAML::Node Binning_ArbitraryHyperRect::to_YAML() const {
  YAML::Node doc = YAML::Load("");
  int ndims = GetNDimensions();
  for (int ibin = 0; ibin < int(BinExtents.size()); ++ibin) {
    YAML::Node bin;
    for (int iax = 0; iax < ndims; ++iax) {
      bin.push_back(BinExtents[ibin][iax]);
    }
    doc["binning"].push_back(bin);
  }
  return doc;
}
