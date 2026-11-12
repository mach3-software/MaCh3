#include "Utils/BinningUtils.h"
namespace M3 {
// ************************************************
std::vector<double> BinRangeToBinEdges(YAML::Node const &bin_range) {
// ************************************************
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

// ************************************************
/// @brief Builds a single dimension's bin edges from YAML::Node
/// @details
/// BinEdges:  [ <dim0bin0lowedge>, <dim0bin1upedge>, <dim0bin2upedge>, ...
/// <dim0binNupedge> ] BinEdges:  { linspace: { nb: 100, low: 0, up: 10} }
/// BinEdges:  { logspace: { nb: 100, low: 1E-1, up: 10} }
/// BinEdges:  [ { linspace: { nb: 100, low: 0, up: 10} }, 10, 15, { logspace: {
/// nb: 5, low: 15, up: 100} } ]
std::vector<double> BuildBinEdgesFromNode(YAML::Node const &bin_edges_node,
                           bool &found_range_specifier) {
// ************************************************
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
        std::stringstream ss;
        ss << "[ ";
        for(auto const & e : edges_builder){
          ss << fmt::format("{:.3g} ", e);
        }
        ss << "]";
        MACH3LOG_ERROR(
            "When parsing binning, found edges that were not monotonically "
            "increasing, problem bin at index: {}:\n{}",
            eb_it, ss.str());
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    edges.push_back(edges_builder[eb_it]);
  }

  return edges;
}

// ***************************************************************************
//KS: Convert TH2Poly into yaml config accepted by MaCh3
YAML::Node PolyToYaml(TH2Poly* Hist, const std::string& YamlName, const std::string& File, const int Line) {
// ***************************************************************************
  if (!Hist) {
    MACH3LOG_ERROR("Null TH2Poly pointer");
    throw MaCh3Exception(File, Line);
  }

  YAML::Node bins(YAML::NodeType::Sequence);
  bins.SetStyle(YAML::EmitterStyle::Flow);
  const int NBins = Hist->GetNumberOfBins();

  for (int j = 1; j <= NBins; j++)
  {
    TH2PolyBin* polybin = static_cast<TH2PolyBin*>(Hist->GetBins()->At(j - 1));

    double xmin = polybin->GetXMin();
    double xmax = polybin->GetXMax();
    double ymin = polybin->GetYMin();
    double ymax = polybin->GetYMax();

    YAML::Node xNode(YAML::NodeType::Sequence);
    xNode.SetStyle(YAML::EmitterStyle::Flow);
    xNode.push_back(xmin);
    xNode.push_back(xmax);

    YAML::Node yNode(YAML::NodeType::Sequence);
    yNode.SetStyle(YAML::EmitterStyle::Flow);
    yNode.push_back(ymin);
    yNode.push_back(ymax);

    YAML::Node bin(YAML::NodeType::Sequence);
    bin.SetStyle(YAML::EmitterStyle::Flow);
    bin.push_back(xNode);
    bin.push_back(yNode);

    bins.push_back(bin);
  }

  YAML::Node result;
  result[YamlName] = bins;

  return result;
}
} //end M3 namespace
