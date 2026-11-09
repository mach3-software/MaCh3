#pragma once

// MaCh3 includes
#include "Manager/Manager.h"
#include "Manager/Random.h"

namespace M3 {
/// @brief Builds a single dimension's bin edges from YAML::Node
/// @param bin_edges_node Yaml node containing binning
/// @param found_range_specifier Whether or not a range specifier (linspace/logspace) was found
std::vector<double> BuildBinEdgesFromNode(YAML::Node const &bin_edges_node, bool &found_range_specifier);

/// @brief Converts a range (linspace/logspace) to a std::vector<double>
/// @param bin_range Yaml node containing binning
std::vector<double> BinRangeToBinEdges(YAML::Node const &bin_range);

/// @brief KS: Convert TH2Poly into yaml config accepted by MaCh3
YAML::Node PolyToYaml(TH2Poly* Hist, const std::string& YamlName, const std::string& File, const int Line);
} //end M3
