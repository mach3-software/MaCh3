#pragma once

#include "TAxis.h"

#include "yaml-cpp/yaml.h"

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

struct Binning {
  constexpr static int const npos = std::numeric_limits<int>::max();

  // Some algorithms can perform a pre-search if given the approximate bin
  virtual int GetBinNumber(std::vector<double> const &values,
                           int gbi_hint = npos) const = 0;
  virtual bool IsInBin(std::vector<double> const &values, int gbi) const = 0;

  virtual int GetNDimensions() const = 0;
  virtual double GetBinHyperVolume(int gbi) const = 0;

  virtual std::string to_string() const = 0;
  virtual YAML::Node to_YAML() const = 0;

  virtual ~Binning() {}
  // factory function for making binnings
  static std::unique_ptr<Binning> MakeBinning(YAML::Node const &);
};

struct Binning_Product1D : public Binning {
private:
  // for each axis tells you how many bins to step to get to the next bin
  std::vector<int> nbins_per_slice;
  std::vector<TAxis> Axes;

  std::vector<int> DecomposeBinNumber(int gbi) const;
  int GetBinNumber(std::vector<int> const &axis_bin_numbers) const;

public:
  Binning_Product1D(YAML::Node const &);
  int GetBinNumber(std::vector<double> const &values,
                   int gbi_hint = npos) const;
  bool IsInBin(std::vector<double> const &values, int gbi) const;

  int GetNDimensions() const { return int(Axes.size()); }
  double GetBinHyperVolume(int gbi) const;

  std::string to_string() const;
  YAML::Node to_YAML() const;
};

struct Binning_ArbitraryHyperRect : public Binning {

  // BinExtents[bin][dimension] = {low, high};
  std::vector<std::vector<std::array<double, 2>>> BinExtents;

public:
  Binning_ArbitraryHyperRect(YAML::Node const &);
  int GetBinNumber(std::vector<double> const &values,
                   int gbi_hint = npos) const;
  bool IsInBin(std::vector<double> const &values, int gbi) const;

  int GetNDimensions() const {
    return int(BinExtents.size() ? BinExtents.front().size() : 0);
  }

  double GetBinHyperVolume(int gbi) const {
    if ((gbi < 0) || (gbi >= int(BinExtents.size()))) {
      return 0xdeadbeef;
    }
    return std::accumulate(
        BinExtents[gbi].begin(), BinExtents[gbi].end(), double(1),
        [](double accum, std::array<double, 2> const &extent) -> double {
          return accum * (extent[1] - extent[0]);
        });
  }

  std::string to_string() const;
  YAML::Node to_YAML() const;
};
