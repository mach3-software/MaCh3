#include "binning/Binning.h"

#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/ranges.h"

#include <iostream>

int main() {
  auto config_2Dproduct = YAML::Load(R"(
binning:
  axes:
    - variable: [0,1,2,3,4,5,6,7,8,9,10]
    - linspace: { n: 10, low: 0, high: 20 }
)");

  auto binning_2dproduct = Binning::MakeBinning(config_2Dproduct["binning"]);

  std::cout << binning_2dproduct->to_string() << std::endl;

  std::vector<std::vector<double>> test_values = {
      {
          0.5,
          0.5,
      },
      {
          1,
          0.5,
      },
      {
          0.5,
          1,
      },
      {-1, -1},
      {10, 5},
      {5, 10},
      {10, 10},
      {100, 100},
  };

  for (auto const &point : test_values) {
    auto gbin = binning_2dproduct->GetBinNumber(point);
    std::cout << fmt::format("{} -> gbin = {}, gbin_HV = {}", point, gbin,
                             binning_2dproduct->GetBinHyperVolume(gbin))
              << std::endl;
  }

  auto config_2Darb = YAML::Load(R"(
binning:
  bins:
    - [ [0,1], [0,2] ]
    - [ [1,1.5], [0,1] ]
    - [ [1,1.5], [1,2] ]
    - [ [1,1.5], [1,2] ]
    - [ [1.5,3], [0,1] ]
    - [ [1.5,3], [1,2] ]
    - [ [0,3], [2,3] ]
)");

  auto binning_2darb = Binning::MakeBinning(config_2Darb["binning"]);

  std::cout << binning_2darb->to_string() << std::endl;

  for (auto const &point : test_values) {
    auto gbin = binning_2darb->GetBinNumber(point);
    std::cout << fmt::format("{} -> gbin = {}, HV = {}", point, gbin,
                             binning_2darb->GetBinHyperVolume(gbin))
              << std::endl;
  }
}