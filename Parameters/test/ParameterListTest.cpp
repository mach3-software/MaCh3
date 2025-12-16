#include "Parameters/ParameterList.h"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cassert>

using namespace Catch::Matchers;

YAML::Node BuildConfig() {
  return YAML::Load(
      R"(
Systematics:
- Systematic: {
    Names: { FancyName: par1 },
    ParameterValues: { PreFitValue: 1 },
    Error: 0.5,
    StepScale: { MCMC: 0.1 },
    ParameterBounds: [ -3, 3 ],
    FlatPrior: False,
    SampleNames: [ "A", "B" ],
    FixParam: False,
    SpecialProposal: {
      CircularBounds: [ -2.5, 2.5 ],
      FlipParameter: 0
    },
    Correlations: [ { par2: 0.5 }, ]
  }
- Systematic: {
    Names: { FancyName: par2 },
    ParameterValues: { PreFitValue: 1 },
    Error: 0.5,
    StepScale: { MCMC: 0.1 },
    ParameterBounds: [ -2, 2 ],
    Correlations: [ { par1: 0.5 }, ]
  }
)");
}

TEST_CASE("MakeFromYAML", "[Construction]") {
  auto parlist = MakeFromYAML(BuildConfig());

  for (int i = 0; i < parlist.params.prefit.size(); ++i) {
    std::cout << parlist.SystematicParameterToString(i) << std::endl;
  }
}
