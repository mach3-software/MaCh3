#include "Parameters/ParameterList.h"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cassert>

using namespace Catch::Matchers;

TEST_CASE("MakeFromYAML", "[Construction]") {
  auto parlist = ParameterList::MakeFromYAML(YAML::Load(
      R"(
Systematics:
- Systematic: {
    Names: { FancyName: par1 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 0.1 },
    ParameterBounds: [ -3, 3 ],
    FlatPrior: False,
    SampleNames: [ "A", "B" ],
    FixParam: False,
    SpecialProposal: {
      CircularBounds: [ -2.5, 2.5 ],
      FlipParameter: 0
    },
    Correlations: [ { par2: 0.49 }, ]
  }
- Systematic: {
    Names: { FancyName: par2 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 0.1 },
    ParameterBounds: [ -2, 2 ],
    FlatPrior: True,
    Correlations: [ { par1: 0.49 }, ]
  }
- Systematic: {
    Names: { FancyName: par3 },
    ParameterValues: { PreFitValue: 0.25 },
    Error: 0.5,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -2, 2 ],
  }
)"));

  REQUIRE(parlist.FindParameter("par1") == 0);
  REQUIRE(parlist.FindParameter("par2") == 1);
  REQUIRE(parlist.FindParameter("par3") == 2);

  REQUIRE_THROWS_AS(parlist.FindParameter("par4"), MaCh3Exception);

  REQUIRE(parlist.params.name[0] == "par1");
  REQUIRE(parlist.params.lowbound[0] == -3);
  REQUIRE(parlist.params.upbound[0] == 3);
  REQUIRE(parlist.params.flip_pivot[0].first == true);
  REQUIRE(parlist.params.flip_pivot[0].second == 0);
  REQUIRE(std::get<0>(parlist.params.circ_bounds[0]) == true);
  REQUIRE(std::get<1>(parlist.params.circ_bounds[0]) == -2.5);
  REQUIRE(std::get<2>(parlist.params.circ_bounds[0]) == 2.5);

  REQUIRE(parlist.params.name[1] == "par2");
  REQUIRE(parlist.params.error[1] == 1);
  REQUIRE(parlist.params.flatprior[1] == 1);

  REQUIRE(parlist.params.name[2] == "par3");
  REQUIRE_THAT(parlist.params.prefit[2], WithinAbs(0.25, 1E-8));
  REQUIRE_THAT(parlist.params.error[2], WithinAbs(0.5, 1E-8));
  REQUIRE_THAT(parlist.params.stepscale[2], WithinAbs(1, 1E-8));

  REQUIRE_THAT(parlist.params.covariance(0, 0), WithinAbs(1, 1E-8));
  REQUIRE_THAT(parlist.params.covariance(1, 1), WithinAbs(1, 1E-8));
  REQUIRE_THAT(parlist.params.covariance(2, 2), WithinAbs(0.25, 1E-8));

  REQUIRE_THAT(parlist.params.covariance(1, 0), WithinAbs(0.49, 1E-8));
  REQUIRE_THAT(parlist.params.covariance(0, 1), WithinAbs(0.49, 1E-8));
  REQUIRE_THAT(parlist.params.covariance(0, 2), WithinAbs(0, 1E-8));
  REQUIRE_THAT(parlist.params.covariance(1, 2), WithinAbs(0, 1E-8));
}

TEST_CASE("MakeFromYAML Bad Correlations", "[Construction]") {
  REQUIRE_THROWS_AS(ParameterList::MakeFromYAML(YAML::Load(
                        R"(
Systematics:
- Systematic: {
    Names: { FancyName: par1 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 0.1 },
    ParameterBounds: [ -3, 3 ],
    Correlations: [ { par2: 0.49 }, ]
  }
- Systematic: {
    Names: { FancyName: par2 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 0.1 },
    ParameterBounds: [ -3, 3 ],
    Correlations: [ { par1: 0.51 }, ]
  }
)")),
                    MaCh3Exception);

  REQUIRE_THROWS_AS(ParameterList::MakeFromYAML(YAML::Load(
                        R"(
Systematics:
- Systematic: {
    Names: { FancyName: par1 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 0.1 },
    ParameterBounds: [ -3, 3 ],
    Correlations: [ { par2: 0.49 }, ]
  }
- Systematic: {
    Names: { FancyName: par2 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 0.1 },
    ParameterBounds: [ -3, 3 ],
    Correlations: [ { par3: 0.51 }, ]
  }
)")),
                    MaCh3Exception);

  REQUIRE_THROWS_AS(ParameterList::MakeFromYAML(YAML::Load(
                        R"(
Systematics:
- Systematic: {
    Names: { FancyName: par1 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
  }
- Systematic: {
    Names: { FancyName: par2 },
    ParameterValues: { PreFitValue: 2 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
    Correlations: [ { par3: 0.9999 }, ]
  }
- Systematic: {
    Names: { FancyName: par3 },
    ParameterValues: { PreFitValue: 3 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
    Correlations: [ { par2: 0.9999 },
                    { par1: 0.1 }, ]
  }
- Systematic: {
    Names: { FancyName: par4 },
    ParameterValues: { PreFitValue: 4 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
  }
)")),
                    MaCh3Exception);
}

TEST_CASE("4param", "[PCA]") {
  auto parlist = ParameterList::MakeFromYAML(YAML::Load(
      R"(
Systematics:
- Systematic: {
    Names: { FancyName: par1 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
  }
- Systematic: {
    Names: { FancyName: par2 },
    ParameterValues: { PreFitValue: 2 },
    Error: 2,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
    Correlations: [ { par3: 0.99999 }, ]
  }
- Systematic: {
    Names: { FancyName: par3 },
    ParameterValues: { PreFitValue: 3 },
    Error: 3,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
    Correlations: [ { par2: 0.99999 }, ]
  }
- Systematic: {
    Names: { FancyName: par4 },
    ParameterValues: { PreFitValue: 4 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
  }
)"));

  parlist.ConstructTruncatedPCA(1E-4, 1, 2);

  REQUIRE(parlist.pca.enabled);
  REQUIRE(parlist.pca.first_index == 1);
  REQUIRE(parlist.pca.last_index == 2);
  REQUIRE(parlist.pca.ntail == 1);
  REQUIRE(parlist.pca.nrotated_syst_parameters() == 2);
  REQUIRE(parlist.pca.npc_parameters() == 1);
  REQUIRE(parlist.NumPCBasisParameters() == 3);

  REQUIRE(parlist.SystematicParameterIndexToPCIndex(0) == 0);
  REQUIRE(parlist.SystematicParameterIndexToPCIndex(1) ==
          ParameterList::ParameterInPCABlock);
  REQUIRE(parlist.SystematicParameterIndexToPCIndex(2) ==
          ParameterList::ParameterInPCABlock);
  REQUIRE(parlist.SystematicParameterIndexToPCIndex(3) == 2);

  Eigen::ArrayXd syst_vals = parlist.params.prefit;
  Eigen::ArrayXd pc_vals = Eigen::ArrayXd::Zero(parlist.NumPCBasisParameters());

  parlist.RotateSystematicParameterValuesToPCBasis(syst_vals, pc_vals);

  // prefit value should get mapped to 0
  REQUIRE(pc_vals[1] == 0);

  parlist.RotatePCParameterValuesToSystematicBasis(pc_vals, syst_vals);

  // prefit value should get mapped and back to prefit
  REQUIRE(syst_vals[1] == 2);
  REQUIRE(syst_vals[2] == 3);

  // the PC is approximately in the par2+par3 direction, so
  //   should be able to guess where pc_vals[1] = 1 sets
  //   syst_vals[1] and syst_vals[2]
  pc_vals[1] = 1;

  parlist.RotatePCParameterValuesToSystematicBasis(pc_vals, syst_vals);

  REQUIRE_THAT(syst_vals[1], WithinAbs(4, 1E-4));
  REQUIRE_THAT(syst_vals[2], WithinAbs(6, 1E-4));

  // check that we can rotate back to the expected value
  parlist.RotateSystematicParameterValuesToPCBasis(syst_vals, pc_vals);
  REQUIRE_THAT(pc_vals[1], WithinAbs(1, 1E-4));
}

TEST_CASE("Out of block correlations", "[PCA]") {
  auto parlist = ParameterList::MakeFromYAML(YAML::Load(
      R"(
Systematics:
- Systematic: {
    Names: { FancyName: par1 },
    ParameterValues: { PreFitValue: 1 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
    Correlations: [ { par3: 0.3 }, ]
  }
- Systematic: {
    Names: { FancyName: par2 },
    ParameterValues: { PreFitValue: 2 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
    Correlations: [ { par3: 0.5 }, ]
  }
- Systematic: {
    Names: { FancyName: par3 },
    ParameterValues: { PreFitValue: 3 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
    Correlations: [ { par2: 0.5 },
                    { par1: 0.3 }, ]
  }
- Systematic: {
    Names: { FancyName: par4 },
    ParameterValues: { PreFitValue: 4 },
    Error: 1,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
  }
)"));
  REQUIRE_THROWS_AS(parlist.ConstructTruncatedPCA(1E-4, 1, 2), MaCh3Exception);
}
