#include "Parameters/ParameterList.h"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cassert>

using namespace Catch::Matchers;

TEST_CASE("proposal reconstruction", "[Proposal]") {
  auto parlist = ParameterList::MakeFromYAML(YAML::Load(
      R"(
Systematics:
- Systematic: {
    Names: { FancyName: par1 },
    ParameterValues: { PreFitValue: 1 },
    Error: 2.5,
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
    Error: 9,
    StepScale: { MCMC: 1 },
    ParameterBounds: [ -10, 10 ],
  }
)"));

  parlist.ConstructTruncatedPCA(1E-4, 1, 2);

  auto proposer = parlist.MakeProposer();

  REQUIRE(proposer.NumParameters() == parlist.NumPCBasisParameters());
  REQUIRE_FALSE(proposer.special_proposal.enabled);

  Eigen::MatrixXd proposal_reconstructed =
      Eigen::MatrixXd::Zero(proposer.NumParameters(), proposer.NumParameters());

  Eigen::ArrayXd syst_param_vals =
      Eigen::ArrayXd::Zero(parlist.NumSystematicBasisParameters());
  Eigen::MatrixXd syst_basis_covariance_reconstructed =
      Eigen::MatrixXd::Zero(parlist.NumSystematicBasisParameters(),
                            parlist.NumSystematicBasisParameters());

  Eigen::ArrayXd mean_proposed = Eigen::ArrayXd::Zero(proposer.NumParameters());
  Eigen::ArrayXd mean_syst_param_vals =
      Eigen::ArrayXd::Zero(parlist.NumSystematicBasisParameters());

  double nThrows = 1E7;
  for (int i = 0; i < nThrows; ++i) {

    // if we never accept a step then they will always we proposed around the
    // prefit
    Eigen::ArrayXd proposed_vals = proposer.Propose();
    mean_proposed += proposed_vals;

    Eigen::VectorXd pc_basis_delta =
        (proposed_vals - proposer.params.current).matrix();

    proposal_reconstructed += pc_basis_delta * pc_basis_delta.transpose();

    parlist.RotatePCParameterValuesToSystematicBasis(proposed_vals,
                                                     syst_param_vals);
    mean_syst_param_vals += syst_param_vals;
    Eigen::VectorXd syst_basis_deltas =
        (syst_param_vals - parlist.params.prefit).matrix();
    syst_basis_covariance_reconstructed +=
        syst_basis_deltas * syst_basis_deltas.transpose();
  }

  mean_proposed.array() /= nThrows;
  mean_syst_param_vals.array() /= nThrows;
  proposal_reconstructed.array() /= nThrows;
  syst_basis_covariance_reconstructed.array() /= nThrows;

  REQUIRE_THAT(mean_proposed[0], WithinAbs(1, 1E-2));
  REQUIRE_THAT(mean_proposed[1], WithinAbs(0, 1E-2));
  REQUIRE_THAT(mean_proposed[2], WithinAbs(4, 1E-2));

  REQUIRE_THAT(mean_syst_param_vals[0], WithinAbs(1, 1E-2));
  REQUIRE_THAT(mean_syst_param_vals[1], WithinAbs(2, 1E-2));
  REQUIRE_THAT(mean_syst_param_vals[2], WithinAbs(3, 1E-2));
  REQUIRE_THAT(mean_syst_param_vals[3], WithinAbs(4, 1E-2));

  Eigen::MatrixXd proposal_diff =
      proposal_reconstructed - proposer.params.proposal;

  for (int i = 0; i < proposal_diff.rows(); ++i) {
    for (int j = 0; j < proposal_diff.rows(); ++j) {
      REQUIRE_THAT(proposal_diff(i, j), WithinAbs(0, 1E-2));
    }
  }

  Eigen::MatrixXd syst_basis_covariance_diff =
      syst_basis_covariance_reconstructed - parlist.params.covariance;

  for (int i = 0; i < syst_basis_covariance_diff.rows(); ++i) {
    for (int j = 0; j < syst_basis_covariance_diff.rows(); ++j) {
      REQUIRE_THAT(syst_basis_covariance_diff(i, j), WithinAbs(0, 1E-2));
    }
  }
}
