#include "Parameters/StepProposer.h"

StepProposer::StepProposer() {
  special_proposal.enabled = false;
  adaptive.enabled = false;

  rng.gaus = std::normal_distribution<double>(0);
  rng.unif = std::uniform_real_distribution<double>(0, 1);
};

StepProposer::StepProposer(Eigen::MatrixXd proposal_matrix,
                           Eigen::ArrayXd current_values)
    : StepProposer() {
  SetProposalMatrix(proposal_matrix);
  SetParameterValues(current_values);
}

Eigen::ArrayXd const &StepProposer::Propose() {

  Eigen::VectorXd random_vector =
      Eigen::VectorXd::NullaryExpr(params.current.size(), [&](int i) {
        return params.isfree[i] ? rng.gaus(rng.e1) : 0;
      });

  params.proposed =
      params.current + ((params.l_proposal * random_vector).array() *
                        params.scale * params.global_scale);

  // if we don't have any special proposal functions, leave early
  if (!special_proposal.enabled) {
    return params.proposed;
  }

  for (auto const &[i, low, up] : special_proposal.circ_bounds) {
    if (params.isfree[i]) {
      continue;
    }
    if (params.proposed[i] > up) {
      params.proposed[i] = low + std::fmod(params.proposed[i] - up, up - low);
    } else if (params.proposed[i] < low) {
      params.proposed[i] = up - std::fmod(low - params.proposed[i], up - low);
    }
  }

  for (auto const &[i, pivot] : special_proposal.flips) {
    if (params.isfree[i]) {
      continue;
    }
    if (rng.unif(rng.e1) < 0.5) {
      params.proposed[i] = 2 * pivot - params.proposed[i];
    }
  }

  return params.proposed;
}

void StepProposer::SetProposalMatrix(Eigen::MatrixXd proposal_matrix) {
  Eigen::LLT<Eigen::MatrixXd> lltproposal(proposal_matrix);
  if (lltproposal.info() != Eigen::Success) {
    MACH3LOG_ERROR("Failed to LLT decompose proposal matrix");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  params.l_proposal = lltproposal.matrixL();
}

void StepProposer::Accept() {}
