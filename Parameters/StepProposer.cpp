#include "Parameters/StepProposer.h"

StepProposer::StepProposer() {
  special_proposal.enabled = false;
  adaptive.enabled = false;

  rng.gaus = std::normal_distribution<double>(0);
  rng.unif = std::uniform_real_distribution<double>(0, 1);

  proposal_basis._propid = 0;
  systematic_basis._propid = std::numeric_limits<size_t>::max();
};

StepProposer::StepProposer(Eigen::MatrixXd proposal_matrix,
                           Eigen::ArrayXd current_values)
    : StepProposer() {
  SetProposalMatrix(proposal_matrix);
  SetSystematicParameterValues(current_values);
}

int StepProposer::GetProposalParameterIndexFromSystematicIndex(int i) const {
  if (i < systematic_basis.pca.first_index) {
    return i;
  } else if ((i >= systematic_basis.pca.first_index) &&
             (i <= systematic_basis.pca.last_index)) {
    return ParameterInPCABlock;
  } else {
    return i - systematic_basis.pca.nrotated_systematic_parameters() +
           systematic_basis.pca.npc_parameters();
  }
}

void StepProposer::SetProposalMatrix(Eigen::MatrixXd proposal_matrix) {
  proposal_basis.proposal = proposal_matrix;
  Eigen::LLT<Eigen::MatrixXd> lltproposal(proposal_basis.proposal);
  if (lltproposal.info() != Eigen::Success) {
    MACH3LOG_ERROR("Failed to LLT decompose proposal matrix");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  proposal_basis.l_proposal = lltproposal.matrixL();
}

void StepProposer::SetSystematicParameterValues(
    Eigen::ArrayXd const &new_values) {
  systematic_basis.proposed = new_values;
  TranformSystematicParametersProposalToBasis(systematic_basis.proposed,
                                              proposal_basis.proposed);
  systematic_basis._propid = 0;
  proposal_basis._propid = 0;
}
void StepProposer::GetSystematicParameterValues(
    Eigen::ArrayXd &proposed_values) {
  if (systematic_basis._propid != proposal_basis._propid) {
    TranformProposalParametersToSystematicBasis(proposal_basis.proposed,
                                                systematic_basis.proposed);
    systematic_basis._propid = proposal_basis._propid;
  }
  proposed_values = systematic_basis.proposed;
}

void StepProposer::TranformProposalParametersToSystematicBasis(
    Eigen::ArrayXd const &proposal_basis_vals,
    Eigen::ArrayXd &systematic_basis_vals) const {
  systematic_basis_vals.head(systematic_basis.pca.first_index) =
      proposal_basis_vals.head(systematic_basis.pca.first_index);

  // 1) rotate the principle component parameter values into 'deltas' in the
  //  systematic basis into the by using the PCA matrix.
  // 2) add on the offset (prefit) values from to these deltas to return to the
  // systematic basis.
  auto syst_param_delta = (systematic_basis.pca.pc_to_syst_rotation *
                           proposal_basis_vals
                               .segment(systematic_basis.pca.first_index,
                                        systematic_basis.pca.npc_parameters())
                               .matrix())
                              .array();

  systematic_basis_vals.segment(
      systematic_basis.pca.first_index,
      systematic_basis.pca.nrotated_systematic_parameters()) =
      systematic_basis.pca.offset.segment(
          systematic_basis.pca.first_index,
          systematic_basis.pca.nrotated_systematic_parameters()) +
      syst_param_delta;

  systematic_basis_vals.tail(systematic_basis.pca.ntail) =
      proposal_basis_vals.tail(systematic_basis.pca.ntail);
}

void StepProposer::TranformSystematicParametersProposalToBasis(
    Eigen::ArrayXd const &systematic_basis_vals,
    Eigen::ArrayXd &proposal_basis_vals) const {
  proposal_basis_vals.head(systematic_basis.pca.first_index) =
      systematic_basis_vals.head(systematic_basis.pca.first_index);

  // 1) subtract the offset (prefit) values from the current parameter values so
  // that the PC basis has a CV of 0
  // 2) rotate the resulting parameter 'deltas' in the systematic basis into the
  //  PC basis by using the transpose of the PCA matrix
  auto syst_param_delta =
      (systematic_basis_vals - systematic_basis.pca.offset)
          .segment(systematic_basis.pca.first_index,
                   systematic_basis.pca.nrotated_systematic_parameters())
          .matrix();

  proposal_basis_vals.segment(systematic_basis.pca.first_index,
                              systematic_basis.pca.npc_parameters()) =
      (systematic_basis.pca.syst_to_pc_rotation * syst_param_delta).array();

  proposal_basis_vals.tail(systematic_basis.pca.ntail) =
      systematic_basis_vals.tail(systematic_basis.pca.ntail);
}

Eigen::ArrayXd const &StepProposer::Propose() {

  Eigen::VectorXd random_vector =
      Eigen::VectorXd::NullaryExpr(proposal_basis.current.size(), [&](int i) {
        return proposal_basis.isfree[i] ? rng.gaus(rng.e1) : 0;
      });

  proposal_basis.proposed =
      proposal_basis.current +
      ((proposal_basis.l_proposal * random_vector).array() *
       proposal_basis.scale * proposal_basis.global_scale);

  proposal_basis._propid++;

  // if we don't have any special proposal functions, leave early
  if (!special_proposal.enabled) {
    return proposal_basis.proposed;
  }

  for (auto const &[i, low, up] : special_proposal.circ_bounds) {
    if (proposal_basis.isfree[i]) {
      continue;
    }
    if (proposal_basis.proposed[i] > up) {
      proposal_basis.proposed[i] =
          low + std::fmod(proposal_basis.proposed[i] - up, up - low);
    } else if (proposal_basis.proposed[i] < low) {
      proposal_basis.proposed[i] =
          up - std::fmod(low - proposal_basis.proposed[i], up - low);
    }
  }

  for (auto const &[i, pivot] : special_proposal.flips) {
    if (proposal_basis.isfree[i]) {
      continue;
    }
    if (rng.unif(rng.e1) < 0.5) {
      proposal_basis.proposed[i] = 2 * pivot - proposal_basis.proposed[i];
    }
  }

  return proposal_basis.proposed;
}
