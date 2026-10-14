#include "Parameters/AdaptiveProposer.h"

void AdaptiveProposer::UpdateRunningEstimates(Eigen::ArrayXd const &params) {
  nsteps++;
  if (nsteps < covariance_estimation_start_step) {
    return;
  }

  nsteps_added++;
  if (nsteps_added == 1 || (nsteps_added == covariance_estimation_reset)) {
    nsteps_added = 1;
    last_parameter_mean = params.matrix();
    accepted_parameter_mean = last_parameter_mean;

    M2 = Eigen::MatrixXd::Zero(params.size(),params.size());
    accepted_parameter_covariance = M2;
    return;
  }

  // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
  accepted_parameter_mean =
  last_parameter_mean +
  (params.matrix() - last_parameter_mean) / double(nsteps_added);

  M2 += (params.matrix() - last_parameter_mean) *
    (params.matrix() - accepted_parameter_mean).transpose();

  last_parameter_mean = accepted_parameter_mean;

  accepted_parameter_covariance = M2 / double(nsteps_added - 1);
}

void AdaptiveProposer::Configure(YAML::Node const &adapt_manager, std::string const &matrix_name_str, int proposal_dimension) {

  matrix_name = matrix_name_str;
  auto const &matrix_config = adapt_manager["AdaptiveOptions"]["Covariance"][matrix_name_str];
  enabled = matrix_config.IsDefined() && GetFromManager<bool>(matrix_config["Enabled"], false);

  if (enabled) {
    covariance_estimation_start_step = GetFromManager<int>(adapt_manager["AdaptiveOptions"]["Settings"]["CovarianceEstimationStartStep"], 0);
    proposal_update_interval = GetFromManager<int>(adapt_manager["AdaptiveOptions"]["Settings"]["ProposalUpdateInterval"], 1000);
    covariance_estimation_reset = GetFromManager<int>(adapt_manager["AdaptiveOptions"]["Settings"]["CovarianceEstimationReset"], 0);
    proposal_matrix_output_file = GetFromManager<std::string>(adapt_manager["AdaptiveOptions"]["Settings"]["ProposalMatrixOutputFile"], "adaptive_proposal_matrix.root");
    proposal_scale = 2.38 / std::sqrt(double(proposal_dimension));
  }
}

void AdaptiveProposer::WriteProposalMatrixToFile(Eigen::MatrixXd const &proposal_matrix) {
  TFile output_file(proposal_matrix_output_file.c_str(), "UPDATE");
  if (output_file.IsZombie()) {
    MACH3LOG_ERROR("Couldn't open {}", proposal_matrix_output_file);
    return;
  }

  TMatrixDSym proposal_to_save(static_cast<int>(proposal_matrix.rows()));
  for (Eigen::Index row = 0; row < proposal_matrix.rows(); ++row) {
    for (Eigen::Index col = 0; col < proposal_matrix.cols(); ++col) {
      proposal_to_save(static_cast<int>(row), static_cast<int>(col)) = proposal_matrix(row, col);
    }
  }

  output_file.cd();
  std::string const matrix_key = std::to_string(nsteps) + "_" + matrix_name + "_throw_matrix";
  std::string const final_matrix_key = matrix_name + "_final_throw_matrix";
  if (proposal_to_save.Write(final_matrix_key.c_str(), TObject::kOverwrite) <= 0) {
    MACH3LOG_ERROR("Failed to write {} to {}", final_matrix_key, proposal_matrix_output_file);
  }
  
  if (proposal_to_save.Write(matrix_key.c_str()) <= 0) {
    MACH3LOG_ERROR("Failed to write {} to {}", matrix_key, proposal_matrix_output_file);
  }
  output_file.Close();
}