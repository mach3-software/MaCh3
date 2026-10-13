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

  auto const &matrix_config = adapt_manager["AdaptiveOptions"]["Covariance"][matrix_name_str];

  enabled = matrix_config.IsDefined() && GetFromManager<bool>(matrix_config["Enabled"], false);

  if (enabled) {
    covariance_estimation_start_step = GetFromManager<int>(adapt_manager["AdaptiveOptions"]["Settings"]["CovarianceEstimationStartStep"], 0);
    proposal_update_interval = GetFromManager<int>(adapt_manager["AdaptiveOptions"]["Settings"]["ProposalUpdateInterval"], 1000);
    covariance_estimation_reset = GetFromManager<int>(adapt_manager["AdaptiveOptions"]["Settings"]["CovarianceEstimationReset"], 0);
    proposal_scale = 2.38 / std::sqrt(double(proposal_dimension));
  }
}