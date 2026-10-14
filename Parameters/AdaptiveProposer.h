#pragma once

#include "Manager/Manager.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Walloca"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wduplicated-branches"
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include "Eigen/Dense"
#pragma GCC diagnostic pop


struct AdaptiveProposer {
    bool enabled = false;
    int nsteps = 0;
    int nsteps_added = 0;
    int proposal_update_interval = 1000;
    int covariance_estimation_start_step = 0;
    int covariance_estimation_reset = 0;
    double proposal_scale = 1.0;
    std::string matrix_name;
    std::string proposal_matrix_output_file;

    Eigen::VectorXd accepted_parameter_mean, last_parameter_mean;
    Eigen::MatrixXd accepted_parameter_covariance;
    Eigen::MatrixXd M2;

    bool ShouldUpdateProposalMatrix(int proposal_dimension) const {
      return proposal_update_interval > 0 &&
      nsteps_added > 0 &&
      nsteps_added % proposal_update_interval == 0 &&
      nsteps_added >= proposal_dimension + 1;
    }

    void UpdateRunningEstimates(Eigen::ArrayXd const &params);
    
    void Configure(YAML::Node const &adapt_manager, std::string const &matrix_name_str, int proposal_dimension);

    void WriteProposalMatrixToFile(Eigen::MatrixXd const &proposal_matrix);
  };