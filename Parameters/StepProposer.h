#pragma once

#include "Manager/Manager.h"

#include "TMatrixDSym.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Walloca"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#include "Eigen/Dense"
#pragma GCC diagnostic pop

#include <random>
#include <utility>
#include <vector>

struct StepProposer {
  struct {
    Eigen::ArrayXd current, proposed, scale;
    Eigen::ArrayXi isfree;
    double global_scale;
    Eigen::MatrixXd l_proposal;
  } params;

  struct {
    bool enabled;

    // param_index, flip pivot value
    std::vector<std::pair<int, double>> flips;
    // param_index, low bound, up bound
    std::vector<std::tuple<int, double, double>> circ_bounds;

  } special_proposal;

  struct {
    std::ranlux48 e1;
    std::normal_distribution<double> gaus;
    std::uniform_real_distribution<double> unif;
  } rng;

  struct {
    TMatrixDSym proposal;
  } root_copies;

  struct {
    bool enabled;
    size_t nsteps;
    Eigen::MatrixXd accepted_parameter_covariance;
  } adaptive;

  StepProposer();
  StepProposer(Eigen::MatrixXd proposal_matrix, Eigen::ArrayXd current_values);

  void SetProposalMatrix(Eigen::MatrixXd proposal_matrix);
  void SetParameterValues(Eigen::ArrayXd current_values) {
    params.current = current_values;
  }

  void EnableAdaption(YAML::Node const &adaption_config);
  void SetAdaptionCovariance(Eigen::MatrixXd parameter_covariance,
                             size_t nsteps) {
    adaptive.nsteps = nsteps;
    adaptive.accepted_parameter_covariance = parameter_covariance;
  }

  Eigen::ArrayXd const &Propose();
  void Accept();
};
