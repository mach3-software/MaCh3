#pragma once

#include "Manager/Manager.h"

#include "TMatrixDSym.h"

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

#include <random>
#include <utility>
#include <vector>

struct StepProposer {
  struct {
    Eigen::ArrayXd current, proposed, scale;
    Eigen::ArrayXi isfree;
    Eigen::MatrixXd proposal;
    Eigen::MatrixXd l_proposal;
    double global_scale;

    // Used to keep track of if we have already transformed the current proposed
    // parameters to the systematic_basis
    size_t _propid;

  } proposal_basis;

  struct {
    Eigen::ArrayXd proposed;

    // Used to keep track of if we have already transformed the current proposed
    // parameters to the systematic_basis
    size_t _propid;

    Eigen::ArrayXi isfree;

    struct {
      int first_index, last_index, ntail;
      Eigen::ArrayXd offset;
      Eigen::MatrixXd pc_to_syst_rotation;
      Eigen::MatrixXd syst_to_pc_rotation;
      int nrotated_systematic_parameters() const {
        return int(pc_to_syst_rotation.rows());
      }
      int npc_parameters() const { return int(pc_to_syst_rotation.cols()); }
    } pca;

  } systematic_basis;

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
    bool enabled;
    int nsteps = 0;
    int nsteps_added = 0;
    int proposal_update_interval = 1000;
    int covariance_estimation_start_step = 0;
    double proposal_scale = 1.0;
    Eigen::VectorXd accepted_parameter_mean, last_parameter_mean;
    Eigen::MatrixXd accepted_parameter_covariance;
    Eigen::MatrixXd M2;
    bool ShouldUpdateProposalMatrix() const {
      return proposal_update_interval > 0 && nsteps_added > 0 && nsteps_added % proposal_update_interval == 0;
    }

    void UpdateRunningEstimates(Eigen::ArrayXd const &params) {
      nsteps++;
      if (nsteps < covariance_estimation_start_step) {
        return;
      }

      nsteps_added++;
      if (nsteps_added == 1) {
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
  } adaptive;

  StepProposer();
  StepProposer(Eigen::MatrixXd proposal_matrix, Eigen::ArrayXd current_values);

  int NumSystematicBasisParameters() const {
    return int(proposal_basis.current.size());
  }
  int NumProposalBasisParameters() const {
    return NumSystematicBasisParameters() -
           systematic_basis.pca.nrotated_systematic_parameters() +
           systematic_basis.pca.npc_parameters();
  }

  void SetProposalMatrix(Eigen::MatrixXd proposal_matrix);

  void SetSystematicParameterValues(Eigen::ArrayXd const &new_values);
  void GetSystematicParameterValues(Eigen::ArrayXd &proposed_values);

  void TranformProposalParametersToSystematicBasis(
      Eigen::ArrayXd const &pc_vals, Eigen::ArrayXd &systematic_vals) const;
  void TranformSystematicParametersProposalToBasis(
      Eigen::ArrayXd const &systematic_vals, Eigen::ArrayXd &pc_vals) const;

  constexpr static int ParameterInPCABlock = std::numeric_limits<int>::max();
  int GetProposalParameterIndexFromSystematicIndex(int i) const;

  void ConfigureAdaptive(YAML::Node const &, const std::string& matrix_name_str);

  Eigen::ArrayXd const &Propose();

  void Accept();

  void Reject();
};
