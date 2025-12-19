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
    size_t nsteps;
    Eigen::MatrixXd accepted_parameter_covariance;
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

  void EnableAdaption(YAML::Node const &) {}
  // void SetAdaptionCovariance(Eigen::MatrixXd parameter_covariance,
  //                            size_t nsteps) {
  //   adaptive.nsteps = nsteps;
  //   adaptive.accepted_parameter_covariance = parameter_covariance;
  // }

  Eigen::ArrayXd const &Propose();
  void Accept() { proposal_basis.current = proposal_basis.proposed; }
};
