#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Walloca"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wduplicated-branches"
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include "Eigen/Dense"
#pragma GCC diagnostic pop

#include <cmath>
#include <vector>

Eigen::MatrixXd CalculateTruncatedPCARotation(Eigen::MatrixXd mx_phyiscs_basis,
                                              double threshold) {

  // a covariance is real symmetric, so self adjoint
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EigenSolver(mx_phyiscs_basis);

  if (EigenSolver.info() != Eigen::Success) {
    MACH3LOG_ERROR("Failed to eigen decompose matrix.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  Eigen::VectorXd eigen_val = EigenSolver.eigenvalues();
  Eigen::MatrixXd eigen_vect = EigenSolver.eigenvectors();
  std::vector<std::pair<int, double>> evals;
  for (int i = 0; i < eigen_val.size(); ++i) {
    evals.emplace_back(i, eigen_val[i]);
  }

  // sorting is not strictly neccessary but orders the orthogonal parameter
  //   basis approximately by uncertainty.
  std::sort(
      evals.begin(), evals.end(),
      [](std::pair<int, double> const &l, std::pair<int, double> const &r) {
        return l.second < r.second;
      });

  double evalsum = 0;
  int northo_parameters = 0;

  for (auto &[i, ev] : evals) {
    if (std::fabs(ev) > threshold) {
      evalsum += std::fabs(ev);
      northo_parameters++;
    } else {
      break;
    }
  }

  // rows of this matrix correspond to eigenvectors of the input matrix
  //   we can go from parameters defined in the orthogonal basis back to the
  //  'physics' parameters by RowVectOfOrthoParamVals * pca.ortho_to_physics
  Eigen::MatrixXd ortho_to_physics =
      Eigen::MatrixXd::Zero(northo_parameters, mx_phyiscs_basis.rows());

  for (int i = 0; i < northo_parameters; ++i) {
    for (int j = 0; j < mx_phyiscs_basis.rows(); ++j) {
      ortho_to_physics.row(i)[j] =
          eigen_vect.col(evals[i].first)[evals[j].first] *
          std::sqrt(evals[i].second);
    }
  }

  return ortho_to_physics;
}
