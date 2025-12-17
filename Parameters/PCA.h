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

inline std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
CalculateTruncatedPCARotation(Eigen::MatrixXd mx_syst_basis, double threshold) {

  int nsyst_params = int(mx_syst_basis.rows());

  // a covariance is real symmetric, so self adjoint
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EigenSolver(mx_syst_basis);

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

  // sorting is not strictly neccessary but orders the pc parameter
  //   basis approximately by uncertainty.
  std::sort(
      evals.rbegin(), evals.rend(),
      [](std::pair<int, double> const &l, std::pair<int, double> const &r) {
        return l.second < r.second;
      });

  double evalsum = 0;
  int npc_parameters = 0;

  double trace = mx_syst_basis.trace();

  for (auto &[i, ev] : evals) {
    if ((ev / trace) > threshold) {
      evalsum += ev;
      npc_parameters++;
    } else {
      break;
    }
  }

  // cols of this matrix correspond to eigenvectors of the input matrix
  //   we can rotate from the systematic basis to the pc basis by:
  //     TruncatedEigenVectorMatrix * ColVecPCParams = ColVecSystParams
  Eigen::MatrixXd TruncatedEigenVectorMatrix =
      Eigen::MatrixXd::Zero(nsyst_params, npc_parameters);
  Eigen::MatrixXd TruncatedsqrtEigenValueMatrix =
      Eigen::MatrixXd::Zero(npc_parameters, npc_parameters);

  for (int i = 0; i < npc_parameters; ++i) {
    TruncatedEigenVectorMatrix.col(i) = eigen_vect.col(evals[i].first);

    // bias eigen vectors towards being in the positive direction relative to
    // other parameters
    if (TruncatedEigenVectorMatrix.col(i).sum() < 0) {
      TruncatedEigenVectorMatrix.col(i).array() *= -1;
    }

    TruncatedsqrtEigenValueMatrix(i, i) = std::sqrt(evals[i].second);
  }

  return {TruncatedEigenVectorMatrix * TruncatedsqrtEigenValueMatrix,
          TruncatedsqrtEigenValueMatrix.inverse() *
              TruncatedEigenVectorMatrix.transpose()};
}
