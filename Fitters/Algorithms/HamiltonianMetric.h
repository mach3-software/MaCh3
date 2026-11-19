#pragma once

#include <string>
#include <vector>

#include "TMatrixDSym.h"

/// @brief One parameter handler's block of the metric a Hamiltonian sampler integrates against
///
/// @details The sampler carries a velocity \f$w = M^{-1} p\f$ rather than a momentum, so every
/// operation it needs is a product with \f$M^{-1}\f$ or with its Cholesky factor \f$L\f$,
/// \f$L L^T = M^{-1}\f$:
/// - a fresh velocity is \f$w = L z\f$ for standard normal \f$z\f$ (ThrowVelocity);
/// - half a leapfrog kick is \f$w \rightarrow w - \frac{h}{2} M^{-1}\nabla U\f$ (Kick);
/// - the kinetic energy \f$\frac{1}{2} w^T M w\f$ is \f$\frac{1}{2}|u|^2\f$ with \f$L u = w\f$, and
///   the momentum \f$p = M w\f$ follows from \f$L^T p = u\f$ (KineticEnergy...).
///
/// \f$M^{-1}\f$ is a handler's throw matrix, optionally shrunk towards its diagonal, with the rows
/// and columns of frozen parameters (fixed, or owned by another kernel) decoupled and given a unit
/// diagonal. Products skip the exact zeros outside each row's first and last non-zero column (its
/// skyline), which leaves every sum unchanged; a fully shrunk block uses O(N) loops throughout.
///
/// @author Henry Wallace
class MetricBlock {
 public:
  /// @brief Build the block
  /// @param Throw The handler's throw matrix, \f$M^{-1}\f$ before shrinkage
  /// @param Frozen Per parameter, non-zero if it must never move
  /// @param Shrinkage \f$\lambda\f$: off-diagonal elements are scaled by \f$1-\lambda\f$
  /// @param PriorWidth Per parameter, the width to fall back on when a diagonal block has no usable variance
  /// @param Name The handler's name, for the log
  void Build(const TMatrixDSym& Throw, const std::vector<char>& Frozen, double Shrinkage,
             const std::vector<double>& PriorWidth, const std::string& Name);

  /// @brief Number of parameters in the block
  int NPars() const { return N; }
  /// @brief Whether the block is diagonal
  bool IsDiagonal() const { return Diagonal; }

  /// @brief \f$w = L z\f$
  /// @param Throws N standard normal throws
  /// @param Velocity N entries to fill
  void ThrowVelocity(const double* Throws, double* Velocity) const;
  /// @brief \f$w \rightarrow w - \frac{h}{2} M^{-1}\nabla U\f$
  /// @param Velocity N entries to update
  /// @param Gradient N entries of \f$\nabla U\f$
  /// @param Step \f$h\f$, negative to integrate backwards
  void Kick(double* Velocity, const double* Gradient, double Step) const;
  /// @brief Add this block's \f$\frac{1}{2} w^T M w\f$ to a running total
  /// @param Velocity N entries
  /// @param Scratch N entries of scratch space
  /// @param Total The kinetic energy of the blocks before this one
  /// @return Total plus this block's, added element by element so the sum is the same over any
  /// split into blocks
  double AddKineticEnergy(const double* Velocity, double* Scratch, double Total) const;
  /// @brief AddKineticEnergy(), and the momentum \f$p = M w\f$
  /// @param Velocity N entries
  /// @param Scratch N entries of scratch space
  /// @param Momentum N entries to fill
  /// @param Total The kinetic energy of the blocks before this one
  double AddKineticEnergyAndMomentum(const double* Velocity, double* Scratch, double* Momentum, double Total) const;
  /// @brief Reflect off a wall normal to parameter i, conserving the kinetic energy
  /// @details Flips the momentum along the wall's normal, \f$p \rightarrow p - 2 (w_i / M^{-1}_{ii}) e_i\f$,
  /// so \f$w \rightarrow w - 2 (w_i / M^{-1}_{ii}) M^{-1} e_i\f$: \f$w_i\f$ changes sign and the parameters
  /// correlated with i change with it, which keeps \f$\frac{1}{2} w^T M w\f$ exactly. Flipping \f$w_i\f$ alone
  /// is only right for a diagonal metric; with a dense one it changes the kinetic energy and the
  /// trajectory diverges. A linear involution, so reversible and volume preserving
  /// @param Velocity N entries to update
  /// @param i The parameter at its wall
  void Reflect(double* Velocity, int i) const;

 private:
  /// @brief Find each row's first and last non-zero columns, see the class description
  void BuildSkyline();

  /// Number of parameters
  int N = 0;
  /// Whether every off-diagonal element is zero, so every loop collapses to O(N)
  bool Diagonal = false;
  /// \f$M^{-1}\f$, row major
  std::vector<double> InvMass;
  /// \f$L\f$, lower triangular, row major
  std::vector<double> Chol;
  /// In row i of Chol, the first column that isn't exactly zero
  std::vector<int> CholRowFirst;
  /// In column i of Chol, the last row that isn't exactly zero
  std::vector<int> CholColLast;
  /// In row i of InvMass, the first column that isn't exactly zero
  std::vector<int> InvRowFirst;
  /// In row i of InvMass, the last column that isn't exactly zero
  std::vector<int> InvRowLast;
};
