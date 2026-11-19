#pragma once

// SpecialProposals.h uses YAML::Node without including it
#include "Manager/YamlHelper.h"
#include "Parameters/SpecialProposals.h"

#include <functional>
#include <vector>

/// @brief The deterministic map one flip group applies, for proposing it as a move of its own
///
/// @details The group's plain reflections send \f$x \to 2p - x\f$ about their FlipParameter points;
/// its functional flips set their target to a formula of the values from before any flip. Fixed
/// parameters are left out. ParameterHandlerBase::ProposeFlips() uses this to propose the flip,
/// check the map is its own inverse (without which the move breaks detailed balance), and find
/// \f$\log|\det J|\f$, which a functional flip's acceptance needs since it need not preserve volume.
///
/// @author Henry Wallace
class FlipMap {
 public:
  /// @param Group The flip group
  /// @param IsFixed Whether a parameter is fixed, and so left out
  FlipMap(const FlipGroup& Group, const std::function<bool(int)>& IsFixed);

  /// @brief The parameters the map changes: the plain reflections, then the functional flips' targets
  const std::vector<int>& Changed() const { return ChangedIndex; }

  /// @brief New values of the changed parameters, in the order of Changed()
  /// @param X A full parameter vector
  std::vector<double> Apply(const std::vector<double>& X) const;

  /// @brief \f$\log|\det J|\f$ of the map over the parameters it changes, at X
  /// @details By central differences. A plain reflection contributes exactly -1 on its diagonal; a
  /// functional flip, say a theta_13 that keeps \f$\sin^2\theta_{23}\sin^2 2\theta_{13}\f$ fixed
  /// while theta_23 is reflected, does not preserve volume
  /// @return NaN if the Jacobian is singular
  double LogAbsJacobian(const std::vector<double>& X) const;

 private:
  /// The group
  const FlipGroup& Group;
  /// Per plain reflection in the group, whether it is used (its parameter isn't fixed)
  std::vector<char> UseReflection;
  /// Per functional flip in the group, whether it is used
  std::vector<char> UseFunctional;
  /// See Changed()
  std::vector<int> ChangedIndex;
};

/// @brief \f$\log|\det A|\f$ of a K by K row-major matrix, by Gaussian elimination with partial pivoting
/// @return NaN if the matrix is singular
double LogAbsDeterminant(std::vector<double> A, std::size_t K);
