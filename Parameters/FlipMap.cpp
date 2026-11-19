#include "Parameters/FlipMap.h"

#include <algorithm>
#include <cmath>
#include <limits>

// ************************************************
FlipMap::FlipMap(const FlipGroup& Group_, const std::function<bool(int)>& IsFixed) : Group(Group_) {
// ************************************************
  for (const int index : Group.FlipParameterIndex) {
    UseReflection.push_back(IsFixed(index) ? 0 : 1);
    if (UseReflection.back()) ChangedIndex.push_back(index);
  }
  for (const auto& functional_flip : Group.FunctionalFlipParameters) {
    UseFunctional.push_back(IsFixed(functional_flip.target_index) ? 0 : 1);
    if (UseFunctional.back()) ChangedIndex.push_back(functional_flip.target_index);
  }
}

// ************************************************
std::vector<double> FlipMap::Apply(const std::vector<double>& X) const {
// ************************************************
  // As in ParameterHandlerBase::FlipParameterGroup(), every functional flip reads the values from
  // before any flip
  std::vector<double> Out(ChangedIndex.size());
  std::size_t k = 0;
  for (std::size_t i = 0; i < Group.FlipParameterIndex.size(); ++i) {
    if (!UseReflection[i]) continue;
    Out[k++] = 2. * Group.FlipParameterPoint[i] - X[static_cast<std::size_t>(Group.FlipParameterIndex[i])];
  }
  for (std::size_t f = 0; f < Group.FunctionalFlipParameters.size(); ++f) {
    if (!UseFunctional[f]) continue;
    const auto& functional_flip = Group.FunctionalFlipParameters[f];
    std::vector<double> Args(functional_flip.argument_indices.size());
    for (std::size_t a = 0; a < Args.size(); ++a) Args[a] = X[static_cast<std::size_t>(functional_flip.argument_indices[a])];
    const double Target = X[static_cast<std::size_t>(functional_flip.target_index)];
    Out[k++] = functional_flip.evaluator->EvalPar(&Target, Args.data());
  }
  return Out;
}

// ************************************************
double FlipMap::LogAbsJacobian(const std::vector<double>& X0) const {
// ************************************************
  const std::size_t K = ChangedIndex.size();
  std::vector<double> J(K * K);
  std::vector<double> X = X0;
  for (std::size_t c = 0; c < K; ++c) {
    const std::size_t index = static_cast<std::size_t>(ChangedIndex[c]);
    const double Centre = X[index];
    const double H = 1e-6 * std::max(std::fabs(Centre), 1e-3);
    X[index] = Centre + H;
    const std::vector<double> Up = Apply(X);
    X[index] = Centre - H;
    const std::vector<double> Down = Apply(X);
    X[index] = Centre;
    for (std::size_t r = 0; r < K; ++r) J[r * K + c] = (Up[r] - Down[r]) / (2. * H);
  }
  return LogAbsDeterminant(std::move(J), K);
}

// ************************************************
double LogAbsDeterminant(std::vector<double> A, const std::size_t K) {
// ************************************************
  double LogDet = 0.;
  for (std::size_t c = 0; c < K; ++c) {
    std::size_t Pivot = c;
    for (std::size_t r = c + 1; r < K; ++r) if (std::fabs(A[r * K + c]) > std::fabs(A[Pivot * K + c])) Pivot = r;
    if (!(std::fabs(A[Pivot * K + c]) > 0.)) return std::numeric_limits<double>::quiet_NaN();
    if (Pivot != c) for (std::size_t k = 0; k < K; ++k) std::swap(A[c * K + k], A[Pivot * K + k]);
    LogDet += std::log(std::fabs(A[c * K + c]));
    for (std::size_t r = c + 1; r < K; ++r) {
      const double F = A[r * K + c] / A[c * K + c];
      for (std::size_t k = c; k < K; ++k) A[r * K + k] -= F * A[c * K + k];
    }
  }
  return LogDet;
}
