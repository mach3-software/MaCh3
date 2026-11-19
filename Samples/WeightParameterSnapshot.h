#pragma once

#include "Parameters/ParameterHandlerGeneric.h"

#include <vector>

/// @brief The values of the parameters a sample's oscillation and spline weights depend on, when
/// those weights were last evaluated
///
/// @details Those weights depend only on the non functional parameters. A move of the shifts alone
/// (every random walk sweep of the hybrid sampler) leaves them as they were, so with a reweight
/// backend in use SampleHandlerBase::EvaluateWeights() skips them when Matches() says so.
///
/// @author Henry Wallace
class WeightParameterSnapshot {
 public:
  /// @brief Whether every watched parameter still has the value it had at the last Take()
  bool Matches(const ParameterHandlerGeneric& Handler) const {
    if (!Valid) return false;
    for (std::size_t k = 0; k < Params.size(); ++k) {
      if (Handler.GetParProp(Params[k]) != Values[k]) return false;
    }
    return true;
  }
  /// @brief Record the current values, after the weights have been evaluated at them
  void Take(const ParameterHandlerGeneric& Handler) {
    if (Params.empty()) {
      for (int i = 0; i < Handler.GetNumParams(); ++i) {
        if (Handler.GetParamType(i) != SystType::kFunc) Params.push_back(i);
      }
    }
    Values.resize(Params.size());
    for (std::size_t k = 0; k < Params.size(); ++k) Values[k] = Handler.GetParProp(Params[k]);
    Valid = true;
  }
  /// @brief Forget the snapshot, so the next evaluation can't be skipped
  void Invalidate() { Valid = false; }

 private:
  /// Indices of the watched parameters: every non functional one
  std::vector<int> Params;
  /// Their values at the last Take()
  std::vector<M3::float_t> Values;
  /// Whether #Values describes the weights currently held
  bool Valid = false;
};
