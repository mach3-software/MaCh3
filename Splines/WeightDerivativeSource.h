#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

#include "Manager/Core.h"

/// @brief Where the derivatives of a block of event weights live
///
/// @details Events hold plain pointers to their weights (EventInfo::total_weight_pointers). Many of
/// those pointers point into storage owned by something that can also give the weight's derivatives:
/// the parameter values (a normalisation weight is the parameter itself), a spline handler's
/// evaluated weights, or an oscillation calculation's probabilities. A source describes one such set
/// of weights, so a weight's derivatives can be found from its address alone, with nothing stored per
/// event.
///
/// Which weights belong to the source:
/// - a contiguous block: weight k is Weights[k], for k < N; or
/// - when Listed isn't empty, exactly the addresses in it (sorted): weight k is Listed[k].
///
/// What each weight depends on:
/// - **one parameter** (MultiParam empty): weight k's parameter is GlobalParam[ParamOfWeight[k]], or
///   simply k when ParamOfWeight is null, and its derivative is Derivatives[k], or 1 when
///   Derivatives is null (a weight that is the parameter value itself);
/// - **several parameters** (MultiParam not empty): every weight depends on the parameters in
///   MultiParam, and its derivative with respect to MultiParam[j] is
///   Derivatives[k * MultiParam.size() + j].
///
/// A weight outside every source (an experiment specific weight, the constant 1) has no derivative.
///
/// @author Henry Wallace
struct WeightDerivativeSource {
  /// First weight of a contiguous block
  const M3::float_t* Weights = nullptr;
  /// Number of weights in a contiguous block
  std::size_t N = 0;
  /// The weights' addresses, sorted, when they aren't one contiguous block
  std::vector<const M3::float_t*> Listed;
  /// The derivatives, laid out as the class description says. Null means every derivative is 1
  const M3::float_t* Derivatives = nullptr;
  /// One parameter form: per weight, the local number of the parameter it depends on. Null means
  /// weight k depends on parameter k
  const short int* ParamOfWeight = nullptr;
  /// One parameter form: local parameter number to index in the parameter handler
  std::vector<int> GlobalParam;
  /// Several parameter form: the parameters every weight depends on
  std::vector<int> MultiParam;

  /// @brief Whether this describes anything
  bool IsValid() const { return (Weights != nullptr && N > 0) || !Listed.empty(); }

  /// @brief Position of a weight in the source, or -1 if it isn't one of its weights
  long IndexOf(const M3::float_t* Weight) const {
    if (!Listed.empty()) {
      const auto Found = std::lower_bound(Listed.begin(), Listed.end(), Weight, std::less<const M3::float_t*>());
      return (Found != Listed.end() && *Found == Weight) ? static_cast<long>(Found - Listed.begin()) : -1;
    }
    // std::less gives a total order over pointers into unrelated arrays, which < doesn't promise
    if (std::less<const M3::float_t*>()(Weight, Weights) || !std::less<const M3::float_t*>()(Weight, Weights + N)) return -1;
    return static_cast<long>(Weight - Weights);
  }
};

/// @brief Visit every derivative a weight has
/// @param Sources Every set of weights whose derivatives are known
/// @param Weight Address of the weight
/// @param Visit Called as Visit(Param, Derivative) for each parameter the weight depends on, with
/// the address of its derivative, or nullptr for a derivative of 1
/// @return false if no source holds the weight, in which case it has no derivatives
template <typename Visitor>
inline bool VisitWeightDerivatives(const std::vector<WeightDerivativeSource>& Sources, const M3::float_t* Weight,
                                   Visitor&& Visit) {
  for (const auto& Source : Sources) {
    const long Index = Source.IndexOf(Weight);
    if (Index < 0) continue;
    const std::size_t k = static_cast<std::size_t>(Index);
    if (!Source.MultiParam.empty()) {
      const std::size_t NPars = Source.MultiParam.size();
      for (std::size_t j = 0; j < NPars; ++j) Visit(Source.MultiParam[j], Source.Derivatives + k * NPars + j);
      return true;
    }
    const int Param = (Source.ParamOfWeight == nullptr) ? static_cast<int>(k)
                                                        : Source.GlobalParam[static_cast<std::size_t>(Source.ParamOfWeight[k])];
    Visit(Param, (Source.Derivatives == nullptr) ? nullptr : Source.Derivatives + k);
    return true;
  }
  return false;
}
