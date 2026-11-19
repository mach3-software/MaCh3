#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Manager/Core.h"

_MaCh3_Safe_Include_Start_
#include "TF1.h"
_MaCh3_Safe_Include_End_

/// @brief Configuration for a formula-based functional flip.
/// @author David Riley
struct FunctionalFlipProposal {
  int target_index = M3::_BAD_INT_;
  /// parameter indices passed to evaluator
  std::vector<int> argument_indices;
  /// Parameter names that are passed to evaluator
  std::vector<std::string> argument_names;
  /// function used to evaluate flipping formula
  std::unique_ptr<TF1> evaluator;
  /// @brief get original formula used for evaluator
  std::string GetFormula() const {return evaluator->GetTitle();}
  FunctionalFlipProposal() = default;
  FunctionalFlipProposal(FunctionalFlipProposal&&) _noexcept_ = default;
  FunctionalFlipProposal& operator=(FunctionalFlipProposal&&) _noexcept_ = default;

  FunctionalFlipProposal(const FunctionalFlipProposal&) = delete;
  FunctionalFlipProposal& operator=(const FunctionalFlipProposal&) = delete;
};

/// @brief Configuration for a pending functional flip awaiting all parameters to be initialised
/// @author David Riley
struct PendingFunctionalFlipProposal {
  int target_index = M3::_BAD_INT_;
  std::string group_name;
  YAML::Node config;
};


/// @brief Struct to hold information about a group of parameters that flip together at the same time
/// @author Charlotte Knight
/// @author Liban Warsame
struct FlipGroup {
   /// Indices of parameters with flip symmetry
   std::vector<int> FlipParameterIndex;  
   /// Central points around which parameters are flipped
   std::vector<double> FlipParameterPoint; 
   /// Formula-driven flips that should be applied when this group flips.
   std::vector<FunctionalFlipProposal> FunctionalFlipParameters;
};
