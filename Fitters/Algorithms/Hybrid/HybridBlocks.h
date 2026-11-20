#pragma once

#include "Algorithms/Hybrid/HybridSettings.h"

#include <vector>

class ParameterHandlerBase;
class SampleHandlerInterface;

/// @brief Which parameters HybridMCMC's random walk owns, and which NUTS owns
///
/// @details A parameter goes to the random walk if its name matches RandomWalkNames, or if it is of a
/// type in RandomWalkTypes; without RandomWalkTypes, if some sample has no gradient for it
/// (SampleHandlerInterface::GetGradientCoverage). Everything else goes to NUTS. Fixed parameters
/// sit in one mask or the other but are counted in neither block.
///
/// The masks are what ParameterHandlerBase::SetProposalMask takes: per handler, per parameter,
/// non-zero for the parameters the kernel moves.
///
/// @author Henry Wallace
class HybridBlocks {
 public:
  /// @brief Assign every parameter to a block and log the split
  /// @param Systematics The fitter's parameter handlers
  /// @param Samples The fitter's samples, which already know their gradient coverage
  /// @param Settings Which types and names go to the random walk
  void Build(const std::vector<ParameterHandlerBase*>& Systematics,
             const std::vector<SampleHandlerInterface*>& Samples, const HybridSettings& Settings);

  /// @brief Whether parameter i of handler s belongs to the random walk
  bool InRandomWalk(const size_t s, const int i) const { return RandomWalkMask[s][static_cast<size_t>(i)] != 0; }
  /// @brief Whether handler s holds any random walk parameter
  bool HandlerHasRandomWalk(const size_t s) const;
  /// @brief Proposal mask of handler s for one of the two kernels
  /// @param RandomWalk True for the random walk's block, false for NUTS's
  const std::vector<char>& Mask(const size_t s, const bool RandomWalk) const {
    return RandomWalk ? RandomWalkMask[s] : NutsMask[s];
  }

  /// @brief Free parameters the random walk moves
  int NRandomWalk() const { return NRandomWalk_; }
  /// @brief Free parameters NUTS moves
  int NNuts() const { return NNuts_; }

  /// @brief Throw unless both blocks have something to move
  /// @param RandomWalkProb The mixture's random walk probability, for the message
  void RequireBothBlocks(double RandomWalkProb) const;
  /// @brief Log how many random walk parameters each handler gives, and (debug) which
  void LogRandomWalkMembers(const std::vector<ParameterHandlerBase*>& Systematics) const;

 private:
  /// @brief Whether a parameter belongs to the random walk, by name, type or gradient coverage
  bool BelongsToRandomWalk(const ParameterHandlerBase& Handler, size_t s, int i, const HybridSettings& Settings) const;
  /// @brief Per handler, per parameter: whether some sample has no gradient for it
  void FindParametersWithoutGradient(const std::vector<ParameterHandlerBase*>& Systematics,
                                     const std::vector<SampleHandlerInterface*>& Samples);

  /// Per handler, per parameter: non-zero for the random walk's parameters
  std::vector<std::vector<char>> RandomWalkMask;
  /// The complement of #RandomWalkMask: NUTS's parameters
  std::vector<std::vector<char>> NutsMask;
  /// Per handler, per parameter: whether some sample can't differentiate with respect to it.
  /// Only filled when the blocks follow gradient coverage
  std::vector<std::vector<char>> NoGradient;
  /// See NRandomWalk()
  int NRandomWalk_ = 0;
  /// See NNuts()
  int NNuts_ = 0;
  /// Fixed parameters, in neither block's count
  int NFixed = 0;
};
