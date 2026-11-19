#pragma once

#include <vector>

#include "Samples/ReweightModel.h"

/// @brief The live events of a ReweightModel summed into groups, one per bin and signature
///
/// @details Inside one bin, events with the same signature differ only by their constant weights,
/// so their total is (sum of constants) times the signature's product. Built from the bin of every
/// live event; valid until the bins change, i.e. until the shifts move.
///
/// A group's events also share the sign of their constants: a group's weight is then positive
/// exactly when each of its events' weights is, which is the test the per event loop applies.
///
/// @author Henry Wallace
class EventGroups {
 public:
  /// @brief Group the live events
  /// @param Model The model the bins refer to
  /// @param LiveBin Global bin of every live event, negative if it contributes to none
  void Build(const ReweightModel& Model, const std::vector<int>& LiveBin);

  /// Bin b's groups are BinStart[b] to BinStart[b+1]-1
  std::vector<int> BinStart;
  /// Signature of every group
  std::vector<int> Sig;
  /// Sum of every group's constant weights
  std::vector<double> Const;
  /// Sum of their squares, for the sum of squared weights
  std::vector<double> Const2;
  /// Bins with at least one group, most groups first, the order threads take them in
  std::vector<int> BinOrder;

  /// @brief Number of groups
  int NGroups() const { return static_cast<int>(Sig.size()); }

 private:
  /// Live events ordered by bin, and within a bin by signature. Kept to save reallocating
  std::vector<int> Ordered;
  /// Next free slot of every bin while Ordered is filled
  std::vector<int> Next;
  /// @brief Add one run of a signature's events as up to two groups, split by the sign of the constant
  void AddRun(const ReweightModel& Model, int Begin, int End);
};
