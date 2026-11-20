#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "Samples/EventInfo.h"

/// @brief A flat description of a sample handler's events, built once, for the fast reweight backends
///
/// @details Every event weight is a product. Some factors are fixed for the whole fit (the ones an
/// experiment adds in AddAdditionalWeightPointers, when it promises they are constant); the rest
/// vary with the parameters (normalisations, splines, oscillation). The varying factors are
/// pointers, and many events share the very same ones. So:
///
/// - the **pool** is every distinct varying weight pointer, once;
/// - a **signature** is the list of pool entries one event multiplies, and events with identical
///   lists share a signature;
/// - each **live** event (one that some setting of the shifts could select) has a signature and
///   the product of its constant factors.
///
/// An event's weight is then its constant times its signature's product, and inside one bin all
/// events of a signature can be summed into one before any reweight (see EventGroups).
///
/// The layout is plain arrays, so the same model can be copied to a GPU unchanged.
///
/// @author Henry Wallace
struct ReweightModel {
  /// Number of events in the handler, live or not
  int NEvents = 0;
  /// Number of bins over all of the handler's samples
  int NBins = 0;

  /// Every distinct varying weight, once
  std::vector<const M3::float_t*> PoolPointer;

  /// Signature s multiplies pool entries SigPool[SigStart[s]] to SigPool[SigStart[s+1]-1], in the
  /// order the events list them
  std::vector<int> SigStart;
  /// See SigStart
  std::vector<int> SigPool;

  /// Handler index of every live event
  std::vector<int> Event;
  /// Signature of every live event
  std::vector<int> Sig;
  /// Product of every live event's constant weights
  std::vector<double> Const;
  /// Live events ordered by signature, so binning them in this order leaves each bin sorted by
  /// signature and its groups contiguous
  std::vector<int> BySignature;

  /// @brief Number of signatures
  int NSignatures() const { return static_cast<int>(SigStart.size()) - 1; }
  /// @brief Number of live events
  int NLive() const { return static_cast<int>(Event.size()); }
  /// @brief A one line summary for the log
  std::string Summary() const;
  /// @brief Memory the model holds, in bytes
  std::size_t Bytes() const {
    return sizeof(const M3::float_t*) * PoolPointer.size() + sizeof(int) * (SigStart.size() + SigPool.size() +
           Event.size() + Sig.size() + BySignature.size()) + sizeof(double) * Const.size();
  }
};

/// @brief Build a ReweightModel from a handler's events
/// @param Events The handler's events
/// @param ConstStart,ConstEnd Per event, the range of total_weight_pointers holding weights that are
/// constant for the whole fit
/// @param NeverSelected Per event, non-zero if a cut no shift can change rejects it
/// @param NBins Number of bins over all of the handler's samples
ReweightModel BuildReweightModel(const std::vector<EventInfo>& Events, const std::vector<std::uint16_t>& ConstStart,
                                 const std::vector<std::uint16_t>& ConstEnd, const std::vector<char>& NeverSelected,
                                 int NBins);
