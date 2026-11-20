#pragma once

#include <string>
#include <vector>

#include "Samples/EventInfo.h"
#include "Splines/WeightDerivativeSource.h"

/// @brief The gradient of a sample handler's likelihood with respect to its parameters
///
/// @details For a binned likelihood \f$-\log\mathcal{L} = \sum_b \ell(D_b, \mu_b)\f$ with
/// \f$\mu_b = \sum_{e \in b} W_e\f$ and every event weight a product \f$W_e = \prod_k w_{ek}\f$,
/// \f[
///   \frac{\partial(-\log\mathcal{L})}{\partial x_p}
///     = \sum_e \frac{\partial\ell}{\partial\mu}\Big|_{b(e)} \, W_e
///       \sum_{k:\,w_{ek}\text{ depends on } x_p} \frac{1}{w_{ek}}\frac{\partial w_{ek}}{\partial x_p}.
/// \f]
/// So the gradient is one pass over the events after the reweight, given each event's bin and total
/// weight and the per bin factor \f$\partial\ell/\partial\mu\f$. It needs no array of bin
/// contents per parameter, and nothing per event weight: a weight's derivative is found from its
/// address, see WeightDerivativeSource.
///
/// The reweight records each event's bin and weight with Record() and MarkUnused(); Compute() then
/// runs the pass.
///
/// @author Henry Wallace
class LikelihoodGradient {
 public:
  /// Bin recorded for an event that contributes nothing: not selected, zero weight, or outside the binning
  static constexpr int kNoBin = -1;

  /// @brief Prepare for a handler's events
  /// @param Sources Every block of weights whose derivatives are known
  /// @param NPars Number of parameters in the handler the gradient is indexed by
  /// @param Events The handler's events, to count how many of their weights have a derivative
  LikelihoodGradient(std::vector<WeightDerivativeSource> Sources, int NPars, const std::vector<EventInfo>& Events);

  /// @brief Record an event's bin and total weight from the reweight
  /// @details Only call with Weight > 0, which guarantees none of its weights is zero
  void Record(const unsigned int iEvent, const int Bin, const double Weight) {
    EventBin[iEvent] = Bin;
    EventWeight[iEvent] = Weight;
  }
  /// @brief Record that an event contributes nothing to the last reweight
  void MarkUnused(const unsigned int iEvent) { EventBin[iEvent] = kNoBin; }

  /// @brief Run the pass
  /// @param Events The handler's events, as given to the constructor
  /// @param BinFactor \f$\partial\ell/\partial\mu_b\f$ for every bin
  /// @param Out Resized to one entry per parameter and overwritten with
  /// \f$\partial(-\log\mathcal{L})/\partial x_p\f$
  /// @details Takes a buffer rather than returning one: this runs on every leapfrog step, so a
  /// fresh vector per call is allocation traffic on the hottest path there is
  void Compute(const std::vector<EventInfo>& Events, const std::vector<double>& BinFactor,
               std::vector<double>& Out) const;

  /// @brief A one line summary of how the weights were classified, for the log
  std::string Summary() const;

 private:
  /// See the constructor
  std::vector<WeightDerivativeSource> Sources;
  /// Number of parameters the gradient is indexed by
  int NPars;
  /// Bin of each event in the last reweight, kNoBin if it contributes nothing
  std::vector<int> EventBin;
  /// Total weight of each event in the last reweight, meaningful only where EventBin isn't kNoBin
  std::vector<double> EventWeight;
  /// Event weights found in a source, and not, when the constructor looked
  std::size_t NResolved = 0;
  /// See NResolved
  std::size_t NUnresolved = 0;
};
