#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Samples/OscillationProbability.h"
#include "Splines/WeightDerivativeSource.h"

namespace YAML { class Node; }

class OscillationGradientDevice;

/// @brief Derivatives of NuOscillator's oscillation weights with respect to the oscillation parameters
///
/// @details NuOscillator returns probabilities only. For its CUDAProb3Linear calculation (constant
/// density, straight baseline) this works out each weight again with the same formulae carried
/// through dual numbers (OscillationProbability.h), which gives the probability and its exact
/// derivative with respect to every oscillation parameter at once. Then the oscillation parameters
/// have a gradient, and a Hamiltonian sampler can move them like any other parameter.
///
/// Every distinct weight pointer is worked out once:
/// - **binned** NuOscillator: one weight per (energy bin, channel), evaluated at the bin centre,
///   which NuOscillator works out with the same ReturnBinCentersFromBinEdges used here. A few
///   thousand evaluations however many events there are;
/// - **unbinned**: one weight per event, at the event's own energy.
///
/// Create() refuses (returns null) unless the formulae reproduce NuOscillator's own weights to
/// rounding, see CompareWithNuOscillator(), so a configuration they don't describe is never
/// differentiated.
///
/// @author Henry Wallace
class OscillationGradient {
 public:
  /// @brief One event's oscillation weight and what it is the probability of
  struct EventWeight {
    /// NuOscillator's weight for the event
    const M3::float_t* Weight;
    /// The event's true energy, as NuOscillator was handed it
    double Energy;
    /// Whether it is an antineutrino
    bool Anti;
    /// Flavour before and after oscillation: 0 electron, 1 muon, 2 tau
    signed char From;
    /// See From
    signed char To;
  };

  /// @brief Set up, if the formulae apply
  /// @param Name The sample handler's name, for the log
  /// @param NuOscConfigFile NuOscillator's config
  /// @param ParIndex Index in the parameter handler of sin2_th12, sin2_th23, sin2_th13, dm2_12, dm2_23,
  /// delta_cp, path_length and matter_density
  /// @param ParValue Where each of those parameters' current value lives
  /// @param Free Which of the eight to differentiate with respect to (the others are fixed)
  /// @param Events Every event with a NuOscillator weight
  /// @return null, with the reason logged, if the formulae don't apply or don't reproduce NuOscillator
  static std::unique_ptr<OscillationGradient> Create(const std::string& Name, const std::string& NuOscConfigFile,
                                                     const std::vector<int>& ParIndex, const std::vector<const M3::float_t*>& ParValue,
                                                     const std::vector<bool>& Free, const std::vector<EventWeight>& Events);

  /// @brief Work out the derivatives at the current parameters, unless they haven't moved
  void Evaluate();
  /// @brief Where the derivatives live, for the gradient
  WeightDerivativeSource Source() const;
  /// @brief Largest relative difference between the formulae and NuOscillator's current weights
  double CompareWithNuOscillator() const;
  /// @brief Number of distinct weights
  std::size_t NWeights() const { return Weights.size(); }

  ~OscillationGradient();
  OscillationGradient(const OscillationGradient&) = delete;
  OscillationGradient& operator=(const OscillationGradient&) = delete;

  /// @brief Whether the derivatives are worked out on the GPU
  bool OnDevice() const { return Device != nullptr; }

 private:
  OscillationGradient();
  /// @brief The current values of the eight parameters
  void CurrentValues(double* Values) const;
  /// @brief The derivatives at constants C, on the CPU, into Derivatives
  void EvaluateOnHost(const M3::Osc::Constants<M3::Osc::OscDual>& C);
  /// @brief Largest difference between EvaluateOnHost's shared transition matrices and Probability one
  /// weight at a time, at constants C (zero: the two do the same operations in the same order)
  double CompareSharedWithSingle(const M3::Osc::Constants<M3::Osc::OscDual>& C);
  /// @brief Move the derivatives to the GPU if there is one and it agrees with the CPU (GPU builds only)
  void TryDevice(const std::string& Name);

  /// @name Steps of Create()
  /// @{
  /// @brief Whether NuOscillator's config is one the formulae describe: CUDAProb3Linear, binned or unbinned
  static bool SupportsConfig(const YAML::Node& Config, const std::string& Name);
  /// @brief Whether all eight CUDAProb3Linear parameters were found
  static bool HasAllParameters(const std::vector<int>& ParIndex, const std::vector<const M3::float_t*>& ParValue,
                               const std::vector<bool>& Free, const std::string& Name);
  /// @brief Keep one entry per distinct weight pointer, with the energy it is evaluated at
  /// @return false if events sharing a weight disagree on its energy or channel
  template <typename EnergyMap>
  bool CollectDistinctWeights(const std::vector<EventWeight>& Events, const EnergyMap& EvaluatedAt, const std::string& Name);
  /// @brief Record where the eight parameters live and which are differentiated with respect to
  void SetParameters(const std::vector<int>& ParIndex, const std::vector<const M3::float_t*>& ParValue, const std::vector<bool>& Free);
  /// @brief Fill SharedStart and SharedWeight
  void GroupBySharedTransition();
  /// @brief Compare the shared transition matrices with one weight at a time, and log the grouping
  void CheckSharedTransition(const std::string& Name);
  /// @}

  /// Works the derivatives out on the GPU, if set (GPU builds only). A shared_ptr, whose deleter is
  /// fixed where it is made, so a CPU build never needs the class's definition
  std::shared_ptr<OscillationGradientDevice> Device;
  /// The GPU's results before they are stored as M3::float_t
  std::vector<double> DeviceOut;

  /// Every distinct weight pointer, sorted
  std::vector<const M3::float_t*> Weights;
  /// Energy each weight is evaluated at
  std::vector<double> Energy;
  /// Antineutrino or not, per weight
  std::vector<char> Anti;
  /// Flavours before and after, per weight
  std::vector<signed char> From;
  /// See From
  std::vector<signed char> To;
  /// Weights sharing an energy and neutrino type, whose channels come from one transition matrix:
  /// group g is weights SharedWeight[SharedStart[g]] onwards
  std::vector<int> SharedStart;
  /// See SharedStart
  std::vector<int> SharedWeight;
  /// Where the eight parameters' values live
  std::vector<const M3::float_t*> ParValue;
  /// Which of the eight are differentiated with respect to, as indices 0..7
  std::vector<int> FreeSlot;
  /// Their indices in the parameter handler, see WeightDerivativeSource::MultiParam
  std::vector<int> FreeParam;
  /// d(weight)/d(free parameter j) at weight * FreeSlot.size() + j
  std::vector<M3::float_t> Derivatives;
  /// Parameter values the derivatives are for
  std::vector<double> LastValues;
};
