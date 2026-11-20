#pragma once

#include <algorithm>
#include <cmath>

/// @brief Book keeping and proposal width tuning for the hybrid sampler's random walk kernel
///
/// @details Counts every random walk proposal (accepted, rejected on the likelihood, or thrown out
/// before the likelihood because something was out of bounds) and the mean likelihood penalty of
/// the evaluated ones, overall and per batch of steps. Tune() walks the proposal width towards a
/// target acceptance, a batch at a time:
///
/// - normally by Robbins-Monro on the log of the width, \f$\log s \leftarrow \log s + \gamma_n
///   (\alpha_n - \alpha^*)\f$ with \f$\gamma_n = 2/\sqrt{n}\f$, capped at a factor of two;
/// - but when a batch accepted almost nothing with a large mean penalty, the acceptance carries no
///   magnitude, so the width jumps instead by \f$\Delta\log L^{-1/p}\f$ for \f$\Delta\log L
///   \propto s^p\f$ (at most a factor of 100). That escapes a stalled block in an update or two.
///
/// Pure arithmetic, with no MaCh3 dependencies.
///
/// @author Henry Wallace
class RandomWalkTuner {
 public:
  /// @brief What one batch of proposals did
  struct Batch {
    /// Proposals in the batch
    unsigned int Steps = 0;
    /// Fraction accepted
    double Rate = 0.;
    /// Fraction thrown out on bounds before the likelihood was evaluated
    double OutOfBounds = 0.;
    /// Number that reached the Metropolis test
    unsigned int NEvaluated = 0;
    /// Their mean \f$\Delta(-\log L)\f$
    double MeanDLogL = 0.;
  };
  /// @brief How Tune() changed the width
  enum class Update { Ordinary, Jump };

  /// @param InitialFactor Starting multiplier of the optimal \f$2.38/\sqrt{d}\f$ width
  /// @param TargetAcceptance The acceptance Tune() steers towards
  /// @param PenaltyExponent \f$p\f$ in \f$\Delta\log L \propto s^p\f$, see the class description
  RandomWalkTuner(const double InitialFactor, const double TargetAcceptance, const double PenaltyExponent)
      : Factor(InitialFactor), Target(TargetAcceptance), Exponent(PenaltyExponent) {}

  /// @brief Record a proposal thrown out on bounds, before any likelihood was evaluated
  void RecordOutOfBounds() { ++Steps; ++OutOfBounds; }
  /// @brief Record a proposal that reached the Metropolis test
  /// @param DeltaLogL Its \f$\Delta(-\log L)\f$
  /// @param Accepted Whether it was accepted
  void RecordEvaluated(const double DeltaLogL, const bool Accepted) {
    ++Steps;
    DLogLSum += DeltaLogL;
    ++NEvaluated;
    if (Accepted) ++Accepted_;
  }

  /// @brief The statistics of the proposals since the last call, which starts the next batch
  Batch CloseBatch() {
    Batch B;
    B.Steps = Steps - Mark.Steps;
    if (B.Steps == 0) return B;
    B.Rate = static_cast<double>(Accepted_ - Mark.Accepted) / static_cast<double>(B.Steps);
    B.OutOfBounds = static_cast<double>(OutOfBounds - Mark.OutOfBounds) / static_cast<double>(B.Steps);
    B.NEvaluated = NEvaluated - Mark.NEvaluated;
    B.MeanDLogL = (B.NEvaluated > 0) ? (DLogLSum - Mark.DLogLSum) / static_cast<double>(B.NEvaluated) : 0.;
    Mark = {Steps, Accepted_, OutOfBounds, NEvaluated, DLogLSum};
    return B;
  }

  /// @brief Move the width factor according to a batch, see the class description
  /// @return Whether it was an ordinary update or a jump out of a stall
  Update Tune(const Batch& B) {
    ++Updates;
    Update Kind = Update::Ordinary;
    if (B.Rate < 0.01 && B.NEvaluated > 0 && std::isfinite(B.MeanDLogL) && B.MeanDLogL > 2.) {
      // A batch which accepted nothing says the width is wrong but not by how much, while its mean
      // penalty does. One jump is still bounded, so a pathological batch can't send it anywhere absurd
      LastLogStep = std::max(-std::log(100.), -std::log(B.MeanDLogL) / Exponent);
      Kind = Update::Jump;
    } else {
      // The gain has to decay for the sequence to settle rather than ring around the target
      const double Gain = 2. / std::sqrt(static_cast<double>(Updates));
      const double MaxLogStep = std::log(2.);
      LastLogStep = std::max(-MaxLogStep, std::min(MaxLogStep, Gain * (B.Rate - Target)));
    }
    Factor *= std::exp(LastLogStep);
    return Kind;
  }

  /// @brief The current width factor
  double GetFactor() const { return Factor; }
  /// @brief The factor the last Tune() multiplied the width by
  double GetLastChange() const { return std::exp(LastLogStep); }
  /// @brief Number of Tune() calls
  int GetUpdates() const { return Updates; }
  /// @brief The acceptance Tune() steers towards
  double GetTarget() const { return Target; }

  /// @name Totals over the whole chain
  /// @{
  unsigned int GetSteps() const { return Steps; }
  unsigned int GetAccepted() const { return Accepted_; }
  unsigned int GetOutOfBounds() const { return OutOfBounds; }
  unsigned int GetNEvaluated() const { return NEvaluated; }
  double GetMeanDLogL() const { return (NEvaluated > 0) ? DLogLSum / static_cast<double>(NEvaluated) : 0.; }
  /// @}

 private:
  /// Totals at the start of the current batch
  struct Totals {
    unsigned int Steps = 0;
    unsigned int Accepted = 0;
    unsigned int OutOfBounds = 0;
    unsigned int NEvaluated = 0;
    double DLogLSum = 0.;
  };
  /// See the constructor
  double Factor;
  /// See the constructor
  double Target;
  /// See the constructor
  double Exponent;
  /// Number of Tune() calls
  int Updates = 0;
  /// Log of the last change Tune() made
  double LastLogStep = 0.;
  /// Running totals
  unsigned int Steps = 0;
  /// See Steps
  unsigned int Accepted_ = 0;
  /// See Steps
  unsigned int OutOfBounds = 0;
  /// See Steps
  unsigned int NEvaluated = 0;
  /// See Steps
  double DLogLSum = 0.;
  /// The totals when the current batch started
  Totals Mark;
};
