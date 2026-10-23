#pragma once

#include <algorithm>
#include <cmath>

/// @brief Dual averaging step size adaptation for Hamiltonian samplers
/// @details Nudges the integrator step size during a warmup phase until the mean Metropolis
/// acceptance over a trajectory sits at #Target, then Freeze() swaps in the running average, which
/// is far steadier than the last raw value. Deliberately free of MaCh3 dependencies so it can be
/// tested on its own.
/// @author Henry Wallace
struct DualAveragingStepSize {
  /// @brief Start adapting from scratch around an initial step size
  /// @param InitialStepSize The step size to shrink the search towards
  void Restart(const double InitialStepSize) {
    // Biased towards larger steps, as the paper suggests, so the search approaches from above
    Mu = std::log(10. * InitialStepSize);
    HBar = 0.;
    LogStepSizeBar = std::log(InitialStepSize);
    StepSize = InitialStepSize;
    Counter = 0;
  }

  /// @brief One adaptation iteration
  /// @param MeanAccept Mean \f$\min(1, e^{H_0 - H})\f$ over the trajectory just taken
  void Update(const double MeanAccept) {
    // A trajectory with nothing to average over tells us nothing
    if (!std::isfinite(MeanAccept)) return;

    const double Accept = std::min(1., std::max(0., MeanAccept));
    const double Count = static_cast<double>(++Counter);

    // Running average of how far off target we are, then a step size that over-corrects early on
    // and settles as sqrt(m) grows
    HBar = (1. - 1. / (Count + T0)) * HBar + (Target - Accept) / (Count + T0);
    const double LogStepSize = Mu - std::sqrt(Count) / Gamma * HBar;

    // The averaging weight starts at 1, so the very first update discards the initial guess
    const double Weight = std::pow(Count, -Kappa);
    LogStepSizeBar = Weight * LogStepSize + (1. - Weight) * LogStepSizeBar;

    StepSize = std::exp(LogStepSize);
    if (StepSize < MinStepSize) {
      StepSize = MinStepSize;
      ++NClamped;
    } else if (StepSize > MaxStepSize) {
      StepSize = MaxStepSize;
      ++NClamped;
    }
  }

  /// @brief Stop adapting and keep the averaged step size rather than the last raw one
  void Freeze() {
    StepSize = std::min(MaxStepSize, std::max(MinStepSize, std::exp(LogStepSizeBar)));
  }

  /// Acceptance statistic we are steering towards
  double Target = 0.8;
  /// How hard the correction pulls, the paper's gamma
  double Gamma = 0.05;
  /// How fast the running average forgets, the paper's kappa
  double Kappa = 0.75;
  /// Damps the early iterations, the paper's t0
  double T0 = 10.;

  /// Smallest step size we will hand back
  double MinStepSize = 1e-10;
  /// Largest step size we will hand back
  double MaxStepSize = 1e3;

  /// Point the search shrinks towards, \f$\log(10\,\epsilon_0)\f$
  double Mu = 0.;
  /// Running average of Target minus the acceptance
  double HBar = 0.;
  /// Running average of the log step size, what Freeze() keeps
  double LogStepSizeBar = 0.;
  /// The step size to use right now
  double StepSize = 0.;
  /// Iterations so far
  long Counter = 0;
  /// How often the step size hit one of the clamps, a sign the geometry is pathological
  long NClamped = 0;
};
