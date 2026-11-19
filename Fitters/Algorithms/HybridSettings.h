#pragma once

#include "Manager/YamlHelper.h"
#include "Parameters/ParameterStructs.h"

#include <string>
#include <vector>

/// @brief Everything General.MCMC.Hybrid can set for HybridMCMC, read and checked in one place
///
/// @details Read() applies the defaults (including FittingAlgorithm Auto's) and throws on any value
/// that can't work, so HybridMCMC can take every field at face value.
///
/// @author Henry Wallace
struct HybridSettings {
  /// @brief Read General.MCMC.Hybrid, which may be missing
  /// @param Hybrid The section, or an empty map if the config has none
  /// @param AutoMode Whether FittingAlgorithm is Auto, which changes some defaults
  /// @param WarmupSteps The NUTS warmup length, which Auto scales the tuning interval to
  static HybridSettings Read(const YAML::Node& Hybrid, bool AutoMode, unsigned int WarmupSteps);

  /// @brief Log which scheme runs and what it implies
  void LogScheme() const;

  /// Scheme: true for composition (sweeps then NUTS, every step), false for mixture (a coin per step)
  bool Compose = true;
  /// RandomWalkProb: chance of the random walk on a mixture step
  double RandomWalkProb = 0.5;
  /// Whether RandomWalkProb was given, to warn when composition ignores it
  bool RandomWalkProbGiven = false;

  /// RandomWalkSweeps: sweeps per composition step
  int Sweeps = 1;
  /// Whether RandomWalkSweeps is "auto", chosen during warmup to meet #TimeShare
  bool AutoSweeps = false;
  /// RandomWalkTimeShare: with automatic sweeps, the share of each step's time on the random walk
  double TimeShare = 0.4;

  /// Whether RandomWalkTypes was left out, so the blocks follow the samples' gradient coverage
  bool AutoTypes = false;
  /// RandomWalkTypes: parameter types that go to the random walk
  std::vector<SystType> Types;
  /// RandomWalkNames: names (wildcards allowed) forced into the random walk
  std::vector<std::string> Names;

  /// RandomWalkScale: starting multiple of the optimal \f$2.38/\sqrt{d}\f$ width
  double ScaleFactor = 1.;
  /// RandomWalkTuneScale: whether to tune the width onto #TargetAcceptance during warmup
  bool TuneScale = true;
  /// RandomWalkTargetAcceptance
  double TargetAcceptance = 0.234;
  /// RandomWalkPenaltyExponent: p in \f$\Delta\log L \propto s^p\f$, see RandomWalkTuner
  double PenaltyExponent = 1.;
  /// RandomWalkCheckEvery: steps per tuning batch, 0 for no batches
  int CheckEvery = 1000;

  /// RandomWalkSurrogateStride: every how many events the delayed acceptance surrogate samples, 0 for off
  int SurrogateStride = 0;
  /// FlipMoves: flip moves per step, -1 until the handlers say whether there is anything to flip
  int FlipMoves = -1;
};
