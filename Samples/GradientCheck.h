#pragma once

#include <string>
#include <vector>

class ParameterHandlerBase;
class SampleHandlerInterface;

/// @brief Settings for CheckLikelihoodGradient()
struct GradientCheckSettings {
  /// Parameter points to check at
  int NPoints = 3;
  /// Parameters checked against finite differences at each point, chosen at random among those
  /// with a gradient. Each costs two reweights. 0 checks all of them
  int ParsPerPoint = 30;
  /// Points are thrown around the prefit values with this many prior widths, clipped to the bounds
  double Width = 0.5;
  /// Finite difference step, as a fraction of each parameter's prior width
  double RelativeStep = 1e-4;
  /// Largest acceptable disagreement, relative to the larger of the two derivatives (or to a floor
  /// of 1e-4 of the largest gradient entry at that point, for derivatives near zero)
  double Tolerance = 1e-5;
  /// Seed for the points and the choice of parameters
  unsigned int Seed = 2026;
};

/// @brief What CheckLikelihoodGradient() found
struct GradientCheckResult {
  /// Largest disagreement over every parameter checked, in the relative terms of the settings
  double Worst = 0.;
  /// Name of the parameter it was found for
  std::string WorstParameter;
  /// Number of parameter and point combinations compared
  int NChecked = 0;
  /// How many of those were over the tolerance
  int NFailed = 0;
  /// Whether everything was within the tolerance
  bool Passed() const { return NFailed == 0 && NChecked > 0; }
};

/// @brief Check a sample handler's likelihood gradient against central finite differences
///
/// @details Switches gradients on, moves the parameters through Settings.NPoints points and at each
/// compares GetGradLikelihood() with \f$(L(x+h) - L(x-h))/2h\f$ of GetLikelihood() for a random
/// selection of the parameters that have a gradient (GetGradientCoverage()). Parameters within a
/// step of a bound are skipped. The prior gradient, ParameterHandlerBase::CalcGradient(), is
/// checked the same way at every point, for every parameter, since it costs no reweights.
///
/// Every disagreement over the tolerance is logged. The parameters are put back and the handler
/// reweighted before returning.
///
/// @param Handler The sample handler, whose parameters are Pars
/// @param Pars Its parameter handler
/// @param Settings See GradientCheckSettings
GradientCheckResult CheckLikelihoodGradient(SampleHandlerInterface& Handler, ParameterHandlerBase& Pars,
                                            const GradientCheckSettings& Settings = {});
