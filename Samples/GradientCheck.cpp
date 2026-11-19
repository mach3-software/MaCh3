#include "Samples/GradientCheck.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <numeric>
#include <random>

#include "Parameters/ParameterHandlerGeneric.h"
#include "Samples/SampleHandlerInterface.h"

namespace {

/// @brief Throw every parameter with a gradient around its prefit value, clipped to its bounds
void ThrowPoint(ParameterHandlerBase& Pars, const std::vector<bool>& Covered, const double Width, std::mt19937& Rng) {
  std::normal_distribution<double> Gaus(0., 1.);
  for (int i = 0; i < Pars.GetNumParams(); ++i) {
    if (!Covered[static_cast<std::size_t>(i)] || Pars.IsParameterFixed(i)) continue;
    const double Value = Pars.GetParPreFit(i) + Width * Pars.GetDiagonalError(i) * Gaus(Rng);
    Pars.SetParProp(i, std::clamp(Value, Pars.GetLowerBound(i), Pars.GetUpperBound(i)));
  }
}

/// @brief Central finite difference of a function of one parameter, which is restored afterwards
template <typename Function>
double FiniteDifference(ParameterHandlerBase& Pars, const int i, const double Step, Function&& Evaluate) {
  const double Value = Pars.GetParProp(i);
  Pars.SetParProp(i, Value + Step);
  const double Up = Evaluate();
  Pars.SetParProp(i, Value - Step);
  const double Down = Evaluate();
  Pars.SetParProp(i, Value);
  return (Up - Down) / (2. * Step);
}

/// @brief Disagreement between an analytic and a finite difference derivative, relative to the
/// larger of the two or to Floor
double Disagreement(const double Analytic, const double Numeric, const double Floor) {
  return std::fabs(Analytic - Numeric) / std::max({std::fabs(Analytic), std::fabs(Numeric), Floor});
}

/// @brief Fold one comparison into the result, logging it if it is over the tolerance
void Tally(GradientCheckResult& Result, const double Difference, const double Tolerance, const std::string& What,
           const double Analytic, const double Numeric) {
  ++Result.NChecked;
  if (Difference > Result.Worst) {
    Result.Worst = Difference;
    Result.WorstParameter = What;
  }
  if (Difference > Tolerance) {
    ++Result.NFailed;
    MACH3LOG_ERROR("Gradient check: {} analytic {:.10g}, finite difference {:.10g}, relative difference {:.3g}",
                   What, Analytic, Numeric, Difference);
  }
}

/// @brief The usable step for parameter i, or 0 if it sits too close to a bound to step both ways
double StepFor(const ParameterHandlerBase& Pars, const int i, const double RelativeStep) {
  const double Step = RelativeStep * Pars.GetDiagonalError(i);
  const double Value = Pars.GetParProp(i);
  if (!(Step > 0.) || Value - Step < Pars.GetLowerBound(i) || Value + Step > Pars.GetUpperBound(i)) return 0.;
  return Step;
}

/// @brief Check the prior gradient for every parameter, at the current point
void CheckPrior(ParameterHandlerBase& Pars, const GradientCheckSettings& Settings, GradientCheckResult& Result) {
  const std::vector<double> Analytic = Pars.CalcGradient();
  double Largest = 0.;
  for (const double g : Analytic) Largest = std::max(Largest, std::fabs(g));
  for (int i = 0; i < Pars.GetNumParams(); ++i) {
    const double Step = StepFor(Pars, i, Settings.RelativeStep);
    if (Step == 0.) continue;
    const double Numeric = FiniteDifference(Pars, i, Step, [&Pars]() { return Pars.CalcLikelihood(); });
    const double g = Analytic[static_cast<std::size_t>(i)];
    Tally(Result, Disagreement(g, Numeric, 1e-4 * Largest), Settings.Tolerance, "prior of " + Pars.GetParFancyName(i), g, Numeric);
  }
}

/// @brief A label for the kind of parameter, so the choice below can take some of each
std::string KindOf(const ParameterHandlerBase& Pars, const int i) {
  const auto* Generic = dynamic_cast<const ParameterHandlerGeneric*>(&Pars);
  return Generic ? Generic->GetParamTypeString(i) : "Parameter";
}

/// @brief Parameters to compare at one point: all of those with a gradient, or HowMany of them
/// @details Taken from each kind of parameter in turn, so a few hundred of one kind (energy
/// normalisations, say) can't crowd the splines out of a small selection
std::vector<int> ChooseParameters(const ParameterHandlerBase& Pars, const std::vector<bool>& Covered,
                                  const int HowMany, std::mt19937& Rng) {
  std::map<std::string, std::vector<int>> ByKind;
  for (int i = 0; i < Pars.GetNumParams(); ++i) {
    if (Covered[static_cast<std::size_t>(i)] && !Pars.IsParameterFixed(i)) ByKind[KindOf(Pars, i)].push_back(i);
  }
  for (auto& [Kind, List] : ByKind) std::shuffle(List.begin(), List.end(), Rng);
  std::vector<int> Chosen;
  for (std::size_t Round = 0; HowMany <= 0 || static_cast<int>(Chosen.size()) < HowMany; ++Round) {
    bool Any = false;
    for (const auto& [Kind, List] : ByKind) {
      if (Round >= List.size()) continue;
      Any = true;
      Chosen.push_back(List[Round]);
      if (HowMany > 0 && static_cast<int>(Chosen.size()) == HowMany) break;
    }
    if (!Any) break;
  }
  std::sort(Chosen.begin(), Chosen.end());
  return Chosen;
}

/// @brief Log the parameters without a gradient, grouped by kind, so a missing one stands out
void ReportUncovered(const ParameterHandlerBase& Pars, const std::vector<bool>& Covered) {
  std::map<std::string, std::vector<std::string>> ByKind;
  for (int i = 0; i < Pars.GetNumParams(); ++i) {
    if (!Covered[static_cast<std::size_t>(i)]) ByKind[KindOf(Pars, i)].push_back(Pars.GetParFancyName(i));
  }
  for (const auto& [Kind, Names] : ByKind) {
    std::string List;
    for (const auto& Name : Names) List += (List.empty() ? "" : ", ") + Name;
    MACH3LOG_INFO("Gradient check: {} {} parameters have no gradient: {}", Names.size(), Kind, List);
  }
}

} // namespace

// ************************************************
GradientCheckResult CheckLikelihoodGradient(SampleHandlerInterface& Handler, ParameterHandlerBase& Pars,
                                            const GradientCheckSettings& Settings) {
// ************************************************
  Handler.SetComputeGradients(true);
  const std::vector<bool> Covered = Handler.GetGradientCoverage();
  if (static_cast<int>(Covered.size()) != Pars.GetNumParams()) {
    MACH3LOG_ERROR("{} reports gradient coverage for {} parameters, but {} has {}", Handler.GetName(),
                   Covered.size(), Pars.GetName(), Pars.GetNumParams());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  ReportUncovered(Pars, Covered);
  std::vector<double> Original(static_cast<std::size_t>(Pars.GetNumParams()));
  for (int i = 0; i < Pars.GetNumParams(); ++i) Original[static_cast<std::size_t>(i)] = Pars.GetParProp(i);

  std::mt19937 Rng(Settings.Seed);
  GradientCheckResult Result;
  const auto SampleLikelihood = [&Handler]() { Handler.Reweight(); return Handler.GetLikelihood(); };
  for (int Point = 0; Point < Settings.NPoints; ++Point) {
    ThrowPoint(Pars, Covered, Settings.Width, Rng);
    CheckPrior(Pars, Settings, Result);

    Handler.Reweight();
    const std::vector<double> Analytic = Handler.GetGradLikelihood();
    double Largest = 0.;
    for (const double g : Analytic) Largest = std::max(Largest, std::fabs(g));

    const std::vector<int> Chosen = ChooseParameters(Pars, Covered, Settings.ParsPerPoint, Rng);
    MACH3LOG_INFO("Gradient check, point {} of {}: {} parameters against finite differences, largest gradient entry {:.4g}",
                  Point + 1, Settings.NPoints, Chosen.size(), Largest);
    for (const int i : Chosen) {
      const double Step = StepFor(Pars, i, Settings.RelativeStep);
      if (Step == 0.) continue;
      const double Numeric = FiniteDifference(Pars, i, Step, SampleLikelihood);
      const double g = Analytic[static_cast<std::size_t>(i)];
      Tally(Result, Disagreement(g, Numeric, 1e-4 * Largest), Settings.Tolerance, Pars.GetParFancyName(i), g, Numeric);
    }
  }

  for (int i = 0; i < Pars.GetNumParams(); ++i) Pars.SetParProp(i, Original[static_cast<std::size_t>(i)]);
  Handler.Reweight();
  MACH3LOG_INFO("Gradient check: {} comparisons, {} over the tolerance of {:.1g}, worst {:.3g} ({})",
                Result.NChecked, Result.NFailed, Settings.Tolerance, Result.Worst, Result.WorstParameter);
  return Result;
}
