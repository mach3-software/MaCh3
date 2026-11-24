#include "Samples/Oscillation/OscillationGradient.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>

#include "Manager/Manager.h"

//HW Plain C++; the class is only defined (OscillationGradientDevice.cu) in a GPU build, and only made there
#include "Samples/Oscillation/OscillationGradientDevice.h"

_MaCh3_Safe_Include_Start_ //{
#include "Constants/OscillatorConstants.h"
_MaCh3_Safe_Include_End_ //}

// ************************************************
OscillationGradient::OscillationGradient() = default;
OscillationGradient::~OscillationGradient() = default;
// ************************************************

namespace {
  // ************************************************
  /// @brief The energy NuOscillator evaluates a weight at: its bin centre (binned) or its own (unbinned)
  class EvaluationEnergy {
  // ************************************************
   public:
    /// @param Config NuOscillator's config
    explicit EvaluationEnergy(const YAML::Node& Config) {
      if (Config["General"]["CalculationType"].as<std::string>() != "Binned") return;
      Edges = ReadBinEdgesFromFile(Config["Binned"]["FileName"].as<std::string>(), Config["Binned"]["EnergyAxisHistName"].as<std::string>());
      Centres = ReturnBinCentersFromBinEdges(Edges);
    }
    /// @brief The energy a weight for true energy E is evaluated at, -1 outside the binning
    double operator()(const double E) const {
      if (Edges.empty()) return E;
      //HW FLOAT_T is double in a double precision NuOscillator build, which makes these casts redundant there only
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wuseless-cast"
      const auto Upper = std::upper_bound(Edges.begin(), Edges.end(), static_cast<FLOAT_T>(E));
      if (Upper == Edges.begin() || Upper == Edges.end()) return -1.;
      return static_cast<double>(Centres[static_cast<std::size_t>(Upper - Edges.begin() - 1)]);
      #pragma GCC diagnostic pop
    }
   private:
    /// Bin edges, empty for unbinned
    std::vector<FLOAT_T> Edges;
    /// Bin centres
    std::vector<FLOAT_T> Centres;
  };
}

// ************************************************
std::unique_ptr<OscillationGradient> OscillationGradient::Create(const std::string& Name, const std::string& NuOscConfigFile,
                                                                 const std::vector<int>& ParIndex,
                                                                 const std::vector<const M3::float_t*>& ParValue,
                                                                 const std::vector<bool>& Free,
                                                                 const std::vector<EventWeight>& Events) {
// ************************************************
  const YAML::Node Config = M3OpenConfig(NuOscConfigFile);
  if (!SupportsConfig(Config, Name) || !HasAllParameters(ParIndex, ParValue, Free, Name)) return nullptr;

  std::unique_ptr<OscillationGradient> Gradient(new OscillationGradient());
  if (!Gradient->CollectDistinctWeights(Events, EvaluationEnergy(Config), Name)) return nullptr;
  Gradient->SetParameters(ParIndex, ParValue, Free);
  Gradient->GroupBySharedTransition();

  // The derivatives are only offered if the formulae reproduce NuOscillator's own weights, since
  // everything downstream differentiates the formulae rather than the calculation being used
  constexpr double Tolerance = 1e-10;
  const double Difference = Gradient->CompareWithNuOscillator();
  if (!(Difference <= Tolerance)) {
    MACH3LOG_WARN("{}: the oscillation formulae differ from NuOscillator's weights by up to {:.3g}, so no oscillation "
                  "derivatives", Name, Difference);
    return nullptr;
  }
  MACH3LOG_INFO("{}: oscillation derivatives on for {} parameters over {} distinct {} weights; the formulae reproduce "
                "NuOscillator to {:.2g}", Name, Gradient->FreeSlot.size(), Gradient->Weights.size(),
                Config["General"]["CalculationType"].as<std::string>(), Difference);
  Gradient->CheckSharedTransition(Name);
  Gradient->TryDevice(Name);
  Gradient->Evaluate();
  return Gradient;
}

// ************************************************
bool OscillationGradient::SupportsConfig(const YAML::Node& Config, const std::string& Name) {
// ************************************************
  const std::string Calculation = Config["General"]["CalculationType"] ? Config["General"]["CalculationType"].as<std::string>() : "";
  const std::string Implementation = Config["OscProbCalcerSetup"]["ImplementationName"] ?
      Config["OscProbCalcerSetup"]["ImplementationName"].as<std::string>() : "";
  if (Implementation == "CUDAProb3Linear" && (Calculation == "Binned" || Calculation == "Unbinned")) return true;
  MACH3LOG_INFO("{}: no oscillation derivatives, which need CUDAProb3Linear, binned or unbinned (this is {} {})",
                Name, Calculation, Implementation);
  return false;
}

// ************************************************
bool OscillationGradient::HasAllParameters(const std::vector<int>& ParIndex, const std::vector<const M3::float_t*>& ParValue,
                                           const std::vector<bool>& Free, const std::string& Name) {
// ************************************************
  if (ParIndex.size() != M3::Osc::kNPars || ParValue.size() != M3::Osc::kNPars || Free.size() != M3::Osc::kNPars) {
    MACH3LOG_INFO("{}: no oscillation derivatives, which need all eight CUDAProb3Linear parameters", Name);
    return false;
  }
  for (int k = 0; k < M3::Osc::kNPars; ++k) {
    if (ParIndex[static_cast<std::size_t>(k)] < 0 || ParValue[static_cast<std::size_t>(k)] == nullptr) {
      MACH3LOG_INFO("{}: no oscillation derivatives, oscillation parameter {} of CUDAProb3Linear's eight is missing", Name, k);
      return false;
    }
  }
  return true;
}

// ************************************************
template <typename EnergyMap>
bool OscillationGradient::CollectDistinctWeights(const std::vector<EventWeight>& Events, const EnergyMap& EvaluatedAt,
                                                 const std::string& Name) {
// ************************************************
  // One entry per distinct weight; every event sharing it has to agree on what it is
  std::map<const M3::float_t*, std::size_t, std::less<const M3::float_t*>> Index;
  std::vector<EventWeight> Distinct;
  for (const auto& Ev : Events) {
    EventWeight W = Ev;
    W.Energy = EvaluatedAt(Ev.Energy);
    const auto Found = Index.emplace(Ev.Weight, Distinct.size());
    if (Found.second) { Distinct.push_back(W); continue; }
    const EventWeight& Seen = Distinct[Found.first->second];
    if (Seen.Energy != W.Energy || Seen.Anti != W.Anti || Seen.From != W.From || Seen.To != W.To) {
      MACH3LOG_WARN("{}: events sharing an oscillation weight disagree on its energy or channel, so no oscillation "
                    "derivatives", Name);
      return false;
    }
  }
  // In pointer order
  for (const auto& [Pointer, i] : Index) {
    const EventWeight& W = Distinct[i];
    Weights.push_back(Pointer);
    Energy.push_back(W.Energy);
    Anti.push_back(W.Anti ? 1 : 0);
    From.push_back(W.From);
    To.push_back(W.To);
  }
  return true;
}

// ************************************************
void OscillationGradient::SetParameters(const std::vector<int>& ParIndex, const std::vector<const M3::float_t*>& ParValue_,
                                        const std::vector<bool>& Free) {
// ************************************************
  ParValue = ParValue_;
  for (int k = 0; k < M3::Osc::kNPars; ++k) {
    if (!Free[static_cast<std::size_t>(k)]) continue;
    FreeSlot.push_back(k);
    FreeParam.push_back(ParIndex[static_cast<std::size_t>(k)]);
  }
  Derivatives.assign(Weights.size() * FreeSlot.size(), 0.);
}

// ************************************************
void OscillationGradient::GroupBySharedTransition() {
// ************************************************
  // Binned oscillation has several channels per energy bin, and they all come out of one transition
  // matrix, so the weights are grouped by energy and neutrino type
  std::vector<int> Order(Weights.size());
  for (std::size_t i = 0; i < Order.size(); ++i) Order[i] = static_cast<int>(i);
  std::stable_sort(Order.begin(), Order.end(), [&](const int a, const int b) {
    const std::size_t i = static_cast<std::size_t>(a), j = static_cast<std::size_t>(b);
    if (Anti[i] != Anti[j]) return Anti[i] < Anti[j];
    return Energy[i] < Energy[j];
  });
  for (std::size_t k = 0; k < Order.size(); ++k) {
    const std::size_t w = static_cast<std::size_t>(Order[k]);
    const std::size_t Previous = (k == 0) ? 0 : static_cast<std::size_t>(Order[k - 1]);
    const bool New = (k == 0) || Energy[w] != Energy[Previous] || Anti[w] != Anti[Previous];
    if (New) SharedStart.push_back(static_cast<int>(k));
    SharedWeight.push_back(Order[k]);
  }
  SharedStart.push_back(static_cast<int>(Order.size()));
}

// ************************************************
void OscillationGradient::CheckSharedTransition(const std::string& Name) {
// ************************************************
  double Values[M3::Osc::kNPars];
  CurrentValues(Values);
  M3::Osc::OscDual Dual[M3::Osc::kNPars];
  for (int k = 0; k < M3::Osc::kNPars; ++k) Dual[k] = M3::Osc::OscDual::Variable(Values[k], k);
  M3::Osc::Constants<M3::Osc::OscDual> C;
  if (!M3::Osc::MakeConstants(Dual, C)) return;
  const double Shared = CompareSharedWithSingle(C);
  if (Shared != 0.) MACH3LOG_WARN("{}: oscillation derivatives from shared transition matrices differ from one weight at a "
                                  "time by up to {:.3g}", Name, Shared);
  MACH3LOG_INFO("{}: {} oscillation weights share {} transition matrices (energy and neutrino type)", Name,
                Weights.size(), SharedStart.size() - 1);
}

// ************************************************
void OscillationGradient::CurrentValues(double* Values) const {
// ************************************************
  for (int k = 0; k < M3::Osc::kNPars; ++k) Values[k] = static_cast<double>(*ParValue[static_cast<std::size_t>(k)]);
}

// ************************************************
double OscillationGradient::CompareWithNuOscillator() const {
// ************************************************
  double Values[M3::Osc::kNPars];
  CurrentValues(Values);
  M3::Osc::Constants<double> C;
  if (!M3::Osc::MakeConstants(Values, C)) return 1e30;
  double Worst = 0.;
  for (std::size_t i = 0; i < Weights.size(); ++i) {
    const double Formula = M3::Osc::Probability(C, Energy[i], Anti[i] != 0, From[i], To[i]);
    const double Expected = static_cast<double>(*Weights[i]);
    Worst = std::max(Worst, std::fabs(Formula - Expected) / std::max(std::fabs(Expected), 1e-12));
  }
  return Worst;
}

// ************************************************
void OscillationGradient::Evaluate() {
// ************************************************
  std::vector<double> Values(M3::Osc::kNPars);
  CurrentValues(Values.data());
  if (Values == LastValues) return;
  LastValues = Values;

  M3::Osc::OscDual Dual[M3::Osc::kNPars];
  for (int k = 0; k < M3::Osc::kNPars; ++k) Dual[k] = M3::Osc::OscDual::Variable(Values[static_cast<std::size_t>(k)], k);
  M3::Osc::Constants<M3::Osc::OscDual> C;
  if (!M3::Osc::MakeConstants(Dual, C)) {
    MACH3LOG_ERROR("Oscillation parameters NuOscillator can't take: a mixing parameter is negative");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  #ifdef MaCh3_CUDA
  //HW M3::float_t is double in most builds, which makes the casts below redundant there only
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuseless-cast"
  if (Device) {
    Device->Evaluate(C, FreeSlot, DeviceOut.data());
    const std::size_t N = Derivatives.size();
    for (std::size_t i = 0; i < N; ++i) Derivatives[i] = static_cast<M3::float_t>(DeviceOut[i]);
    return;
  }
  #pragma GCC diagnostic pop
  #endif
  EvaluateOnHost(C);
}

// ************************************************
void OscillationGradient::EvaluateOnHost(const M3::Osc::Constants<M3::Osc::OscDual>& C) {
// ************************************************
  // One transition matrix per energy and neutrino type, then each of its weights' channels. Dynamic
  // scheduling, since cores differ in speed (Apple's efficiency cores); each group writes only its
  // own weights, so the result doesn't depend on which thread took it
  const int NGroups = static_cast<int>(SharedStart.size()) - 1;
  const std::size_t NFree = FreeSlot.size();
  #ifdef MULTITHREAD
  #pragma omp parallel for schedule(dynamic, 16)
  #endif
  for (int g = 0; g < NGroups; ++g) {
    const std::size_t First = static_cast<std::size_t>(SharedWeight[static_cast<std::size_t>(SharedStart[static_cast<std::size_t>(g)])]);
    M3::Osc::Complex<M3::Osc::OscDual> A[3][3];
    M3::Osc::AllAmplitudes(C, Energy[First], Anti[First] != 0, A);
    for (int k = SharedStart[static_cast<std::size_t>(g)]; k < SharedStart[static_cast<std::size_t>(g) + 1]; ++k) {
      const std::size_t w = static_cast<std::size_t>(SharedWeight[static_cast<std::size_t>(k)]);
      const M3::Osc::OscDual P = M3::Osc::ProbabilityFromAmplitude(A[To[w]][From[w]]);
      for (std::size_t j = 0; j < NFree; ++j) {
        Derivatives[w * NFree + j] = static_cast<M3::float_t>(P.d[FreeSlot[j]]);
      }
    }
  }
}

// ************************************************
double OscillationGradient::CompareSharedWithSingle(const M3::Osc::Constants<M3::Osc::OscDual>& C) {
// ************************************************
  //HW M3::float_t is double in most builds, which makes the cast below redundant there only
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuseless-cast"
  EvaluateOnHost(C);
  double Worst = 0.;
  const std::size_t NFree = FreeSlot.size();
  for (std::size_t w = 0; w < Weights.size(); ++w) {
    const M3::Osc::OscDual P = M3::Osc::Probability(C, Energy[w], Anti[w] != 0, From[w], To[w]);
    for (std::size_t j = 0; j < NFree; ++j) {
      Worst = std::max(Worst, std::fabs(static_cast<double>(Derivatives[w * NFree + j]) - P.d[FreeSlot[j]]));
    }
  }
  #pragma GCC diagnostic pop
  return Worst;
}

// ************************************************
void OscillationGradient::TryDevice(const std::string& Name) {
// ************************************************
  #ifdef MaCh3_CUDA
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuseless-cast"
  if (!OscillationGradientDevice::Available() || Weights.empty() || FreeSlot.empty()) return;
  double Values[M3::Osc::kNPars];
  CurrentValues(Values);
  M3::Osc::OscDual Dual[M3::Osc::kNPars];
  for (int k = 0; k < M3::Osc::kNPars; ++k) Dual[k] = M3::Osc::OscDual::Variable(Values[k], k);
  M3::Osc::Constants<M3::Osc::OscDual> C;
  if (!M3::Osc::MakeConstants(Dual, C)) return;
  try {
    auto Candidate = std::make_shared<OscillationGradientDevice>(Energy, Anti, From, To);
    std::vector<double> Out(Derivatives.size());
    Candidate->Evaluate(C, FreeSlot, Out.data());
    EvaluateOnHost(C);
    //HW Relative to each parameter's largest derivative: the GPU's sin, cos and sqrt may round
    //differently in the last place, so this is a rounding-level comparison, not bit for bit
    const std::size_t NFree = FreeSlot.size();
    std::vector<double> Scale(NFree, 0.);
    for (std::size_t i = 0; i < Out.size(); ++i) Scale[i % NFree] = std::max(Scale[i % NFree], std::fabs(static_cast<double>(Derivatives[i])));
    double Worst = 0.;
    for (std::size_t i = 0; i < Out.size(); ++i) {
      const double Ref = static_cast<double>(Derivatives[i]);
      Worst = std::max(Worst, std::fabs(Out[i] - Ref) / std::max(Scale[i % NFree], 1e-300));
    }
    constexpr double Tolerance = 1e-12;
    if (!(Worst <= Tolerance)) {
      MACH3LOG_WARN("{}: oscillation derivatives on the GPU differ from the CPU by up to {:.3g}, so they stay on the CPU", Name, Worst);
      return;
    }
    Device = std::move(Candidate);
    DeviceOut.assign(Derivatives.size(), 0.);
    LastValues.clear();
    MACH3LOG_INFO("{}: oscillation derivatives on the GPU, matching the CPU to {:.2g}", Name, Worst);
  } catch (const std::exception& Error) {
    MACH3LOG_WARN("{}: oscillation derivatives stay on the CPU, the GPU failed: {}", Name, Error.what());
  }
  #pragma GCC diagnostic pop
  #else
  (void)Name;
  #endif
}

// ************************************************
WeightDerivativeSource OscillationGradient::Source() const {
// ************************************************
  WeightDerivativeSource S;
  S.Listed = Weights;
  S.Derivatives = Derivatives.data();
  S.MultiParam = FreeParam;
  return S;
}
