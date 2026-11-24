// SampleHandlerBase: choosing and checking the reweight backend, the oscillation and spline weight
// cache, and the likelihood gradient. Kept apart from SampleHandlerBase.cpp, which deals with
// loading and describing samples.
#include "Samples/SampleHandlerBase.h"
#include "Samples/Reweight/LinearShiftBuilder.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace {

/// @brief The largest difference between two arrays, relative to each reference entry, or to
/// Floor times the largest reference entry for entries near zero
double LargestRelativeDifference(const std::vector<double>& Reference, const std::vector<double>& Other,
                                 const double Floor, int& Where) {
  double Largest = 0.;
  for (const double r : Reference) Largest = std::max(Largest, std::fabs(r));
  const double Scale = std::max(Floor * Largest, 1e-300);
  double Worst = 0.;
  Where = -1;
  for (std::size_t i = 0; i < Reference.size(); ++i) {
    if (Other[i] == Reference[i]) continue;
    const double Difference = std::fabs(Other[i] - Reference[i]) / std::max(std::fabs(Reference[i]), Scale);
    if (Difference > Worst) { Worst = Difference; Where = static_cast<int>(i); }
  }
  return Worst;
}

/// @brief Everything one reweight produced, for comparing two of them
struct ReweightSnapshot {
  std::vector<double> MC;
  std::vector<double> W2;
  std::vector<double> Gradient;
  double Likelihood = 0.;
};

} // namespace

// ************************************************
void SampleHandlerBase::SetupReweightBackend() {
// ************************************************
  const ReweightOptions Options = ReweightOptions::Read(SampleManager->raw(), GetName());
  if (Options.Backend == "legacy") {
    MACH3LOG_INFO("{}: reweighting with the plain event loop, as configured", GetName());
    return;
  }
  if (!AdditionalWeightsAreConstant()) {
    if (Options.Backend != "auto") {
      MACH3LOG_ERROR("{}: Reweight.Backend is {}, but the handler doesn't declare its additional weights "
                     "constant (AdditionalWeightsAreConstant), which that backend relies on", GetName(), Options.Backend);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    MACH3LOG_INFO("{}: reweighting with the plain event loop, since the handler doesn't declare its "
                  "additional weights constant (AdditionalWeightsAreConstant)", GetName());
    return;
  }

  // Build the backend's description of the events, leaving out any that can never be selected
  const std::vector<char> NeverSelected = FindNeverSelected();
  ReweightModel Model = BuildReweightModel(MCEvents, ConstWeightStart, ConstWeightEnd, NeverSelected, Binning->GetNBins());
  // How the backend shifts, cuts and bins events: the experiment's own code, or the faster pointer
  // and table lookups if the config asks for them. Nothing is checked here; the start-up self-check
  // at the end of this function is what compares the result with the plain loop
  BackendAccess = std::make_unique<BackendEventAccess>(*this);
  ReweightEventAccess* Access = BackendAccess.get();
  if (Options.FastAccess) {
    SetupFastEventAccess(NeverSelected);
    if (Options.LinearShifts) SetupLinearShifts(NeverSelected);
    Access = FastAccess.get();
  } else if (Options.LinearShifts) {
    MACH3LOG_WARN("{}: Reweight.LinearShifts needs FastAccess, so is ignored", GetName());
  }
  Backend = MakeReweightBackend(Options.Backend, std::move(Model), *Access);
  MACH3LOG_INFO("{}: reweighting with the {} backend. {}", GetName(), Backend->Name(), Backend->Summary());
  // Gradients switched on before the backend existed have to be re-attached to it
  if (Gradient) SetComputeGradients(true);
  ReportReweightMemory();
  SelfCheckReweightBackend(Options);
}

// ************************************************
void SampleHandlerBase::SelfCheckReweightBackend(const ReweightOptions& Options) {
// ************************************************
  if (Options.SelfCheckPoints <= 0) return;
  const double Tolerance = Backend->Precision();
  const double Worst = CheckReweightBackend(Options.SelfCheckPoints, Tolerance);
  if (Worst <= Tolerance) return;
  if (Options.Backend != "auto") {
    MACH3LOG_ERROR("{}: the {} backend differs from the plain event loop by up to {:.3g}", GetName(), Backend->Name(), Worst);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  MACH3LOG_WARN("{}: the {} backend differs from the plain event loop by up to {:.3g}, so falling back to "
                "the plain loop. Please report this", GetName(), Backend->Name(), Worst);
  Backend.reset();
  FastAccess.reset();
  BackendAccess.reset();
}

// ************************************************
void SampleHandlerBase::SetupFastEventAccess(const std::vector<char>& NeverSelected) {
// ************************************************
  FastAccess = std::make_unique<FastEventAccess>(*BackendAccess);
  const std::size_t NEvents = GetNEvents();
  std::vector<int> Sample(NEvents);
  for (std::size_t e = 0; e < NEvents; ++e) Sample[e] = MCEvents[e].NominalSample;

  // Cuts on shifted variables, through pointers. A live event already passes every other cut
  std::vector<std::vector<KinematicCut>> ShiftedCuts(static_cast<std::size_t>(GetNSamples()));
  for (int s = 0; s < GetNSamples(); ++s) {
    for (const auto& Cut : StoredSelection[static_cast<std::size_t>(s)]) {
      if (IsKinematicParameterShifted(Cut.ParamToCutOnIt)) ShiftedCuts[static_cast<std::size_t>(s)].push_back(Cut);
    }
  }
  std::vector<std::size_t> Start(NEvents, 0);
  std::vector<const double*> Value;
  bool CutsUsable = true;
  try {
    for (std::size_t e = 0; e < NEvents && CutsUsable; ++e) {
      Start[e] = Value.size();
      if (NeverSelected[e]) continue;
      for (const auto& Cut : ShiftedCuts[static_cast<std::size_t>(Sample[e])]) {
        const double* Pointer = GetPointerToKinematicParameter(Cut.ParamToCutOnIt, static_cast<int>(e));
        const double Expected = ReturnKinematicParameter(Cut.ParamToCutOnIt, static_cast<int>(e));
        //HW A variable that only exists as a computed value can't be read through a pointer
        if (Pointer == nullptr || !(*Pointer == Expected || (std::isnan(*Pointer) && std::isnan(Expected)))) {
          CutsUsable = false;
          break;
        }
        Value.push_back(Pointer);
      }
    }
  } catch (const std::exception&) {
    CutsUsable = false;
  }
  if (CutsUsable) FastAccess->UseCutPointers(Sample, std::move(ShiftedCuts), std::move(Start), std::move(Value));

  // Uniform one dimensional binning, through each sample binning's own FindBin
  bool BinsUsable = true;
  std::vector<const SampleBinningInfo*> SampleBinning(static_cast<std::size_t>(GetNSamples()));
  for (int s = 0; s < GetNSamples(); ++s) {
    const SampleBinningInfo& Info = Binning->GetSampleBinningInfo(s);
    BinsUsable = BinsUsable && Info.Uniform && Info.BinEdges.size() == 1;
    SampleBinning[static_cast<std::size_t>(s)] = &Info;
  }
  std::vector<const double*> Var(NEvents, nullptr);
  std::vector<int> Nominal(NEvents, M3::UnderOverFlowBin);
  for (std::size_t e = 0; e < NEvents && BinsUsable; ++e) {
    if (MCEvents[e].KinVar.size() != 1 || MCEvents[e].NomBin.size() != 1) { BinsUsable = false; break; }
    Var[e] = MCEvents[e].KinVar[0];
    Nominal[e] = MCEvents[e].NomBin[0];
  }
  // Check it against FindGlobalBin for every live event before using it
  for (std::size_t e = 0; e < NEvents && BinsUsable; ++e) {
    if (NeverSelected[e]) continue;
    const int Flat = FastEventAccess::FlatBin(*SampleBinning[static_cast<std::size_t>(Sample[e])], *Var[e], Nominal[e]);
    if (Flat != BackendAccess->GlobalBinOf(static_cast<int>(e))) {
      MACH3LOG_WARN("{}: the flat bin lookup disagrees with FindGlobalBin for event {}, so it isn't used", GetName(), e);
      BinsUsable = false;
    }
  }
  if (BinsUsable) FastAccess->UseFlatBinning(std::move(Sample), std::move(SampleBinning), std::move(Var), std::move(Nominal));
  MACH3LOG_INFO("{}: shifted events are cut {} and binned {}", GetName(),
                FastAccess->UsesCutPointers() ? "through pointers" : "through ReturnKinematicParameter",
                FastAccess->UsesFlatBinning() ? "with the flat lookup" : "with FindGlobalBin");
}

// ************************************************
void SampleHandlerBase::SetupLinearShifts(const std::vector<char>& NeverSelected) {
// ************************************************
  if (!FastAccess->UsesCutPointers() || !FastAccess->UsesFlatBinning()) {
    MACH3LOG_WARN("{}: linear shifts need the cut pointers and the flat bin lookup, which aren't in use, so "
                  "they are off", GetName());
    return;
  }
  functional.update_vals();
  PrepFunctionalParameters();
  RefreshShiftTable();
  if (!functional.table_ready) {
    MACH3LOG_WARN("{}: linear shifts need every shift registered as a plain function of one parameter, so "
                  "they are off", GetName());
    return;
  }
  const std::size_t NShifts = functional.raw_table.size();
  if (NShifts > 65535 || GetNEvents() > 4000000000u) {
    MACH3LOG_WARN("{}: too many shifts or events for the linear shift table, so linear shifts are off", GetName());
    return;
  }
  // Each shift is stepped by its parameter's prior width, which is the scale it is linear over.
  // A shift that isn't a function of exactly one parameter, or whose parameter has no prior width,
  // keeps a step of 1
  std::vector<double> Width(NShifts, 1.);
  for (std::size_t s = 0; s < NShifts; ++s) {
    const auto& Indices = functional.shifts[s].par_indices;
    if (Indices.size() == 1 && ParHandler->GetDiagonalError(Indices[0]) > 0.) Width[s] = ParHandler->GetDiagonalError(Indices[0]);
  }
  LinearShiftTable Linear = LinearShiftBuilder(functional, *BackendAccess, *FastAccess, NeverSelected, std::move(Width)).Build();

  std::size_t NLinear = 0, NLive = 0;
  for (std::size_t e = 0; e < GetNEvents(); ++e) { NLinear += Linear.Linear[e] ? 1 : 0; NLive += NeverSelected[e] ? 0 : 1; }
  MACH3LOG_INFO("{}: linear shifts hold for {} of {} live events ({:.1f} MB of coefficients); the rest go "
                "through the experiment's shifts", GetName(), NLinear, NLive, static_cast<double>(Linear.Bytes()) / 1e6);
  if (NLinear == 0) return;
  FastAccess->UseLinearShifts(std::move(Linear));
}

// ************************************************
void SampleHandlerBase::SyncShiftedKinematics() {
// ************************************************
  if (!KinematicsStale) return;
  const int NEvents = static_cast<int>(GetNEvents());
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int iEvent = 0; iEvent < NEvents; ++iEvent) ApplyShifts(iEvent);
  KinematicsStale = false;
}

// ************************************************
void SampleHandlerBase::ReportReweightMemory() const {
// ************************************************
  const auto MB = [](const std::size_t Bytes) { return static_cast<double>(Bytes) / 1e6; };
  MACH3LOG_INFO("{}: reweight memory: backend {:.1f} MB, fast lookups {:.1f} MB, shift lists {:.1f} MB ({} distinct), "
                "event gradient record {:.1f} MB", GetName(), MB(Backend ? Backend->Bytes() : 0),
                MB(FastAccess ? FastAccess->Bytes() : 0), MB(functional.event_shifts.Bytes()),
                functional.event_shifts.NLists(), MB(Gradient ? GetNEvents() * 12 : 0));
}

// ************************************************
std::vector<char> SampleHandlerBase::FindNeverSelected() {
// ************************************************
  std::vector<char> Never(GetNEvents(), 0);
  std::size_t NNever = 0;
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    for (const auto& Cut : StoredSelection[static_cast<std::size_t>(MCEvents[iEvent].NominalSample)]) {
      if (IsKinematicParameterShifted(Cut.ParamToCutOnIt)) continue;
      const double Value = ReturnKinematicParameter(Cut.ParamToCutOnIt, static_cast<int>(iEvent));
      if (Value < Cut.LowerBound || Value >= Cut.UpperBound) {
        Never[iEvent] = 1;
        ++NNever;
        break;
      }
    }
  }
  if (NNever > 0) {
    MACH3LOG_INFO("{}: {} of {} events fail a cut no shift can change, so are never selected", GetName(), NNever, GetNEvents());
  }
  return Never;
}

// ************************************************
std::vector<double> SampleHandlerBase::CurrentShiftKey() const {
// ************************************************
  std::vector<double> Key;
  for (const auto& Shift : functional.shifts) Key.insert(Key.end(), Shift.par_vals.begin(), Shift.par_vals.end());
  return Key;
}

// ************************************************
void SampleHandlerBase::InvalidateReweightCache() {
// ************************************************
  if (Backend) Backend->Invalidate();
  WeightSnapshot.Invalidate();
}

// ************************************************
void SampleHandlerBase::EvaluateWeights() {
// ************************************************
  //HW Oscillation and spline weights depend only on non functional parameters. With a backend in
  //use, skip them when none of those has moved (every random walk sweep over the shifts). The
  //plain loop always evaluates, as it always has
  const bool MaySkip = Backend && !UseLegacyLoop && ParHandler != nullptr;
  if (MaySkip && WeightSnapshot.Matches(*ParHandler)) return;

  if (Oscillator) Oscillator->Evaluate();
  if (OscGradient) OscGradient->Evaluate();
  if (SplineHandler) SplineHandler->Evaluate();
  if (ParHandler != nullptr) WeightSnapshot.Take(*ParHandler);
}

// ************************************************
double SampleHandlerBase::CheckReweightBackend(const int NPoints, const double Tolerance, const double Width,
                                               const unsigned int Seed) {
// ************************************************
  if (!Backend || ParHandler == nullptr) return 0.;
  const int NPars = ParHandler->GetNumParams();
  std::vector<double> Original(static_cast<std::size_t>(NPars));
  for (int i = 0; i < NPars; ++i) Original[static_cast<std::size_t>(i)] = ParHandler->GetParProp(i);

  //HW W2 is normally filled only on the first reweight (unless UpdateW2), so force it on for the
  //comparison and put back what was there afterwards
  const bool FirstTimeW2Before = FirstTimeW2;
  const std::vector<double> W2Before = SampleHandler_array_w2;
  const auto Take = [this]() {
    ReweightSnapshot Snap;
    FirstTimeW2 = true;
    Reweight();
    Snap.MC = SampleHandler_array;
    Snap.W2 = SampleHandler_array_w2;
    Snap.Likelihood = GetLikelihood();
    if (Gradient) Snap.Gradient = GetGradLikelihood();
    return Snap;
  };

  std::mt19937 Rng(Seed);
  std::normal_distribution<double> Gaus(0., 1.);
  const char* const Moves[3] = {"shifts only", "everything but the shifts", "everything"};
  double Worst = 0.;
  for (int Point = 0; Point < NPoints; ++Point) {
    // Cycle through moving the shifts, the rest, and both, which takes the backend through its full
    // pass with the weights unchanged, its cached pass, and its full pass with new weights
    const int Move = Point % 3;
    for (int i = 0; i < NPars; ++i) {
      if (ParHandler->IsParameterFixed(i)) continue;
      const bool IsShift = (ParHandler->GetParamType(i) == SystType::kFunc);
      if ((Move == 0 && !IsShift) || (Move == 1 && IsShift)) continue;
      const double Value = ParHandler->GetParPreFit(i) + Width * ParHandler->GetDiagonalError(i) * Gaus(Rng);
      ParHandler->SetParProp(i, std::clamp(Value, ParHandler->GetLowerBound(i), ParHandler->GetUpperBound(i)));
    }
    const ReweightSnapshot Fast = Take();
    UseLegacyLoop = true;
    const ReweightSnapshot Reference = Take();
    UseLegacyLoop = false;

    int Where = -1, WhereW2 = -1, WhereGrad = -1;
    const double MCDifference = LargestRelativeDifference(Reference.MC, Fast.MC, 1e-6, Where);
    const double W2Difference = LargestRelativeDifference(Reference.W2, Fast.W2, 1e-6, WhereW2);
    const double LLHDifference = std::fabs(Fast.Likelihood - Reference.Likelihood) / std::max(std::fabs(Reference.Likelihood), 1e-300);
    // The gradient sums the same terms in a different order, so compare to its largest entry
    const double GradDifference = Reference.Gradient.empty() ? 0. :
        LargestRelativeDifference(Reference.Gradient, Fast.Gradient, 1., WhereGrad);
    const double PointWorst = std::max({MCDifference, W2Difference, LLHDifference, GradDifference});
    Worst = std::max(Worst, PointWorst);
    const std::string Line = fmt::format("{}: {} backend against the plain loop, point {} ({}): MC {:.2g}, W2 {:.2g}, "
                                         "-lnL {:.2g}{}", GetName(), Backend->Name(), Point + 1, Moves[Move], MCDifference,
                                         W2Difference, LLHDifference, Reference.Gradient.empty() ? std::string() :
                                         fmt::format(", gradient {:.2g}", GradDifference));
    if (PointWorst > Tolerance) MACH3LOG_ERROR("{}", Line);
    else MACH3LOG_INFO("{}", Line);
  }

  for (int i = 0; i < NPars; ++i) ParHandler->SetParProp(i, Original[static_cast<std::size_t>(i)]);
  FirstTimeW2 = FirstTimeW2Before;
  Reweight();
  SampleHandler_array_w2 = W2Before;
  return Worst;
}

// ************************************************
void SampleHandlerBase::SetComputeGradients(const bool Enable) {
// ************************************************
  if (!Enable) {
    OscGradient.reset();
    Gradient.reset();
    GradientCurrent = false;
    GradientCoverage.clear();
    if (SplineHandler) SplineHandler->SetComputeDerivatives(false);
    return;
  }
  if (ParHandler == nullptr) {
    MACH3LOG_ERROR("{} has no parameter handler to take gradients with respect to", GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  // Fail now rather than at the first GetGradLikelihood, hours into a fit
  GetTestStatDerivative(1., 1.);

  //HW Normalisation weights point at the parameter values themselves, so their derivative is 1
  std::vector<WeightDerivativeSource> Sources;
  WeightDerivativeSource Norms;
  Norms.Weights = ParHandler->RetPointer(0);
  Norms.N = static_cast<std::size_t>(ParHandler->GetNumParams());
  Sources.push_back(std::move(Norms));
  bool SplineDerivatives = false;
  if (SplineHandler) {
    SplineDerivatives = SplineHandler->SetComputeDerivatives(true);
    if (SplineDerivatives) {
      Sources.push_back(SplineHandler->GetWeightDerivativeSource());
    } else {
      MACH3LOG_WARN("{}: its splines give no derivatives, so spline parameters have no gradient", GetName());
    }
  }
  SetupOscillationGradient();
  if (OscGradient) Sources.push_back(OscGradient->Source());
  if (Backend) Backend->AttachDerivatives(Sources, ParHandler->GetNumParams());
  Gradient = std::make_unique<LikelihoodGradient>(std::move(Sources), ParHandler->GetNumParams(), MCEvents);
  GradientCurrent = false;
  // The spline derivatives only exist from the next evaluation on, so it mustn't be skipped
  WeightSnapshot.Invalidate();

  //HW Decided by the kind of parameter, not by which ones events happen to use: a normalisation no
  //event depends on has a derivative of exactly zero, which is what the gradient says. Functional
  //shifts move events between bins, so the binned likelihood is piecewise constant in them.
  //Oscillation parameters have one when OscillationGradient can reproduce NuOscillator
  GradientCoverage.assign(static_cast<std::size_t>(ParHandler->GetNumParams()), false);
  int NCovered = 0;
  for (int i = 0; i < ParHandler->GetNumParams(); ++i) {
    const SystType Type = ParHandler->GetParamType(i);
    const bool Covered = (Type == SystType::kNorm) || (Type == SystType::kSpline && SplineDerivatives) ||
                         (Type == SystType::kOsc && OscGradient != nullptr);
    GradientCoverage[static_cast<std::size_t>(i)] = Covered;
    NCovered += Covered ? 1 : 0;
  }
  MACH3LOG_INFO("{}: likelihood gradients on for {} of {} parameters. {}", GetName(), NCovered,
                ParHandler->GetNumParams(), Gradient->Summary());
}

// ************************************************
void SampleHandlerBase::SetupOscillationGradient() {
// ************************************************
  OscGradient.reset();
  if (!Oscillator || ParHandler == nullptr) return;
  if (!GetFromManager<bool>(SampleManager->raw()["NuOsc"]["AddOscPointer"], true, __FILE__ , __LINE__)) return;

  // CUDAProb3Linear's eight parameters, in its order, by the names NuOscillator knows them by
  const std::vector<std::string> Names = {"sin2_th12", "sin2_th23", "sin2_th13", "dm2_12", "dm2_23",
                                          "delta_cp", "path_length", "matter_density"};
  std::vector<int> Index(Names.size(), -1);
  for (const auto& Par : ParHandler->GetOscParsFromSampleName(GetSampleName(0))) {
    for (std::size_t k = 0; k < Names.size(); ++k) if (Par.NuOscName == Names[k]) Index[k] = Par.index;
  }
  std::vector<const M3::float_t*> Value(Names.size(), nullptr);
  std::vector<bool> Free(Names.size(), false);
  for (std::size_t k = 0; k < Names.size(); ++k) {
    if (Index[k] < 0) continue;
    Value[k] = ParHandler->RetPointer(Index[k]);
    Free[k] = !ParHandler->IsParameterFixed(Index[k]);
  }

  // Every event with a NuOscillator weight: neutral current events have none
  const auto Flavour = [](const int Pdg) -> signed char {
    switch (std::abs(Pdg)) { case 12: return 0; case 14: return 1; case 16: return 2; default: return -1; }
  };
  std::vector<OscillationGradient::EventWeight> Events;
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    const EventInfo& Ev = MCEvents[iEvent];
    if (Ev.isNC) continue;
    const signed char From = Flavour(Ev.nupdgUnosc), To = Flavour(Ev.nupdg);
    if (From < 0 || To < 0) continue;
    // NuOscillator is handed the energy as its FLOAT_T, see GetNuOscillatorPointers
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wuseless-cast"
    const double Energy = static_cast<double>(FLOAT_T(Ev.enu_true));
    #pragma GCC diagnostic pop
    Events.push_back({GetNuOscillatorPointers(static_cast<int>(iEvent)), Energy, Ev.nupdgUnosc < 0, From, To});
  }
  // The formulae are checked against NuOscillator's weights, which have to be up to date for that
  Oscillator->Evaluate();
  WeightSnapshot.Invalidate();
  const auto ConfigFile = Get<std::string>(SampleManager->raw()["NuOsc"]["NuOscConfigFile"], __FILE__ , __LINE__);
  OscGradient = OscillationGradient::Create(GetName(), ConfigFile, Index, Value, Free, Events);
}

// ************************************************
std::vector<double> SampleHandlerBase::GetGradLikelihood() const {
// ************************************************
  std::vector<double> Out;
  GetGradLikelihood(Out);
  return Out;
}

// ************************************************
void SampleHandlerBase::GetGradLikelihood(std::vector<double>& Out) const {
// ************************************************
  if (!Gradient) {
    MACH3LOG_ERROR("Asked {} for the likelihood gradient without SetComputeGradients(true)", GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  if (!GradientCurrent) {
    MACH3LOG_ERROR("Asked {} for the likelihood gradient, but the last reweight had gradients paused "
                   "or happened before they were switched on", GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  const int NBins = Binning->GetNBins();
  GradBinFactor.resize(static_cast<std::size_t>(NBins));
  for (int iBin = 0; iBin < NBins; ++iBin) {
    GradBinFactor[static_cast<std::size_t>(iBin)] = GetTestStatDerivative(SampleHandler_data[static_cast<std::size_t>(iBin)],
                                                                       SampleHandler_array[static_cast<std::size_t>(iBin)]);
  }
  if (Backend && !UseLegacyLoop) Backend->Gradient(GradBinFactor, Out);
  else Gradient->Compute(MCEvents, GradBinFactor, Out);
}

// ************************************************
const ParameterHandlerBase* SampleHandlerBase::GetGradParHandler() const {
// ************************************************
  return ParHandler;
}

// ************************************************
std::vector<bool> SampleHandlerBase::GetGradientCoverage() const {
// ************************************************
  return GradientCoverage;
}

// ************************************************
bool SampleHandlerBase::EnableSurrogate(const int Stride) {
// ************************************************
  Surrogate.On = (Backend != nullptr) && !UseLegacyLoop && Backend->EnableSurrogate(Stride);
  if (!Surrogate.On) {
    MACH3LOG_WARN("{}: no surrogate for shift moves, the {} doesn't provide one", GetName(),
                  Backend ? Backend->Name() + " backend" : std::string("plain event loop"));
    return false;
  }
  const std::size_t NBins = static_cast<std::size_t>(Binning->GetNBins());
  Surrogate.DeltaMC.assign(NBins, 0.);
  Surrogate.DeltaW2.assign(NBins, 0.);
  Surrogate.MoveW2 = UpdateW2;
  MACH3LOG_INFO("{}: surrogate for shift moves from every {}th live event, W2 {}", GetName(), Stride,
                Surrogate.MoveW2 ? "moving with the events" : "held at its first fill");
  return true;
}

// ************************************************
void SampleHandlerBase::AnchorSurrogate() {
// ************************************************
  if (!Surrogate.On) return;
  Surrogate.AnchorMC = SampleHandler_array;
  Surrogate.AnchorW2 = SampleHandler_array_w2;
  Backend->AnchorSurrogate();
}

// ************************************************
void SampleHandlerBase::ReweightSurrogate() {
// ************************************************
  if (!Surrogate.On) {
    MACH3LOG_ERROR("{}: ReweightSurrogate without a surrogate, EnableSurrogate first", GetName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  // As Reweight(), for the shifts only: the weights are the anchor's
  functional.update_vals();
  PrepFunctionalParameters();
  RefreshShiftTable();
  Selection = StoredSelection;

  std::fill(Surrogate.DeltaMC.begin(), Surrogate.DeltaMC.end(), 0.);
  std::fill(Surrogate.DeltaW2.begin(), Surrogate.DeltaW2.end(), 0.);
  Backend->SurrogateDelta(Surrogate.MoveW2, Surrogate.DeltaMC.data(), Surrogate.DeltaW2.data());

  //HW Clamped at zero, the same on the way back (ReverseSurrogateLikelihood), so the surrogate is
  //one function of its centre and the point it is evaluated at
  const int NBins = Binning->GetNBins();
  for (int b = 0; b < NBins; ++b) {
    const std::size_t i = static_cast<std::size_t>(b);
    SampleHandler_array[i] = std::max(0., Surrogate.AnchorMC[i] + Surrogate.DeltaMC[i]);
    if (Surrogate.MoveW2) SampleHandler_array_w2[i] = std::max(0., Surrogate.AnchorW2[i] + Surrogate.DeltaW2[i]);
  }
  // Events off the linear table were shifted to the proposal, and the prediction isn't exact
  KinematicsStale = true;
  GradientCurrent = false;
}

// ************************************************
double SampleHandlerBase::ReverseSurrogateLikelihood() const {
// ************************************************
  const int NBins = Binning->GetNBins();
  const double* Data = SampleHandler_data.data();
  const double* MC = SampleHandler_array.data();
  const double* W2 = SampleHandler_array_w2.data();
  const double* D = Surrogate.DeltaMC.data();
  const double* D2 = Surrogate.DeltaW2.data();
  const bool MoveW2 = Surrogate.MoveW2;
  double negLogL = 0.;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:negLogL)
  #endif
  for (int b = 0; b < NBins; ++b) {
    const double Pred = std::max(0., MC[b] - D[b]);
    const double Var = MoveW2 ? std::max(0., W2[b] - D2[b]) : W2[b];
    negLogL += GetTestStatLLH(Data[b], Pred, Var);
  }
  return negLogL;
}
