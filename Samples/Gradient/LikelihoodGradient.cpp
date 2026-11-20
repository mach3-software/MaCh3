#include "Samples/Gradient/LikelihoodGradient.h"

#include <sstream>

//HW GCC implements OpenMP's array reductions with alloca
#pragma GCC diagnostic ignored "-Walloca"

// ************************************************
LikelihoodGradient::LikelihoodGradient(std::vector<WeightDerivativeSource> Sources_, const int NPars_,
                                       const std::vector<EventInfo>& Events)
    : Sources(std::move(Sources_)), NPars(NPars_) {
// ************************************************
  const std::size_t NEvents = Events.size();
  EventBin.assign(NEvents, kNoBin);
  EventWeight.assign(NEvents, 0.);
  for (const auto& Event : Events) {
    for (const M3::float_t* Weight : Event.total_weight_pointers) {
      if (VisitWeightDerivatives(Sources, Weight, [](int, const M3::float_t*) {})) ++NResolved;
      else ++NUnresolved;
    }
  }
}

// ************************************************
void LikelihoodGradient::Compute(const std::vector<EventInfo>& Events,
                                 const std::vector<double>& BinFactor, std::vector<double>& Out) const {
// ************************************************
  Out.assign(static_cast<std::size_t>(NPars), 0.);
  double* Grad = Out.data();
  const int NEvents = static_cast<int>(Events.size());
  const int N = NPars;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:Grad[:N])
  #endif
  for (int iEvent = 0; iEvent < NEvents; ++iEvent) {
    const int Bin = EventBin[static_cast<std::size_t>(iEvent)];
    if (Bin == kNoBin) continue;
    const double Factor = BinFactor[static_cast<std::size_t>(Bin)];
    if (Factor == 0.) continue;
    const double FactorTimesWeight = Factor * EventWeight[static_cast<std::size_t>(iEvent)];
    for (const M3::float_t* Weight : Events[static_cast<std::size_t>(iEvent)].total_weight_pointers) {
      // Product rule: d(W)/dx = W * sum over its weights w of (dw/dx) / w. The event's weight is
      // positive, so none of its weights is zero
      VisitWeightDerivatives(Sources, Weight, [&](const int Param, const M3::float_t* DerivativePointer) {
        const double Derivative = (DerivativePointer == nullptr) ? 1. : static_cast<double>(*DerivativePointer);
        if (Derivative != 0.) Grad[Param] += FactorTimesWeight * Derivative / static_cast<double>(*Weight);
      });
    }
  }
}

// ************************************************
std::string LikelihoodGradient::Summary() const {
// ************************************************
  std::ostringstream Out;
  Out << NResolved << " event weights have a derivative and " << NUnresolved
      << " don't (oscillation, experiment specific or constant weights)";
  return Out.str();
}
