#include "Samples/LinearShiftBuilder.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace {
  /// Events measured at a time, see LinearShiftBuilder
  constexpr std::size_t BlockSize = 65536;
  /// Random points the model is checked at
  constexpr int NCheckPoints = 5;
  /// Relative agreement required at each of them
  constexpr double CheckTolerance = 1e-9;
}

// ************************************************
LinearShiftBuilder::LinearShiftBuilder(M3::detail::Functional& Shifts_, ReweightEventAccess& Generic_,
                                       const FastEventAccess& Fast_, const std::vector<char>& NeverSelected_,
                                       std::vector<double> Width_)
    : Shifts(Shifts_), Generic(Generic_), Fast(Fast_), NeverSelected(NeverSelected_), Width(std::move(Width_)),
      NShifts(Shifts_.raw_table.size()), NEvents(NeverSelected_.size()) {
// ************************************************
  Table.Shifts = &Shifts;
  Table.BaseValue.resize(NShifts);
  for (std::size_t s = 0; s < NShifts; ++s) Table.BaseValue[s] = Shifts.raw_table[s].value;
}

// ************************************************
LinearShiftTable LinearShiftBuilder::Build() {
// ************************************************
  Table.VarStart.assign(NEvents + 1, 0);
  for (std::size_t e = 0; e < NEvents; ++e) Table.VarStart[e + 1] = Table.VarStart[e] + static_cast<std::uint32_t>(NVariables(e));
  Table.Base.assign(Table.VarStart[NEvents], 0.);

  CountTerms();
  FillTerms();
  Dense.clear();
  Dense.shrink_to_fit();

  CheckAtRandomPoints();
  RestoreBase();
  return std::move(Table);
}

// ************************************************
std::size_t LinearShiftBuilder::NVariables(const std::size_t e) const {
// ************************************************
  return NeverSelected[e] ? std::size_t{0} : 1 + Fast.NShiftedCuts(static_cast<int>(e));
}

// ************************************************
void LinearShiftBuilder::MeasureBlock(const std::size_t First, const std::size_t Last) {
// ************************************************
  const int FirstInt = static_cast<int>(First), LastInt = static_cast<int>(Last);
  // Each event's variables at the base point
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = FirstInt; i < LastInt; ++i) {
    const std::size_t e = static_cast<std::size_t>(i);
    if (NeverSelected[e]) continue;
    Generic.ShiftEvent(i);
    for (std::size_t k = 0; k < NVariables(e); ++k) Table.Base[Table.VarStart[e] + k] = Read(e, k);
  }
  // One shift at a time, by one width: each variable's change is that shift's coefficient
  const std::size_t BlockVars = Table.VarStart[Last] - Table.VarStart[First];
  Dense.assign(BlockVars * NShifts, 0.);
  for (std::size_t sh = 0; sh < NShifts; ++sh) {
    Shifts.raw_table[sh].value = Table.BaseValue[sh] + Width[sh];
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = FirstInt; i < LastInt; ++i) {
      const std::size_t e = static_cast<std::size_t>(i);
      if (NeverSelected[e]) continue;
      bool Applies = false;
      for (const int iShift : Shifts.event_shifts[e]) Applies = Applies || (iShift == static_cast<int>(sh));
      if (!Applies) continue;
      Generic.ShiftEvent(i);
      for (std::size_t k = 0; k < NVariables(e); ++k) {
        const std::size_t Var = Table.VarStart[e] + k;
        Dense[(Var - Table.VarStart[First]) * NShifts + sh] = (Read(e, k) - Table.Base[Var]) / Width[sh];
      }
    }
    Shifts.raw_table[sh].value = Table.BaseValue[sh];
  }
}

// ************************************************
void LinearShiftBuilder::CountTerms() {
// ************************************************
  const std::size_t NVars = Table.VarStart[NEvents];
  Table.TermStart.assign(NVars + 1, 0);
  for (std::size_t First = 0; First < NEvents; First += BlockSize) {
    const std::size_t Last = std::min(NEvents, First + BlockSize);
    MeasureBlock(First, Last);
    for (std::size_t v = 0; v < Table.VarStart[Last] - Table.VarStart[First]; ++v) {
      std::uint32_t Count = 0;
      for (std::size_t sh = 0; sh < NShifts; ++sh) Count += (Dense[v * NShifts + sh] != 0.) ? 1 : 0;
      Table.TermStart[Table.VarStart[First] + v + 1] = Count;
    }
  }
  for (std::size_t v = 0; v < NVars; ++v) Table.TermStart[v + 1] += Table.TermStart[v];
}

// ************************************************
void LinearShiftBuilder::FillTerms() {
// ************************************************
  const std::size_t NTerms = Table.TermStart[Table.VarStart[NEvents]];
  Table.TermShift.resize(NTerms);
  Table.TermCoef.resize(NTerms);
  for (std::size_t First = 0; First < NEvents; First += BlockSize) {
    const std::size_t Last = std::min(NEvents, First + BlockSize);
    MeasureBlock(First, Last);
    for (std::size_t v = 0; v < Table.VarStart[Last] - Table.VarStart[First]; ++v) {
      std::uint32_t t = Table.TermStart[Table.VarStart[First] + v];
      for (std::size_t sh = 0; sh < NShifts; ++sh) {
        const double c = Dense[v * NShifts + sh];
        if (c == 0.) continue;
        Table.TermShift[t] = static_cast<std::uint16_t>(sh);
        Table.TermCoef[t] = c;
        ++t;
      }
    }
  }
}

// ************************************************
void LinearShiftBuilder::CheckAtRandomPoints() {
// ************************************************
  Table.Linear.assign(NEvents, 0);
  for (std::size_t e = 0; e < NEvents; ++e) Table.Linear[e] = NeverSelected[e] ? 0 : 1;
  std::mt19937 Rng(2026);
  std::normal_distribution<double> Gaus(0., 1.);
  const int NEventsInt = static_cast<int>(NEvents);
  for (int Point = 0; Point < NCheckPoints; ++Point) {
    for (std::size_t s = 0; s < NShifts; ++s) Shifts.raw_table[s].value = Table.BaseValue[s] + 2. * Width[s] * Gaus(Rng);
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < NEventsInt; ++i) {
      const std::size_t e = static_cast<std::size_t>(i);
      if (!Table.Linear[e]) continue;
      Generic.ShiftEvent(i);
      for (std::size_t k = 0; k < NVariables(e); ++k) {
        const double Expected = Read(e, k);
        if (!(std::fabs(Table.Value(e, k) - Expected) <= CheckTolerance * std::max(1., std::fabs(Expected)))) {
          Table.Linear[e] = 0;
          break;
        }
      }
    }
  }
}

// ************************************************
void LinearShiftBuilder::RestoreBase() {
// ************************************************
  for (std::size_t s = 0; s < NShifts; ++s) Shifts.raw_table[s].value = Table.BaseValue[s];
  const int NEventsInt = static_cast<int>(NEvents);
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < NEventsInt; ++i) Generic.ShiftEvent(i);
}
