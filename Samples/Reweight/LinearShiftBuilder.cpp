#include "Samples/Reweight/LinearShiftBuilder.h"

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
  // Where each event's variables start in the flat arrays: a prefix sum of how many each has
  Table.VarStart.assign(NEvents + 1, 0);
  for (std::size_t e = 0; e < NEvents; ++e) Table.VarStart[e + 1] = Table.VarStart[e] + static_cast<std::uint32_t>(NVariables(e));
  Table.Base.assign(Table.VarStart[NEvents], 0.);

  // Measure the coefficients, then check the model holds
  MeasureTerms();
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
      // Events this shift doesn't touch keep a coefficient of zero
      bool Applies = false;
      for (const int iShift : Shifts.event_shifts[e]) {
        if (iShift == static_cast<int>(sh)) { Applies = true; break; }
      }
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
void LinearShiftBuilder::MeasureTerms() {
// ************************************************
  const std::size_t NVars = Table.VarStart[NEvents];
  Table.TermStart.assign(NVars + 1, 0);
  Table.TermShift.clear();
  Table.TermCoef.clear();
  //HW One term per variable is the smallest a live variable can be, so this is the right order of
  //magnitude and usually the only allocation each array needs
  Table.TermShift.reserve(NVars);
  Table.TermCoef.reserve(NVars);
  for (std::size_t First = 0; First < NEvents; First += BlockSize) {
    const std::size_t Last = std::min(NEvents, First + BlockSize);
    MeasureBlock(First, Last);
    // Blocks ascend and so do the variables within one, so appending visits the terms in exactly
    // the order their offsets describe: a variable's terms start wherever the arrays have reached
    for (std::size_t v = 0; v < Table.VarStart[Last] - Table.VarStart[First]; ++v) {
      Table.TermStart[Table.VarStart[First] + v] = static_cast<std::uint32_t>(Table.TermShift.size());
      for (std::size_t sh = 0; sh < NShifts; ++sh) {
        const double c = Dense[v * NShifts + sh];
        if (c == 0.) continue;
        Table.TermShift.push_back(static_cast<std::uint16_t>(sh));
        Table.TermCoef.push_back(c);
      }
    }
  }
  // The last variable's end, which no later start fills in
  Table.TermStart[NVars] = static_cast<std::uint32_t>(Table.TermShift.size());
  // Growing the arrays leaves slack that the old count-then-fill pair didn't
  Table.TermShift.shrink_to_fit();
  Table.TermCoef.shrink_to_fit();
}

// ************************************************
void LinearShiftBuilder::CheckAtRandomPoints() {
// ************************************************
  // Start by trusting the model for every live event, then drop any it gets wrong at some point
  Table.Linear.assign(NEvents, 0);
  for (std::size_t e = 0; e < NEvents; ++e) Table.Linear[e] = NeverSelected[e] ? 0 : 1;
  std::mt19937 Rng(2026);
  std::normal_distribution<double> Gaus(0., 1.);
  const int NEventsInt = static_cast<int>(NEvents);
  for (int Point = 0; Point < NCheckPoints; ++Point) {
    // Every shift at once, a typical two widths from its base value
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
