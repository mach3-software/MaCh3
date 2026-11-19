#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "Samples/SampleHandlerFunctional.h"

/// @brief Shifted kinematic variables as linear functions of the shift parameters
///
/// @details Many functional shifts are additive and linear in their parameter: an energy scale adds
/// the parameter times a fixed part of the event's energy, a resolution shift the parameter times
/// a fixed difference. Then every shifted variable of an event is
/// \f[
///   v = v_0 + \sum_k (x_k - x_k^0)\, c_k,
/// \f]
/// with \f$v_0\f$ its value at the base point \f$x^0\f$ and a coefficient \f$c_k\f$ for each shift
/// that moves it. A full pass can then work out the variables it cuts and bins on as plain
/// arithmetic, without calling the experiment's shift functions or writing to the event.
///
/// Only the non-zero terms are kept: most shifts move only some of an event's variables.
///
/// The table is measured from the experiment's own shifts (SampleHandlerBase::SetupLinearShifts) and
/// checked against them at random points; an event for which the model doesn't hold keeps going
/// through the experiment's code.
///
/// The variables of an event are its binning variable first, then the variables of its sample's
/// cuts on shifted variables, in the order FastEventAccess keeps them.
///
/// @author Henry Wallace
struct LinearShiftTable {
  /// The handler's shifts, whose raw_table holds their current values
  const M3::detail::Functional* Shifts = nullptr;
  /// Per shift, the value the base point was taken at
  std::vector<double> BaseValue;
  /// Per event, its first variable; one entry per event plus one
  std::vector<std::uint32_t> VarStart;
  /// Every variable at the base point
  std::vector<double> Base;
  /// Per variable, its first term; one entry per variable plus one
  std::vector<std::uint32_t> TermStart;
  /// The shift of every term
  std::vector<std::uint16_t> TermShift;
  /// The coefficient of every term
  std::vector<double> TermCoef;
  /// Per event, whether the linear model holds; the others go through the experiment's shifts
  std::vector<char> Linear;

  /// @brief Variable k of event e at the current shift values
  double Value(const std::size_t e, const std::size_t k) const {
    const std::size_t Var = VarStart[e] + k;
    const auto* Table = Shifts->raw_table.data();
    double v = Base[Var];
    for (std::uint32_t t = TermStart[Var]; t < TermStart[Var + 1]; ++t) {
      const std::uint16_t s = TermShift[t];
      v += (Table[s].value - BaseValue[s]) * TermCoef[t];
    }
    return v;
  }

  /// @brief Memory the table holds, in bytes
  std::size_t Bytes() const {
    return sizeof(double) * (BaseValue.size() + Base.size() + TermCoef.size()) +
           sizeof(std::uint32_t) * (VarStart.size() + TermStart.size()) +
           sizeof(std::uint16_t) * TermShift.size() + Linear.size();
  }
};
