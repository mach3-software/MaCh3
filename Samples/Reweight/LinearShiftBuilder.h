#pragma once

#include "Samples/Reweight/LinearShifts.h"
#include "Samples/Reweight/FastEventAccess.h"

#include <vector>

/// @brief Measures a LinearShiftTable from the experiment's own shifts, and checks it
///
/// @details Each shift is moved by its parameter's prior width with the others at their base values,
/// and every live event's variables are read before and after: the change per width is that shift's
/// coefficient. The events are measured a block at a time, so the dense coefficients exist for one
/// block only before the zeros are dropped, in two passes (count the terms, then fill them) so the
/// term arrays are allocated once at their final size.
///
/// The model is then compared with the experiment's shifts at random points a few widths out; an
/// event it doesn't hold for (a clamp, say) is marked to keep going through the experiment's code.
///
/// Build() leaves the shifts at their base values and every event shifted to them, as it found them.
///
/// @author Henry Wallace
class LinearShiftBuilder {
 public:
  /// @param Shifts The handler's shifts; their raw_table is moved during the measurement and put back
  /// @param Generic Shifts an event with the experiment's own code
  /// @param Fast Reads an event's variables (binning variable, then shifted cut variables)
  /// @param NeverSelected Per event, whether it can never be selected and so is left out
  /// @param Width Per shift, the step to measure it over: its parameter's prior width
  LinearShiftBuilder(M3::detail::Functional& Shifts, ReweightEventAccess& Generic, const FastEventAccess& Fast,
                     const std::vector<char>& NeverSelected, std::vector<double> Width);

  /// @brief Measure the table and check it, see the class description
  LinearShiftTable Build();

 private:
  /// @brief Number of variables of event e: none for an event never selected
  std::size_t NVariables(std::size_t e) const;
  /// @brief Variable k of event e as the event currently holds it
  double Read(std::size_t e, std::size_t k) const { return *Fast.VariablePointer(static_cast<int>(e), k); }
  /// @brief Fill #Dense for the events in [First, Last): each variable's change per width of each shift
  void MeasureBlock(std::size_t First, std::size_t Last);
  /// @brief Pass one: how many non-zero terms each variable has, into Table.TermStart
  void CountTerms();
  /// @brief Pass two: the terms themselves, in shift order
  void FillTerms();
  /// @brief Mark the events the linear model doesn't reproduce at random points
  void CheckAtRandomPoints();
  /// @brief Every shift back to its base value, with the events shifted to it
  void RestoreBase();

  /// The shifts, moved during the measurement
  M3::detail::Functional& Shifts;
  /// Shifts an event with the experiment's code
  ReweightEventAccess& Generic;
  /// Reads an event's variables
  const FastEventAccess& Fast;
  /// Per event, whether it can never be selected
  const std::vector<char>& NeverSelected;
  /// Per shift, the step it is measured over
  const std::vector<double> Width;
  /// Number of shifts
  const std::size_t NShifts;
  /// Number of events
  const std::size_t NEvents;
  /// The table being built
  LinearShiftTable Table;
  /// One block's coefficients, variable-major: Dense[v * NShifts + shift]
  std::vector<double> Dense;
};
