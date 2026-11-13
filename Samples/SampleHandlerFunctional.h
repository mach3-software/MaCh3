#pragma once

#include <vector>
#include <functional>

namespace M3::detail {

/// @brief Helper object for storing/updating information related to functional
/// shift parameters
/// @author Luke Pickering
///
/// @details For each shift (which may consume multiple parameters) this class
///          stores the function that applies the shift to a given event index
///          given the current vector of parameter values
///          (Functional::Shift::apply), which is supplied by the an
///          experiment sample calling RegisterIndividualFunctionalParameter.
///          It also stores the last parameter values for all parameters
///          consumed by a given shift (Functional::Shift::par_vals) and
///          a function that will get the current parameter values from the
///          ParHandler object (Functional::Shift::update_vals).
///          For each event managed by the SampleHandlerBase instance, a
///          vector of indexes of functional shifts that are applicable is
///          stored (Functional::event_shifts). i.e. to iterate over the
///          relevant Functional::Shift instance for event ev_it you might do
///          for(auto & fs_it : functional.event_shifts[ev_it]){
///            // do something with functional.shifts[fs_it]
///          }
///          functional.update_vals() will update the current parameter values
///          for all registered shifts.
struct Functional {
  struct Shift {
    std::vector<double> par_vals;
    std::function<void(std::vector<double> &)> update_vals;
    std::function<void(std::vector<double> const &, const int)> apply;
    /// HW The experiment's shift as a plain function of one parameter, type erased, or nullptr
    /// when it was registered as anything else. See apply_all
    void (*raw)() = nullptr;
  };
  std::vector<Shift> shifts;
  /// HW Applies every shift of one event through the raw functions in a single typed loop, set by
  /// the first single parameter registration, so one indirect call per event replaces a
  /// std::function call per shift per event. Only used while every shift has a raw function and
  /// all of them were registered on the same event vector
  void (*apply_all)(const Functional &, void *, int) = nullptr;
  /// HW The experiment's event vector apply_all indexes into
  void *apply_all_events = nullptr;
  /// HW Cleared by any registration apply_all can't serve
  bool apply_all_possible = true;
  /// HW What apply_all reads, laid out flat so a shift costs one table entry rather than a walk
  /// through Shift, its std::functions and its par_vals vector
  struct RawShift {
    /// The shift's raw function
    void (*raw)() = nullptr;
    /// Its parameter's current value, copied from par_vals[0] before each full pass
    double value = 0.;
  };
  /// HW One entry per shift, same order as shifts
  std::vector<RawShift> raw_table;
  /// HW event_shifts flattened: event i's shifts are flat_index[flat_start[i]] to flat_index[flat_start[i+1]-1]
  std::vector<int> flat_start;
  /// HW See flat_start
  std::vector<int> flat_index;
  void update_vals() {
    for (auto &shift : shifts) {
      shift.update_vals(shift.par_vals);
    }
  }
  std::vector<std::vector<int>> event_shifts;
};

} // namespace M3::detail
