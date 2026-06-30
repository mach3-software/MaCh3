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
  };
  std::vector<Shift> shifts;
  void update_vals() {
    for (auto &shift : shifts) {
      shift.update_vals(shift.par_vals);
    }
  }
  std::vector<std::vector<int>> event_shifts;
};

} // namespace M3::detail
