#pragma once

#include <cstddef>
#include <functional>
#include <map>
#include <vector>

namespace M3::detail {

/// @brief The functional shifts that apply to each event
/// @author Henry Wallace
///
/// @details Read like a vector of per event lists, as before: size(), and operator[] giving a list
/// to iterate over or index, with size() and empty(). Built with Add() while shifts are registered,
/// then Compact() stores each distinct list once and gives every event the index of its list.
/// Which shifts apply depends only on what kind of event it is, so a few hundred distinct lists
/// cover millions of events: for the DUNE FD that is about 11 MB in place of 394 MB of per event
/// vectors. Add() after Compact() is still allowed, and expands the lists again first.
class EventShiftLists {
 public:
  /// @brief One event's shifts, read only
  class Range {
   public:
    Range(const int* Begin_, const int* End_) : Begin(Begin_), End(End_) {}
    const int* begin() const { return Begin; }
    const int* end() const { return End; }
    std::size_t size() const { return static_cast<std::size_t>(End - Begin); }
    bool empty() const { return Begin == End; }
    int operator[](const std::size_t k) const { return Begin[k]; }
   private:
    const int* Begin;
    const int* End;
  };

  /// @brief Number of events
  std::size_t size() const { return Compacted ? ListOf.size() : Building.size(); }
  /// @brief One event's shifts
  Range operator[](const std::size_t iEvent) const {
    if (!Compacted) {
      const auto& L = Building[iEvent];
      return Range(L.data(), L.data() + L.size());
    }
    const std::size_t List = static_cast<std::size_t>(ListOf[iEvent]);
    return Range(ListShifts.data() + ListStart[List], ListShifts.data() + ListStart[List + 1]);
  }
  /// @brief Set the number of events, while registering
  void resize(const std::size_t NEvents) { Expand(); Building.resize(NEvents); }
  /// @brief Make shift iShift apply to event iEvent
  void Add(const std::size_t iEvent, const int iShift) { Expand(); Building[iEvent].push_back(iShift); }

  /// @brief Store each distinct list once, see the class description
  void Compact() {
    if (Compacted) return;
    std::map<std::vector<int>, int> Index;
    ListOf.resize(Building.size());
    ListStart.assign(1, 0);
    ListShifts.clear();
    for (std::size_t e = 0; e < Building.size(); ++e) {
      const auto Found = Index.emplace(Building[e], static_cast<int>(ListStart.size()) - 1);
      if (Found.second) {
        ListShifts.insert(ListShifts.end(), Building[e].begin(), Building[e].end());
        ListStart.push_back(static_cast<int>(ListShifts.size()));
      }
      ListOf[e] = Found.first->second;
    }
    std::vector<std::vector<int>>().swap(Building);
    Compacted = true;
  }
  /// @brief Number of distinct lists, once compacted
  std::size_t NLists() const { return Compacted ? ListStart.size() - 1 : Building.size(); }
  /// @brief Memory the lists hold, in bytes
  std::size_t Bytes() const {
    std::size_t B = sizeof(int) * (ListOf.capacity() + ListStart.capacity() + ListShifts.capacity());
    for (const auto& L : Building) B += sizeof(L) + sizeof(int) * L.capacity();
    return B;
  }

 private:
  /// @brief Go back to one vector per event, so lists can be added to again
  void Expand() {
    if (!Compacted) return;
    Building.resize(ListOf.size());
    for (std::size_t e = 0; e < ListOf.size(); ++e) {
      const Range R = (*this)[e];
      Building[e].assign(R.begin(), R.end());
    }
    ListOf.clear(); ListStart.clear(); ListShifts.clear();
    Compacted = false;
  }
  /// Whether the compact form is in use
  bool Compacted = false;
  /// One list per event, while registering
  std::vector<std::vector<int>> Building;
  /// Each event's list, once compacted
  std::vector<int> ListOf;
  /// List l's shifts are ListShifts[ListStart[l]] to ListShifts[ListStart[l+1]-1]
  std::vector<int> ListStart;
  /// See ListStart
  std::vector<int> ListShifts;
};

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
    /// HW The parameter handler indices of the parameters this shift consumes, in par_vals order
    std::vector<int> par_indices;
  };
  std::vector<Shift> shifts;
  /// HW Applies all of one event's shifts in a single typed loop over plain function pointers,
  /// set by the first single parameter registration. One indirect call per event then replaces a
  /// std::function call per shift. Used only while every shift has a plain function and all were
  /// registered on the same event vector, see apply_all_possible
  void (*apply_all)(const Functional &, void *, int) = nullptr;
  /// HW The experiment's event vector apply_all indexes into
  void *apply_all_events = nullptr;
  /// HW Cleared by any registration apply_all can't serve
  bool apply_all_possible = true;
  /// HW One shift's plain function and its parameter's current value, as apply_all reads them
  struct RawShift {
    /// The shift's plain function
    void (*raw)() = nullptr;
    /// Its parameter's current value, copied from par_vals[0] before each reweight
    double value = 0.;
  };
  /// HW One entry per shift, in the order of shifts, see RefreshShiftTable
  std::vector<RawShift> raw_table;
  /// HW Whether raw_table can be used, see RefreshShiftTable
  bool table_ready = false;
  void update_vals() {
    for (auto &shift : shifts) {
      shift.update_vals(shift.par_vals);
    }
  }
  /// The shifts that apply to each event, see EventShiftLists
  EventShiftLists event_shifts;
};

} // namespace M3::detail
