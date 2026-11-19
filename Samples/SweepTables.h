#pragma once

#include <cstdint>
#include <vector>

/// @brief Everything a full pass needs to cut and bin the live events itself, as plain arrays
///
/// @details Filled by FastEventAccess::ExportSweepTables when every live event's shifted variables
/// are linear in the shifts (LinearShiftTable), its shifted cuts are plain ranges and its sample's
/// binning is uniform and one dimensional. A GPU backend then does a full pass without calling the
/// experiment's code: variable k of live event i is
///   Base[VarStart[i] + k] + sum over its terms t of Delta[TermShift[t]] * TermCoef[t],
/// with Delta the shifts' current offsets from the table's base point, added in the table's order;
/// variable 0 is the binning variable and variable 1 + c the variable of cut c of its sample.
///
/// Plain standard types only, so device code can use it without MaCh3's headers.
///
/// @author Henry Wallace
struct SweepTables {
  /// Number of live events
  int NLive = 0;
  /// Number of shifts the terms index
  int NShifts = 0;

  /// Sample of each live event
  std::vector<int> Sample;
  /// Live event i's variables are Base[VarStart[i]] onwards; NLive + 1 entries
  std::vector<std::uint32_t> VarStart;
  /// Every variable at the table's base point
  std::vector<double> Base;
  /// Variable v's terms are TermStart[v] onwards; one entry per variable plus one
  std::vector<std::uint32_t> TermStart;
  /// The shift of every term
  std::vector<std::uint16_t> TermShift;
  /// The coefficient of every term
  std::vector<double> TermCoef;

  /// Sample s's cuts on shifted variables are CutLo/CutHi[CutStart[s]] onwards: pass if Lo <= v < Hi
  std::vector<int> CutStart;
  /// See CutStart
  std::vector<double> CutLo;
  /// See CutStart
  std::vector<double> CutHi;

  /// Sample s's bin edges are Edges[EdgeStart[s]] to Edges[EdgeStart[s+1]-1]
  std::vector<int> EdgeStart;
  /// See EdgeStart
  std::vector<double> Edges;
  /// Global bin = bin within the sample * Stride[s] + Offset[s]
  std::vector<int> Stride;
  /// See Stride
  std::vector<int> Offset;
};
