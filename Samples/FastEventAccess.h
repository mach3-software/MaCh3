#pragma once

#include <vector>

#include "Samples/ReweightBackend.h"
#include "Samples/SampleStructs.h"
#include "Samples/LinearShifts.h"
#include "Samples/SweepTables.h"

/// @brief Cuts and binning of shifted events read through plain arrays
///
/// @details Wraps the generic ReweightEventAccess (which goes through the experiment's
/// ReturnKinematicParameter and the BinningHandler) and replaces two of its three calls when it can:
///
/// - **Cuts.** A live event passes every cut on a variable no shift changes; that is what made it
///   live. So only the cuts on shifted variables need checking again, and those are read through a
///   pointer per event and cut, asked for once from the experiment and checked against
///   ReturnKinematicParameter for every event.
/// - **Binning.** For a sample with uniform one dimensional binning, each event's bin is found with
///   the sample binning's own FindBin, from a pointer to its variable and its nominal bin, without
///   the per event vectors FindGlobalBin goes through. Checked against FindGlobalBin for every event.
///
/// The shifts themselves are still applied by the experiment's own code, so the shifted variables
/// stored in the events are exactly what they would be after the plain event loop.
///
/// Each replacement is used only if it was set up and checked; otherwise the generic call is made.
///
/// @author Henry Wallace
class FastEventAccess final : public ReweightEventAccess {
 public:
  /// @param Generic The access to fall back on
  explicit FastEventAccess(ReweightEventAccess& Generic_) : Generic(Generic_) {}

  /// @copydoc ReweightEventAccess::ShiftEvent
  void ShiftEvent(const int iEvent) override { Generic.ShiftEvent(iEvent); }

  /// @copydoc ReweightEventAccess::PassesCuts
  /// @details Only valid for live events, see the class description
  bool PassesCuts(const int iEvent) override {
    if (!CutsReady) return Generic.PassesCuts(iEvent);
    const std::size_t e = static_cast<std::size_t>(iEvent);
    const std::vector<KinematicCut>& Cuts = ShiftedCuts[static_cast<std::size_t>(EventSample[e])];
    const double* const* Value = CutValue.data() + CutStart[e];
    const std::size_t NCuts = Cuts.size();
    for (std::size_t k = 0; k < NCuts; ++k) {
      const double v = *Value[k];
      if (v < Cuts[k].LowerBound || v >= Cuts[k].UpperBound) return false;
    }
    return true;
  }

  /// @copydoc ReweightEventAccess::GlobalBinOf
  int GlobalBinOf(const int iEvent) override {
    if (!BinsReady) return Generic.GlobalBinOf(iEvent);
    const std::size_t e = static_cast<std::size_t>(iEvent);
    return FlatBin(*SampleBinning[static_cast<std::size_t>(EventSample[e])], *BinVar[e], NomBin[e]);
  }

  /// @brief Global bin of a value in a uniform one dimensional sample binning, from its nominal bin
  static int FlatBin(const SampleBinningInfo& Info, const double Var, const int Nominal) {
    const int Bin = Info.FindBin(0, Var, Nominal);
    if (Bin < 0) return M3::UnderOverFlowBin;
    return Bin * Info.Strides[0] + Info.GlobalOffset;
  }

  /// @copydoc ReweightEventAccess::ShiftCutAndBin
  /// @details With linear shifts in use, an event the linear model holds for isn't shifted at all:
  /// its variables are worked out from the table, and the event itself is left as it was
  int ShiftCutAndBin(const int iEvent) override {
    const std::size_t e = static_cast<std::size_t>(iEvent);
    if (LinearReady && Linear.Linear[e]) return LinearBin(e);
    ShiftEvent(iEvent);
    if (!PassesCuts(iEvent)) return -1;
    const int Bin = GlobalBinOf(iEvent);
    return (Bin > M3::UnderOverFlowBin) ? Bin : -1;
  }

  /// @brief Number of cuts on shifted variables an event's sample has
  std::size_t NShiftedCuts(const int iEvent) const {
    return ShiftedCuts[static_cast<std::size_t>(EventSample[static_cast<std::size_t>(iEvent)])].size();
  }
  /// @brief Variable k of an event, as LinearShiftTable numbers them: its binning variable for k = 0,
  /// then its shifted cuts' variables. Needs both the cut pointers and the flat binning
  const double* VariablePointer(const int iEvent, const std::size_t k) const {
    const std::size_t e = static_cast<std::size_t>(iEvent);
    return (k == 0) ? BinVar[e] : CutValue[CutStart[e] + k - 1];
  }

  /// @brief Work out the shifted variables of the events the table holds for, instead of shifting them
  void UseLinearShifts(LinearShiftTable Table) {
    Linear = std::move(Table);
    LinearReady = true;
  }
  /// @brief Whether linear shifts are in use
  bool UsesLinearShifts() const { return LinearReady; }

  /// @brief Use pointers for the cuts on shifted variables
  /// @param Sample Every event's sample
  /// @param SampleShiftedCuts Per sample, its cuts on shifted variables, in the order the pointers follow
  /// @param Start Event e's pointers are Value[Start[e]] onwards, one per cut of its sample
  /// @param Value See Start
  void UseCutPointers(std::vector<int> Sample, std::vector<std::vector<KinematicCut>> SampleShiftedCuts,
                      std::vector<std::size_t> Start, std::vector<const double*> Value) {
    EventSample = std::move(Sample);
    ShiftedCuts = std::move(SampleShiftedCuts);
    CutStart = std::move(Start);
    CutValue = std::move(Value);
    CutsReady = true;
  }

  /// @brief Use the flat lookup for the binning
  /// @param Sample Every event's sample, as given to UseCutPointers
  /// @param Binning Per sample, its binning, which has to be uniform and one dimensional
  /// @param Var Every event's binning variable
  /// @param Nominal Every event's nominal bin
  void UseFlatBinning(std::vector<int> Sample, std::vector<const SampleBinningInfo*> Binning,
                      std::vector<const double*> Var, std::vector<int> Nominal) {
    EventSample = std::move(Sample);
    SampleBinning = std::move(Binning);
    BinVar = std::move(Var);
    NomBin = std::move(Nominal);
    BinsReady = true;
  }

  /// @brief The tables a GPU full pass needs, see SweepTables
  /// @param LiveEvents Handler index of every live event, in the backend's order
  /// @param Out Filled on success
  /// @return false, leaving Out unspecified, unless linear shifts, cut pointers and the flat binning
  /// are all in use and the linear model holds for every live event
  bool ExportSweepTables(const std::vector<int>& LiveEvents, SweepTables& Out) const {
    if (!LinearReady || !CutsReady || !BinsReady) return false;
    Out = SweepTables();
    Out.NLive = static_cast<int>(LiveEvents.size());
    Out.NShifts = static_cast<int>(Linear.BaseValue.size());
    Out.VarStart.push_back(0);
    Out.TermStart.push_back(0);
    for (const int iEvent : LiveEvents) {
      const std::size_t e = static_cast<std::size_t>(iEvent);
      if (!Linear.Linear[e]) return false;
      const int Sample = EventSample[e];
      Out.Sample.push_back(Sample);
      const std::size_t NVars = 1 + ShiftedCuts[static_cast<std::size_t>(Sample)].size();
      for (std::size_t k = 0; k < NVars; ++k) {
        const std::size_t Var = Linear.VarStart[e] + k;
        Out.Base.push_back(Linear.Base[Var]);
        for (std::uint32_t t = Linear.TermStart[Var]; t < Linear.TermStart[Var + 1]; ++t) {
          Out.TermShift.push_back(Linear.TermShift[t]);
          Out.TermCoef.push_back(Linear.TermCoef[t]);
        }
        Out.TermStart.push_back(static_cast<std::uint32_t>(Out.TermShift.size()));
      }
      Out.VarStart.push_back(static_cast<std::uint32_t>(Out.Base.size()));
    }
    Out.CutStart.push_back(0);
    Out.EdgeStart.push_back(0);
    for (std::size_t s = 0; s < ShiftedCuts.size(); ++s) {
      for (const KinematicCut& Cut : ShiftedCuts[s]) {
        Out.CutLo.push_back(Cut.LowerBound);
        Out.CutHi.push_back(Cut.UpperBound);
      }
      Out.CutStart.push_back(static_cast<int>(Out.CutLo.size()));
      const SampleBinningInfo* Info = (s < SampleBinning.size()) ? SampleBinning[s] : nullptr;
      if (Info == nullptr || Info->BinEdges.size() != 1 || Info->Strides.empty()) return false;
      Out.Edges.insert(Out.Edges.end(), Info->BinEdges[0].begin(), Info->BinEdges[0].end());
      Out.EdgeStart.push_back(static_cast<int>(Out.Edges.size()));
      Out.Stride.push_back(Info->Strides[0]);
      Out.Offset.push_back(Info->GlobalOffset);
    }
    return true;
  }
  /// @brief The shifts' current offsets from the linear table's base point, as Value() uses them
  void ShiftDeltas(std::vector<double>& Out) const {
    const std::size_t N = Linear.BaseValue.size();
    Out.resize(N);
    const auto* Table = Linear.Shifts->raw_table.data();
    for (std::size_t s = 0; s < N; ++s) Out[s] = Table[s].value - Linear.BaseValue[s];
  }

  /// @brief Memory the lookups hold, in bytes (the linear shift table's included)
  std::size_t Bytes() const {
    return sizeof(int) * (EventSample.size() + NomBin.size()) + sizeof(std::size_t) * CutStart.size() +
           sizeof(const double*) * (CutValue.size() + BinVar.size()) + (LinearReady ? Linear.Bytes() : 0);
  }
  /// @brief Whether the cut pointers are in use
  bool UsesCutPointers() const { return CutsReady; }
  /// @brief Whether the flat binning is in use
  bool UsesFlatBinning() const { return BinsReady; }

 private:
  /// See the constructor
  ReweightEventAccess& Generic;
  /// Every event's sample
  std::vector<int> EventSample;

  /// Whether the cut pointers are set up
  bool CutsReady = false;
  /// Per sample, its cuts on shifted variables
  std::vector<std::vector<KinematicCut>> ShiftedCuts;
  /// Where every event's pointers start in CutValue
  std::vector<std::size_t> CutStart;
  /// Pointers to the variables of every event's shifted cuts
  std::vector<const double*> CutValue;

  /// Whether the flat binning is set up
  bool BinsReady = false;
  /// Per sample, its binning
  std::vector<const SampleBinningInfo*> SampleBinning;
  /// Every event's binning variable
  std::vector<const double*> BinVar;
  /// Every event's nominal bin
  std::vector<int> NomBin;

  /// Whether linear shifts are in use
  bool LinearReady = false;
  /// See UseLinearShifts
  LinearShiftTable Linear;

  /// @brief Cut and bin an event from the variables the linear model gives
  int LinearBin(const std::size_t e) const {
    const std::vector<KinematicCut>& Cuts = ShiftedCuts[static_cast<std::size_t>(EventSample[e])];
    for (std::size_t k = 0; k < Cuts.size(); ++k) {
      const double v = Linear.Value(e, k + 1);
      if (v < Cuts[k].LowerBound || v >= Cuts[k].UpperBound) return -1;
    }
    const int Bin = FlatBin(*SampleBinning[static_cast<std::size_t>(EventSample[e])], Linear.Value(e, 0), NomBin[e]);
    return (Bin > M3::UnderOverFlowBin) ? Bin : -1;
  }
};
