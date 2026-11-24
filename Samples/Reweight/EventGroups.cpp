#include "Samples/Reweight/EventGroups.h"

#include <algorithm>

// ************************************************
void EventGroups::Build(const ReweightModel& Model, const std::vector<int>& LiveBin) {
// ************************************************
  const int NBins = Model.NBins;
  // Counting sort by bin, visiting events in signature order so each bin comes out sorted by signature
  EventStart.assign(static_cast<std::size_t>(NBins) + 1, 0);
  for (const int Bin : LiveBin) if (Bin >= 0) ++EventStart[static_cast<std::size_t>(Bin) + 1];
  for (int b = 0; b < NBins; ++b) EventStart[static_cast<std::size_t>(b) + 1] += EventStart[static_cast<std::size_t>(b)];
  Ordered.resize(static_cast<std::size_t>(EventStart[static_cast<std::size_t>(NBins)]));
  Next.assign(EventStart.begin(), EventStart.end() - 1);
  for (const int i : Model.BySignature) {
    const int Bin = LiveBin[static_cast<std::size_t>(i)];
    if (Bin >= 0) Ordered[static_cast<std::size_t>(Next[static_cast<std::size_t>(Bin)]++)] = i;
  }

  BinStart.assign(static_cast<std::size_t>(NBins) + 1, 0);
  Sig.clear();
  Const.clear();
  Const2.clear();
  for (int b = 0; b < NBins; ++b) {
    BinStart[static_cast<std::size_t>(b)] = NGroups();
    int Begin = EventStart[static_cast<std::size_t>(b)];
    const int End = EventStart[static_cast<std::size_t>(b) + 1];
    while (Begin < End) {
      const int RunSig = Model.Sig[static_cast<std::size_t>(Ordered[static_cast<std::size_t>(Begin)])];
      int RunEnd = Begin + 1;
      while (RunEnd < End && Model.Sig[static_cast<std::size_t>(Ordered[static_cast<std::size_t>(RunEnd)])] == RunSig) ++RunEnd;
      AddRun(Model, Begin, RunEnd);
      Begin = RunEnd;
    }
  }
  BinStart[static_cast<std::size_t>(NBins)] = NGroups();

  BinOrder.clear();
  for (int b = 0; b < NBins; ++b) {
    if (BinStart[static_cast<std::size_t>(b) + 1] > BinStart[static_cast<std::size_t>(b)]) BinOrder.push_back(b);
  }
  std::stable_sort(BinOrder.begin(), BinOrder.end(), [this](const int a, const int b) {
    return (BinStart[static_cast<std::size_t>(a) + 1] - BinStart[static_cast<std::size_t>(a)]) >
           (BinStart[static_cast<std::size_t>(b) + 1] - BinStart[static_cast<std::size_t>(b)]);
  });
}

// ************************************************
void EventGroups::AddRun(const ReweightModel& Model, const int Begin, const int End) {
// ************************************************
  double Positive = 0., Positive2 = 0., Negative = 0., Negative2 = 0.;
  bool AnyPositive = false, AnyNegative = false;
  for (int k = Begin; k < End; ++k) {
    const double C = Model.Const[static_cast<std::size_t>(Ordered[static_cast<std::size_t>(k)])];
    // A zero constant makes the event's weight zero, which the per event loop skips anyway
    if (C > 0.) { Positive += C; Positive2 += C * C; AnyPositive = true; }
    else if (C < 0.) { Negative += C; Negative2 += C * C; AnyNegative = true; }
  }
  const int RunSig = Model.Sig[static_cast<std::size_t>(Ordered[static_cast<std::size_t>(Begin)])];
  if (AnyPositive) { Sig.push_back(RunSig); Const.push_back(Positive); Const2.push_back(Positive2); }
  if (AnyNegative) { Sig.push_back(RunSig); Const.push_back(Negative); Const2.push_back(Negative2); }
}
