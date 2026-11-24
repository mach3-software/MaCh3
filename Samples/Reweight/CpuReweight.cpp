#include "Samples/Reweight/CpuReweight.h"

#include <algorithm>
#include <sstream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Walloca"

// ************************************************
CpuReweight::CpuReweight(ReweightModel Model_, ReweightEventAccess& Access_)
    : Model(std::move(Model_)), Access(Access_) {
// ************************************************
  PoolValue.assign(Model.PoolPointer.size(), 0.);
  SigProduct.assign(static_cast<std::size_t>(Model.NSignatures()), 0.);
  LiveBin.assign(static_cast<std::size_t>(Model.NLive()), -1);
}

// ************************************************
std::string CpuReweight::Summary() const {
// ************************************************
  return Model.Summary();
}

// ************************************************
std::size_t CpuReweight::Bytes() const {
// ************************************************
  const auto Size = [](const auto& v) { return v.size() * sizeof(v[0]); };
  return Model.Bytes() + Size(PoolValue) + Size(SigProduct) + Size(LiveBin) + Size(Groups.BinStart) + Size(Groups.Sig) +
         Size(Groups.Const) + Size(Groups.Const2) + Size(Groups.BinOrder) + Size(PoolDerivStart) + Size(PoolDerivParam) +
         Size(PoolDerivative) + Size(SigDerivStart) + Size(SigDerivPool);
}

// ************************************************
void CpuReweight::Fill(const std::vector<double>& ShiftKey, const bool FillW2, double* MC, double* W2) {
// ************************************************
  ComputeSignatureProducts();
  //HW Exact comparison on purpose: a shift a sampler holds still comes back bit for bit the same,
  //and anything else has moved. A NaN never matches, so it always gets a full pass
  if (BinsValid && ShiftKey == BinKey) {
    CachedPass(FillW2, MC, W2);
    return;
  }
  BinKey = ShiftKey;
  FullPass(FillW2, MC, W2);
  BinsValid = true;
  GroupsValid = false;
}

// ************************************************
void CpuReweight::ComputeSignatureProducts() {
// ************************************************
  const int NPool = static_cast<int>(PoolValue.size());
  const M3::float_t* const* Pointer = Model.PoolPointer.data();
  double* Value = PoolValue.data();
  #ifdef MULTITHREAD
  #pragma omp parallel for simd
  #endif
  for (int k = 0; k < NPool; ++k) Value[k] = static_cast<double>(*Pointer[k]);

  const int NSig = Model.NSignatures();
  const int* Start = Model.SigStart.data();
  const int* Pool = Model.SigPool.data();
  double* Product = SigProduct.data();
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int s = 0; s < NSig; ++s) {
    double P = 1.;
    for (int k = Start[s]; k < Start[s + 1]; ++k) P *= Value[Pool[k]];
    Product[s] = P;
  }
}

// ************************************************
void CpuReweight::FullPass(const bool FillW2, double* MC, double* W2) {
// ************************************************
  const int NLive = Model.NLive();
  const int NBins = Model.NBins;
  const int* Event = Model.Event.data();
  const int* Sig = Model.Sig.data();
  const double* Const = Model.Const.data();
  const double* Product = SigProduct.data();
  int* Bin = LiveBin.data();
  //HW Small chunks dealt round robin: events differ in how many shifts they carry, and one block
  //per thread leaves threads waiting on the slowest. Static rather than dynamic, so each thread
  //always gets the same events in the same order, which makes each thread's own partial sum
  //repeatable. The reduction that combines them is not ordered, so this does not make the total,
  //or the chain, repeatable at more than one thread: see Doc/GradientSampling.md
  #ifdef MULTITHREAD
  #pragma omp parallel for schedule(static, 2048) reduction(+:MC[:NBins], W2[:NBins])
  #endif
  for (int i = 0; i < NLive; ++i) {
    //HW The bin is kept even for an event whose weight is zero now: it needn't be on the next pass
    const int GlobalBin = Access.ShiftCutAndBin(Event[i]);
    Bin[i] = GlobalBin;
    const double Weight = Const[i] * Product[Sig[i]];
    if (GlobalBin < 0 || Weight <= 0.) continue;
    MC[GlobalBin] += Weight;
    if (FillW2) W2[GlobalBin] += Weight * Weight;
  }
}

// ************************************************
bool CpuReweight::EnableSurrogate(const int Stride) {
// ************************************************
  SurStride = std::max(1, Stride);
  SurLive.clear();
  for (int i = 0; i < Model.NLive(); i += SurStride) SurLive.push_back(i);
  SurBin.assign(SurLive.size(), -1);
  SurWeight.assign(SurLive.size(), 0.);
  return true;
}

// ************************************************
void CpuReweight::AnchorSurrogate() {
// ************************************************
  //HW The products are recomputed rather than trusted, since a GPU backend's host copy may be
  //stale; LiveBin is the last fill's, a full pass's or the one a cached pass reused
  ComputeSignatureProducts();
  SyncLiveBin();
  const int N = static_cast<int>(SurLive.size());
  const int* Live = SurLive.data();
  const int* Bin = LiveBin.data();
  const int* Sig = Model.Sig.data();
  const double* Const = Model.Const.data();
  const double* Product = SigProduct.data();
  int* OutBin = SurBin.data();
  double* OutWeight = SurWeight.data();
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int k = 0; k < N; ++k) {
    const int i = Live[k];
    OutBin[k] = Bin[i];
    OutWeight[k] = Const[i] * Product[Sig[i]];
  }
}

// ************************************************
void CpuReweight::SurrogateDelta(const bool FillW2, double* DMC, double* DW2) {
// ************************************************
  const int N = static_cast<int>(SurLive.size());
  const int NBins = Model.NBins;
  const double Scale = static_cast<double>(SurStride);
  const int* Live = SurLive.data();
  const int* Old = SurBin.data();
  const double* W = SurWeight.data();
  const int* Event = Model.Event.data();
  #ifdef MULTITHREAD
  #pragma omp parallel for schedule(static, 512) reduction(+:DMC[:NBins], DW2[:NBins])
  #endif
  for (int k = 0; k < N; ++k) {
    if (W[k] <= 0.) continue;
    const int New = Access.ShiftCutAndBin(Event[Live[k]]);
    if (New == Old[k]) continue;
    const double w = Scale * W[k];
    const double w2 = Scale * W[k] * W[k];
    if (Old[k] >= 0) { DMC[Old[k]] -= w; if (FillW2) DW2[Old[k]] -= w2; }
    if (New >= 0)    { DMC[New] += w;    if (FillW2) DW2[New] += w2; }
  }
}

// ************************************************
void CpuReweight::EnsureGroups() {
// ************************************************
  if (GroupsValid) return;
  SyncLiveBin();
  Groups.Build(Model, LiveBin);
  GroupsValid = true;
  ++GroupsVersion;
}

// ************************************************
void CpuReweight::CachedPass(const bool FillW2, double* MC, double* W2) {
// ************************************************
  EnsureGroups();
  const int NUsed = static_cast<int>(Groups.BinOrder.size());
  const int* Order = Groups.BinOrder.data();
  const int* Start = Groups.BinStart.data();
  const int* Sig = Groups.Sig.data();
  const double* Const = Groups.Const.data();
  const double* Const2 = Groups.Const2.data();
  const double* Product = SigProduct.data();
  //HW One bin per iteration, so each bin is written by exactly one thread: no reduction copies.
  //Dynamic, since bins differ a lot in how many groups they hold, and the busiest are handed out first
  #ifdef MULTITHREAD
  #pragma omp parallel for schedule(dynamic, 1)
  #endif
  for (int n = 0; n < NUsed; ++n) {
    const int b = Order[n];
    double Sum = 0., Sum2 = 0.;
    for (int g = Start[b]; g < Start[b + 1]; ++g) {
      const double P = Product[Sig[g]];
      const double Weight = Const[g] * P;
      if (Weight <= 0.) continue;
      Sum += Weight;
      Sum2 += Const2[g] * P * P;
    }
    MC[b] += Sum;
    if (FillW2) W2[b] += Sum2;
  }
}

// ************************************************
void CpuReweight::AttachDerivatives(const std::vector<WeightDerivativeSource>& Sources, const int NPars_) {
// ************************************************
  NPars = NPars_;
  const std::size_t NPool = Model.PoolPointer.size();
  PoolDerivStart.assign(1, 0);
  PoolDerivParam.clear();
  PoolDerivative.clear();
  for (std::size_t k = 0; k < NPool; ++k) {
    VisitWeightDerivatives(Sources, Model.PoolPointer[k], [&](const int Param, const M3::float_t* Derivative) {
      PoolDerivParam.push_back(Param);
      PoolDerivative.push_back(Derivative);
    });
    PoolDerivStart.push_back(static_cast<int>(PoolDerivParam.size()));
  }
  const int NSig = Model.NSignatures();
  SigDerivStart.assign(1, 0);
  SigDerivPool.clear();
  for (int s = 0; s < NSig; ++s) {
    for (int k = Model.SigStart[static_cast<std::size_t>(s)]; k < Model.SigStart[static_cast<std::size_t>(s) + 1]; ++k) {
      const int Entry = Model.SigPool[static_cast<std::size_t>(k)];
      if (PoolDerivStart[static_cast<std::size_t>(Entry) + 1] > PoolDerivStart[static_cast<std::size_t>(Entry)]) SigDerivPool.push_back(Entry);
    }
    SigDerivStart.push_back(static_cast<int>(SigDerivPool.size()));
  }
}

// ************************************************
void CpuReweight::Gradient(const std::vector<double>& BinFactor, std::vector<double>& Out) {
// ************************************************
  EnsureGroups();
  // Each signature's groups, rebuilt with the groups (once per full pass, not per gradient)
  if (IndexedGroupsVersion != GroupsVersion || GradSigGroupStart.empty()) {
    IndexGroups(GradGroupBin, GradSigGroupStart, GradSigGroup);
    IndexedGroupsVersion = GroupsVersion;
  }
  // One pass over the signatures, in parallel: each signature's share of the bin factors (over its
  // groups that contribute, weight > 0), then its derivative weights, d(W)/dx = W * (dw/dx) / w for
  // the weight w of x. Each signature reads only its own groups, so there is no serial pass and no
  // race, and the static schedule keeps the sums the same from run to run
  const int NSig = Model.NSignatures();
  Out.assign(static_cast<std::size_t>(NPars), 0.);
  double* Grad = Out.data();
  const int N = NPars;
  const int* Start = SigDerivStart.data();
  const int* Pool = SigDerivPool.data();
  const int* DerivStart = PoolDerivStart.data();
  const int* Param = PoolDerivParam.data();
  const M3::float_t* const* Derivative = PoolDerivative.data();
  const double* Value = PoolValue.data();
  const double* Product = SigProduct.data();
  const int* GroupStart = GradSigGroupStart.data();
  const int* SigGroups = GradSigGroup.data();
  const int* BinOf = GradGroupBin.data();
  const double* Const = Groups.Const.data();
  const double* Factor = BinFactor.data();
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:Grad[:N])
  #endif
  for (int s = 0; s < NSig; ++s) {
    if (Start[s] == Start[s + 1]) continue;
    const double P = Product[s];
    double Share = 0.;
    for (int k = GroupStart[s]; k < GroupStart[s + 1]; ++k) {
      const int g = SigGroups[k];
      const double Weight = Const[g] * P;
      if (Weight > 0.) Share += Factor[BinOf[g]] * Weight;
    }
    if (Share == 0.) continue;
    for (int k = Start[s]; k < Start[s + 1]; ++k) {
      const int Entry = Pool[k];
      for (int j = DerivStart[Entry]; j < DerivStart[Entry + 1]; ++j) {
        const double D = (Derivative[j] == nullptr) ? 1. : static_cast<double>(*Derivative[j]);
        if (D != 0.) Grad[Param[j]] += Share * D / Value[Entry];
      }
    }
  }
}

// ************************************************
void CpuReweight::IndexGroups(std::vector<int>& GroupBin, std::vector<int>& SigGroupStart, std::vector<int>& SigGroup) const {
// ************************************************
  const int NGroups = Groups.NGroups();
  GroupBin.assign(static_cast<std::size_t>(NGroups), 0);
  for (int b = 0; b < Model.NBins; ++b) {
    for (int g = Groups.BinStart[static_cast<std::size_t>(b)]; g < Groups.BinStart[static_cast<std::size_t>(b) + 1]; ++g) {
      GroupBin[static_cast<std::size_t>(g)] = b;
    }
  }
  const int NSig = Model.NSignatures();
  SigGroupStart.assign(static_cast<std::size_t>(NSig) + 1, 0);
  for (const int s : Groups.Sig) ++SigGroupStart[static_cast<std::size_t>(s) + 1];
  for (int s = 0; s < NSig; ++s) SigGroupStart[static_cast<std::size_t>(s) + 1] += SigGroupStart[static_cast<std::size_t>(s)];
  std::vector<int> Next(SigGroupStart.begin(), SigGroupStart.end() - 1);
  SigGroup.assign(static_cast<std::size_t>(NGroups), 0);
  for (int g = 0; g < NGroups; ++g) {
    SigGroup[static_cast<std::size_t>(Next[static_cast<std::size_t>(Groups.Sig[static_cast<std::size_t>(g)])]++)] = g;
  }
}

// ************************************************
void CpuReweight::ParameterTerms(std::vector<int>& ParStart, std::vector<int>& TermSig, std::vector<int>& TermDeriv,
                                 std::vector<int>& TermPool) const {
// ************************************************
  // Count each parameter's terms, then place them, signatures in order within each parameter
  ParStart.assign(static_cast<std::size_t>(NPars) + 1, 0);
  const int NSig = Model.NSignatures();
  const auto ForEachTerm = [&](auto&& Visit) {
    for (int s = 0; s < NSig; ++s) {
      for (int k = SigDerivStart[static_cast<std::size_t>(s)]; k < SigDerivStart[static_cast<std::size_t>(s) + 1]; ++k) {
        const int Entry = SigDerivPool[static_cast<std::size_t>(k)];
        for (int j = PoolDerivStart[static_cast<std::size_t>(Entry)]; j < PoolDerivStart[static_cast<std::size_t>(Entry) + 1]; ++j) {
          Visit(PoolDerivParam[static_cast<std::size_t>(j)], s, j, Entry);
        }
      }
    }
  };
  ForEachTerm([&](const int Par, int, int, int) { ++ParStart[static_cast<std::size_t>(Par) + 1]; });
  for (int p = 0; p < NPars; ++p) ParStart[static_cast<std::size_t>(p) + 1] += ParStart[static_cast<std::size_t>(p)];
  const std::size_t NTerms = static_cast<std::size_t>(ParStart[static_cast<std::size_t>(NPars)]);
  TermSig.assign(NTerms, 0);
  TermDeriv.assign(NTerms, 0);
  TermPool.assign(NTerms, 0);
  std::vector<int> Next(ParStart.begin(), ParStart.end() - 1);
  ForEachTerm([&](const int Par, const int s, const int j, const int Entry) {
    const std::size_t t = static_cast<std::size_t>(Next[static_cast<std::size_t>(Par)]++);
    TermSig[t] = s;
    TermDeriv[t] = j;
    TermPool[t] = Entry;
  });
}

#pragma GCC diagnostic pop
