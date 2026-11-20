#include "Samples/Reweight/ReweightModel.h"

#include <sstream>
#include <unordered_map>

namespace {

/// @brief FNV-1a over a list of pool indices
struct IndexListHash {
  std::size_t operator()(const std::vector<int>& List) const noexcept {
    std::size_t Hash = 1469598103934665603ULL;
    for (const int Index : List) {
      Hash ^= static_cast<std::size_t>(static_cast<unsigned int>(Index));
      Hash *= 1099511628211ULL;
    }
    return Hash;
  }
};

/// @brief Order the live events by signature with a counting sort
std::vector<int> OrderBySignature(const std::vector<int>& Sig, const int NSignatures) {
  std::vector<int> Start(static_cast<std::size_t>(NSignatures) + 1, 0);
  for (const int s : Sig) ++Start[static_cast<std::size_t>(s) + 1];
  for (int s = 0; s < NSignatures; ++s) Start[static_cast<std::size_t>(s) + 1] += Start[static_cast<std::size_t>(s)];
  std::vector<int> Ordered(Sig.size());
  for (std::size_t i = 0; i < Sig.size(); ++i) Ordered[static_cast<std::size_t>(Start[static_cast<std::size_t>(Sig[i])]++)] = static_cast<int>(i);
  return Ordered;
}

} // namespace

// ************************************************
std::string ReweightModel::Summary() const {
// ************************************************
  std::ostringstream Out;
  Out << NLive() << " of " << NEvents << " events can be selected; they share " << NSignatures()
      << " signatures over " << PoolPointer.size() << " distinct varying weights";
  return Out.str();
}

// ************************************************
ReweightModel BuildReweightModel(const std::vector<EventInfo>& Events, const std::vector<std::uint16_t>& ConstStart,
                                 const std::vector<std::uint16_t>& ConstEnd, const std::vector<char>& NeverSelected,
                                 const int NBins) {
// ************************************************
  ReweightModel Model;
  Model.NEvents = static_cast<int>(Events.size());
  Model.NBins = NBins;
  Model.SigStart.push_back(0);

  std::unordered_map<const M3::float_t*, int> PoolIndex;
  std::unordered_map<std::vector<int>, int, IndexListHash> SignatureIndex;
  std::vector<int> Varying;
  for (int iEvent = 0; iEvent < Model.NEvents; ++iEvent) {
    if (NeverSelected[static_cast<std::size_t>(iEvent)]) continue;
    const auto& Pointers = Events[static_cast<std::size_t>(iEvent)].total_weight_pointers;
    const int First = ConstStart[static_cast<std::size_t>(iEvent)];
    const int Last = ConstEnd[static_cast<std::size_t>(iEvent)];
    double Const = 1.;
    Varying.clear();
    for (int k = 0; k < static_cast<int>(Pointers.size()); ++k) {
      const M3::float_t* Pointer = Pointers[static_cast<std::size_t>(k)];
      if (k >= First && k < Last) {
        Const *= static_cast<double>(*Pointer);
        continue;
      }
      const auto Found = PoolIndex.emplace(Pointer, static_cast<int>(Model.PoolPointer.size()));
      if (Found.second) Model.PoolPointer.push_back(Pointer);
      Varying.push_back(Found.first->second);
    }
    const auto Found = SignatureIndex.emplace(Varying, Model.NSignatures());
    if (Found.second) {
      Model.SigPool.insert(Model.SigPool.end(), Varying.begin(), Varying.end());
      Model.SigStart.push_back(static_cast<int>(Model.SigPool.size()));
    }
    Model.Event.push_back(iEvent);
    Model.Sig.push_back(Found.first->second);
    Model.Const.push_back(Const);
  }
  Model.BySignature = OrderBySignature(Model.Sig, Model.NSignatures());
  return Model;
}
