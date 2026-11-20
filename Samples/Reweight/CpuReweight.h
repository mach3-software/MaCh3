#pragma once

#include "Samples/Reweight/EventGroups.h"
#include "Samples/Reweight/ReweightBackend.h"
#include "Samples/Reweight/ReweightModel.h"

/// @brief The reweight on the CPU, from a ReweightModel
///
/// @details Two kinds of pass, chosen by whether the shifts moved since the last full one:
/// - **Full pass** (the shifts moved): every live event is shifted, cut and binned by the
///   experiment's own code, and weighed as its constant times its signature's product.
/// - **Cached pass** (they didn't): every event's bin is what the last full pass found, so the
///   events are summed in their groups, bin by bin, one thread per bin. Nothing is shifted, cut or
///   binned, and nothing runs per event.
///
/// Every pass starts by reading the pool and multiplying out each signature once.
///
/// The gradient is contracted per signature: the bin factors are summed over each signature's
/// groups first, then each signature's derivative weights are visited once.
///
/// @author Henry Wallace
class CpuReweight : public ReweightBackend {
 public:
  /// @param Model The handler's events, see ReweightModel
  /// @param Access How to shift, cut and bin one event
  CpuReweight(ReweightModel Model, ReweightEventAccess& Access);

  /// @copydoc ReweightBackend::Name
  std::string Name() const override { return "cpu"; }
  /// @copydoc ReweightBackend::Fill
  void Fill(const std::vector<double>& ShiftKey, bool FillW2, double* MC, double* W2) override;
  /// @copydoc ReweightBackend::AttachDerivatives
  void AttachDerivatives(const std::vector<WeightDerivativeSource>& Sources, int NPars) override;
  /// @copydoc ReweightBackend::Gradient
  void Gradient(const std::vector<double>& BinFactor, std::vector<double>& Out) override;
  /// @copydoc ReweightBackend::Invalidate
  void Invalidate() override { BinsValid = false; GroupsValid = false; }
  /// @copydoc ReweightBackend::EnableSurrogate
  bool EnableSurrogate(int Stride) override;
  /// @copydoc ReweightBackend::AnchorSurrogate
  void AnchorSurrogate() override;
  /// @copydoc ReweightBackend::SurrogateDelta
  void SurrogateDelta(bool FillW2, double* DMC, double* DW2) override;
  /// @copydoc ReweightBackend::Summary
  std::string Summary() const override;
  /// @copydoc ReweightBackend::Bytes
  std::size_t Bytes() const override;

 protected:
  /// @brief Whether the shifts are where the last full pass left them, so a cached pass will do
  bool BinsCurrent(const std::vector<double>& ShiftKey) const { return BinsValid && ShiftKey == BinKey; }
  /// @brief Read every pool entry and multiply out every signature
  void ComputeSignatureProducts();
  /// @brief Shift, cut and bin every live event, and sum the weights
  void FullPass(bool FillW2, double* MC, double* W2);
  /// @brief Sum the groups of every bin
  void CachedPass(bool FillW2, double* MC, double* W2);
  /// @brief Build the groups from the current bins, if they aren't already
  void EnsureGroups();
  /// @brief Make LiveBin current before the groups are built from it; a backend whose full pass
  /// leaves the bins elsewhere (on a GPU) fetches them here
  virtual void SyncLiveBin() {}
  /// @brief The groups indexed the other way round, for a GPU backend: each group's bin, and each
  /// signature's groups (signature s's are SigGroup[SigGroupStart[s]] onwards)
  void IndexGroups(std::vector<int>& GroupBin, std::vector<int>& SigGroupStart, std::vector<int>& SigGroup) const;
  /// @brief The gradient's terms grouped by parameter, for a GPU backend: parameter p's terms are
  /// ParStart[p] onwards, each a signature, a derivative (index into PoolDerivative) and a pool entry
  void ParameterTerms(std::vector<int>& ParStart, std::vector<int>& TermSig, std::vector<int>& TermDeriv,
                      std::vector<int>& TermPool) const;

  /// See the constructor
  ReweightModel Model;
  /// See the constructor
  ReweightEventAccess& Access;

  /// Current value of every pool entry
  std::vector<double> PoolValue;
  /// Current product of every signature's pool entries
  std::vector<double> SigProduct;

  /// Bin of every live event at the last full pass, negative if it contributes to none
  std::vector<int> LiveBin;
  /// Shift values the last full pass was made at
  std::vector<double> BinKey;
  /// Whether LiveBin belongs to BinKey
  bool BinsValid = false;
  /// The groups of LiveBin
  EventGroups Groups;
  /// Whether Groups was built from the current LiveBin
  bool GroupsValid = false;
  /// Bumped every time Groups is rebuilt, so a copy of them elsewhere knows it is out of date
  unsigned int GroupsVersion = 0;

  /// @name Gradient, set by AttachDerivatives
  /// @{
  /// Number of parameters
  int NPars = 0;
  /// Pool entry k's derivatives are PoolDerivParam / PoolDerivative[PoolDerivStart[k]] to [PoolDerivStart[k+1]-1]
  std::vector<int> PoolDerivStart;
  /// The parameter of each derivative
  std::vector<int> PoolDerivParam;
  /// Each derivative's address, nullptr for a derivative of 1
  std::vector<const M3::float_t*> PoolDerivative;
  /// Signature s's pool entries with a derivative are SigDerivPool[SigDerivStart[s]] to
  /// SigDerivPool[SigDerivStart[s+1]-1]
  std::vector<int> SigDerivStart;
  /// See SigDerivStart
  std::vector<int> SigDerivPool;
  /// Each group's bin, and signature s's groups GradSigGroup[GradSigGroupStart[s]] onwards (IndexGroups),
  /// for the gradient; built from the groups of version IndexedGroupsVersion
  std::vector<int> GradGroupBin;
  /// See GradGroupBin
  std::vector<int> GradSigGroupStart;
  /// See GradGroupBin
  std::vector<int> GradSigGroup;
  /// See GradGroupBin
  unsigned int IndexedGroupsVersion = 0;
  /// @}

  /// @name Surrogate, set by EnableSurrogate and AnchorSurrogate
  /// @{
  /// Every SurStride-th live event is sampled
  int SurStride = 0;
  /// The sampled live events
  std::vector<int> SurLive;
  /// Their bins at the anchor
  std::vector<int> SurBin;
  /// Their weights at the anchor
  std::vector<double> SurWeight;
  /// @}
};
