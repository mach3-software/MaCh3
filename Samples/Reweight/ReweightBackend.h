#pragma once

#include <string>
#include <vector>

#include "Splines/WeightDerivativeSource.h"
#include "Manager/Core.h"

/// @brief What a reweight backend needs to ask the sample handler about one event
///
/// @details The shifts, cuts and binning stay the experiment's own code. A backend only calls them,
/// and only when the shifts have moved; between those calls it reuses what they said.
class ReweightEventAccess {
 public:
  virtual ~ReweightEventAccess() = default;
  /// @brief Apply the functional shifts to an event, at the current parameter values
  virtual void ShiftEvent(int iEvent) = 0;
  /// @brief Whether a shifted event passes its sample's stored cuts
  virtual bool PassesCuts(int iEvent) = 0;
  /// @brief Global bin of a shifted event, M3::UnderOverFlowBin if it is outside the binning
  virtual int GlobalBinOf(int iEvent) = 0;
  /// @brief Shift an event, apply its cuts and find its bin, in one call
  /// @return Its global bin, or -1 if it is either not selected or outside the binning
  /// @details A backend's full pass calls only this, so an access can do the three together, or
  /// without touching the event at all (see FastEventAccess's linear shifts)
  virtual int ShiftCutAndBin(const int iEvent) {
    ShiftEvent(iEvent);
    if (!PassesCuts(iEvent)) return -1;
    const int Bin = GlobalBinOf(iEvent);
    return (Bin > M3::UnderOverFlowBin) ? Bin : -1;
  }
};

/// @brief A way of computing a sample handler's MC prediction and its likelihood gradient
///
/// @details Strategy for SampleHandlerBase::FillArray. Every backend gives the same sums as the
/// plain event loop, to rounding: SampleHandlerBase::CheckReweightBackend compares them.
///
/// Fill() is told the current values of the functional (shift) parameters. When those haven't
/// moved since the last full pass, every event's selection and bin is what that pass found, and a
/// backend may reuse them.
///
/// @author Henry Wallace
class ReweightBackend {
 public:
  virtual ~ReweightBackend() = default;

  /// @brief Short name for the log and the config, e.g. "cpu"
  virtual std::string Name() const = 0;

  /// @brief Fill the MC prediction (and sum of squared weights) at the current parameters
  /// @param ShiftKey The current values of every functional parameter
  /// @param FillW2 Whether to fill W2 as well
  /// @param MC NBins entries to add the prediction to (already zeroed)
  /// @param W2 NBins entries to add the squared weights to, only touched with FillW2
  virtual void Fill(const std::vector<double>& ShiftKey, bool FillW2, double* MC, double* W2) = 0;

  /// @brief Prepare to compute gradients
  /// @param Sources Where the derivatives of the event weights live
  /// @param NPars Number of parameters the gradient is indexed by
  virtual void AttachDerivatives(const std::vector<WeightDerivativeSource>& Sources, int NPars) = 0;

  /// @brief \f$\partial(-\log\mathcal{L})/\partial x_p\f$ at the last Fill()
  /// @param BinFactor \f$\partial(-\log\mathcal{L})/\partial\mu_b\f$ for every bin
  /// @param Out Resized to one entry per parameter and overwritten
  /// @details Takes a buffer rather than returning one: a sampler asks for this every leapfrog step
  virtual void Gradient(const std::vector<double>& BinFactor, std::vector<double>& Out) = 0;

  /// @brief Forget everything reused between fills, so the next one starts from scratch
  virtual void Invalidate() = 0;

  /// @name Surrogate for moves of the shifts alone, see SampleHandlerBase::ReweightSurrogate
  /// @{
  /// @brief Use every Stride-th live event for the surrogate's bin migrations
  /// @return Whether this backend can give the surrogate at all
  virtual bool EnableSurrogate([[maybe_unused]] const int Stride) { return false; }
  /// @brief Record the sampled events' bins and weights at the last Fill(), which is the anchor
  virtual void AnchorSurrogate() {}
  /// @brief At the current shifts, the sampled events' migrations from their anchor bins, times
  /// Stride: added to DMC (and their squared weights to DW2, if FillW2), both zeroed by the caller
  /// @details Only the shifts may have moved since the anchor, so every event keeps its weight
  virtual void SurrogateDelta([[maybe_unused]] const bool FillW2, [[maybe_unused]] double* DMC,
                              [[maybe_unused]] double* DW2) {}
  /// @}

  /// @brief A one line description for the log
  virtual std::string Summary() const = 0;
  /// @brief Memory the backend holds, in bytes
  virtual std::size_t Bytes() const = 0;
  /// @brief How closely this backend reproduces the plain event loop, relative
  /// @details Rounding for a double precision backend; a single precision one says so here, and is
  /// checked against that instead
  virtual double Precision() const { return 1e-9; }
};
