#pragma once

#include <memory>

#include "Samples/Reweight/CpuReweight.h"

/// @brief The reweight on an Apple GPU, in single precision, for testing on Macs
///
/// @details A CpuReweight whose cached pass (the one every leapfrog makes, with the shifts where
/// they were) runs on the GPU through Metal: multiplying out the signatures, summing each bin's
/// groups, and contracting the gradient. Full passes, where the shifts moved, stay on the CPU, since
/// they run the experiment's own shift and cut code.
///
/// Apple GPUs have no double precision, so the sums are in float and agree with the CPU to about
/// 1e-6 relative rather than to rounding. That is fine for testing, not for production, so this
/// backend is only used when asked for (Reweight: Backend: metal), never by "auto".
///
/// Every kernel is arranged so each output has exactly one writer: no atomics, so the results are
/// the same from run to run. The GPU and CPU share memory on Apple silicon, so the inputs are
/// written straight into the GPU's buffers.
///
/// The Metal code lives in MetalReweight.mm, which only builds on Apple with MaCh3_METAL_ENABLED.
///
/// @author Henry Wallace
class MetalReweight : public CpuReweight {
 public:
  /// @copydoc CpuReweight::CpuReweight
  MetalReweight(ReweightModel Model, ReweightEventAccess& Access);
  ~MetalReweight() override;

  /// @copydoc ReweightBackend::Name
  std::string Name() const override { return "metal"; }
  /// @copydoc ReweightBackend::Fill
  void Fill(const std::vector<double>& ShiftKey, bool FillW2, double* MC, double* W2) override;
  /// @copydoc ReweightBackend::AttachDerivatives
  void AttachDerivatives(const std::vector<WeightDerivativeSource>& Sources, int NPars) override;
  /// @copydoc ReweightBackend::Gradient
  void Gradient(const std::vector<double>& BinFactor, std::vector<double>& Out) override;
  /// @copydoc ReweightBackend::Summary
  std::string Summary() const override;
  /// @copydoc ReweightBackend::Precision
  double Precision() const override { return 1e-5; }

  /// @brief Whether a Metal device is available to run on
  static bool Available();

 private:
  /// Everything Metal, kept out of this header so it compiles as plain C++
  struct Device;
  /// See Device
  std::unique_ptr<Device> GPU;
  /// Whether the last Fill() ran on the GPU, so its gradient has to as well
  bool LastFillOnGPU = false;
};
