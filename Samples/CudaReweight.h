#pragma once

#include <memory>

#include "Samples/CpuReweight.h"

class CudaReweightDevice;

/// @brief The reweight on an NVIDIA GPU, in double precision
///
/// @details A CpuReweight whose cached pass (the one every leapfrog makes, with the shifts where
/// they were) and gradient run on the GPU: multiplying out the signatures, summing each bin's
/// groups, and contracting the gradient, with the same arithmetic as the CPU, so the two agree to
/// rounding.
///
/// Full passes, where the shifts moved (every random walk sweep), run on the GPU too when the
/// shifts' effect is known in closed form: linear shifts (Reweight.LinearShifts) holding for every
/// live event, cut pointers and uniform one dimensional binning (see SweepTables). Each event's
/// variables then come out bit for bit as on the CPU, so it lands in the same bin. The bins stay on
/// the GPU and are fetched only when a cached pass needs groups built from them. Otherwise full
/// passes stay on the CPU, through the experiment's own shift and cut code.
///
/// Each pass reads the pool (and, for the gradient, the derivatives) on the CPU into pinned memory,
/// copies it over once, runs the kernels on one stream, and copies the bins (or the gradient) back.
/// The kernels are in CudaReweightDevice.cu; this class is plain C++.
///
/// @author Henry Wallace
class CudaReweight : public CpuReweight {
 public:
  /// @copydoc CpuReweight::CpuReweight
  CudaReweight(ReweightModel Model, ReweightEventAccess& Access);
  ~CudaReweight() override;

  /// @copydoc ReweightBackend::Name
  std::string Name() const override { return "cuda"; }
  /// @copydoc ReweightBackend::Fill
  void Fill(const std::vector<double>& ShiftKey, bool FillW2, double* MC, double* W2) override;
  /// @copydoc ReweightBackend::AttachDerivatives
  void AttachDerivatives(const std::vector<WeightDerivativeSource>& Sources, int NPars) override;
  /// @copydoc ReweightBackend::Gradient
  std::vector<double> Gradient(const std::vector<double>& BinFactor) override;
  /// @copydoc ReweightBackend::Summary
  std::string Summary() const override;
  /// @copydoc ReweightBackend::Bytes
  std::size_t Bytes() const override;
  /// @copydoc ReweightBackend::Invalidate
  void Invalidate() override { CpuReweight::Invalidate(); BinsOnDevice = false; }

  /// @brief Whether a CUDA device is available to run on
  static bool Available();

 private:
  /// The device side
  std::unique_ptr<CudaReweightDevice> GPU;
  /// The GroupsVersion last uploaded
  unsigned int UploadedGroups = 0;
  /// Whether the last Fill() ran on the GPU, so its gradient has to as well
  bool LastFillOnGPU = false;
  /// The fast access, if full passes can run on the GPU (it gives the shifts' offsets)
  const class FastEventAccess* Sweep = nullptr;
  /// Whether the GPU holds bins from a full pass that LiveBin doesn't have yet
  bool BinsOnDevice = false;
  /// Scratch for the shifts' offsets
  std::vector<double> Delta;
  /// @brief Copy the GPU's bins into LiveBin, if it has newer ones
  void SyncLiveBin() override;
  /// @brief Read the pool into the GPU's pinned buffer
  void ReadPool();
  /// Scratch for the bins the GPU returns
  std::vector<double> DeviceMC;
  /// See DeviceMC
  std::vector<double> DeviceW2;
};
