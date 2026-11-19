#pragma once

#include <memory>
#include <vector>

#include "Samples/OscillationProbability.h"

/// @brief OscillationGradient's derivatives worked out on an NVIDIA GPU
///
/// @details The same formulae (OscillationProbability.h) as the CPU, with one thread per weight and
/// free parameter. Each thread carries a single tangent (M3::Osc::OscDual1), which is bit for bit the
/// matching component of the eight-tangent calculation, so the dual numbers stay in registers. The
/// constants are worked out once on the host and sliced into one set per free parameter.
///
/// Plain C++ with standard types, so the CUDA code (OscillationGradientDevice.cu) never compiles
/// MaCh3's headers. Only built with MaCh3_GPU_ENABLED.
///
/// @author Henry Wallace
class OscillationGradientDevice {
 public:
  /// @brief Upload what each weight is the probability of
  /// @param Energy,Anti,From,To One entry per weight, as OscillationGradient keeps them
  OscillationGradientDevice(const std::vector<double>& Energy, const std::vector<char>& Anti,
                            const std::vector<signed char>& From, const std::vector<signed char>& To);
  ~OscillationGradientDevice();
  OscillationGradientDevice(const OscillationGradientDevice&) = delete;
  OscillationGradientDevice& operator=(const OscillationGradientDevice&) = delete;

  /// @brief Every weight's derivative with respect to each free parameter
  /// @param C The constants at the current parameters, as eight-tangent duals
  /// @param FreeSlot Which of the eight tangents to work out
  /// @param Out Derivative j of weight w at Out[w * FreeSlot.size() + j]
  void Evaluate(const M3::Osc::Constants<M3::Osc::OscDual>& C, const std::vector<int>& FreeSlot, double* Out);

  /// @brief Whether a CUDA device is available
  static bool Available();

 private:
  struct Impl;
  std::unique_ptr<Impl> D;
};
