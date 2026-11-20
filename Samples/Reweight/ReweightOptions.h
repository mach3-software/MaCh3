#pragma once

#include "Manager/Manager.h"
#include "Samples/Reweight/ReweightBackend.h"
#include "Samples/Reweight/ReweightModel.h"

#include <memory>
#include <string>

/// @brief A sample config's Reweight section, read and checked
///
/// @details
/// ```yaml
/// Reweight:
///   Backend: auto        # auto | legacy | cpu | cuda | metal
///   SelfCheckPoints: 3   # points compared with the plain event loop at start-up, 0 to skip
///   FastAccess: true     # read shifted cuts and bins through pointers
///   LinearShifts: false  # work out linear shifts as arithmetic
/// ```
/// See Doc/GradientSampling.md.
///
/// @author Henry Wallace
struct ReweightOptions {
  /// @brief Read the section, which may be missing, and throw on a backend that can't be used
  /// @param Config The sample config
  /// @param SampleName For the messages
  static ReweightOptions Read(const YAML::Node& Config, const std::string& SampleName);

  /// Which backend to use
  std::string Backend = "auto";
  /// Points compared with the plain event loop at start-up
  int SelfCheckPoints = 3;
  /// Whether to read shifted cuts and bins through pointers, see FastEventAccess
  bool FastAccess = true;
  /// Whether to work out linear shifts as arithmetic, see LinearShiftTable
  bool LinearShifts = false;
};

/// @brief Build the backend a choice asks for
/// @param Choice "auto", "cpu", "cuda" or "metal" (not "legacy", which has no backend)
/// @param Model What every backend is built from
/// @param Access How the backend shifts, cuts and bins events
/// @details "auto" takes the CUDA backend in a CUDA build with a GPU present, otherwise the CPU one
std::unique_ptr<ReweightBackend> MakeReweightBackend(const std::string& Choice, ReweightModel Model, ReweightEventAccess& Access);
