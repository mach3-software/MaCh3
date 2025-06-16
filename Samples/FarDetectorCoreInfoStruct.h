#pragma once
#include <vector>
#include <string>

#include "Samples/SampleStructs.h"

/// @brief KS: Store info about used osc channels
struct OscChannelInfo {
  /// Name of osc channel
  std::string flavourName;
  /// Fancy channel name (e.g., LaTeX formatted)
  std::string flavourName_Latex;

  /// PDG of initial flavour
  int InitPDG;
  /// PDG of oscillated/final flavour
  int FinalPDG;

  /// In case experiment specific would like to have pointer to channel after using @GetOscChannel, they can set pointer to this
  double ChannelIndex;
};

/// @brief KS: Get Osc Channel Index based on initial and final PDG codes
/// @param OscChannel The vector of available oscillation channels
/// @param InitFlav Initial flavour PDG code
/// @param FinalFlav Final flavour PDG code
/// @return Index in OscChannel vector
inline int GetOscChannel(const std::vector<OscChannelInfo>& OscChannel, const int InitFlav, const int FinalFlav) {
  for (size_t i = 0; i < OscChannel.size(); ++i) {
    if (InitFlav == OscChannel[i].InitPDG && FinalFlav == OscChannel[i].FinalPDG) {
      return static_cast<int>(OscChannel[i].ChannelIndex);
    }
  }

  MACH3LOG_ERROR("Didn't find Osc channel for InitFlav = {}, FinalFlav = {}", InitFlav, FinalFlav);
  throw MaCh3Exception(__FILE__, __LINE__);
}

/// @brief constructors are same for all three so put in here
struct FarDetectorCoreInfo {
  FarDetectorCoreInfo(){}
  FarDetectorCoreInfo(FarDetectorCoreInfo const &other) = delete;
  FarDetectorCoreInfo(FarDetectorCoreInfo &&other) = default;
  FarDetectorCoreInfo& operator=(FarDetectorCoreInfo const &other) = delete;
  FarDetectorCoreInfo& operator=(FarDetectorCoreInfo &&other) = delete;

  ~FarDetectorCoreInfo(){}

  const int* Target = 0; ///< target the interaction was on
  const int* nupdg  = 0;
  const int* nupdgUnosc = 0;

  //THe x_var and y_vars that you're binning in
  const double* x_var = &M3::Unity_D;
  const double* y_var = &M3::Unity_D;
  const double* rw_etru = &M3::_BAD_DOUBLE_;
  const double* rw_truecz = &M3::_BAD_DOUBLE_;

  /// Pointers to normalisation weights which are being taken from Parameter Handler
  std::vector<const double*> xsec_norm_pointers;
  /// Pointers to spline weights which are being calculated by Splines Handler
  std::vector<const M3::float_t*> xsec_spline_pointers;
  /// Total weight of norm and spline parameters
  M3::float_t xsec_w = 1.;
  /// pointer to oscillation weight which is being calculated by Oscillation Handler
  const M3::float_t* osc_w_pointer = &M3::Unity;
  /// Pointers to @xsec_w, @osc_w_pointer and remaining experiment specific weights
  std::vector<const M3::float_t*> total_weight_pointers;

  //M3::float_t total_w  = M3::_BAD_INT_;

  int NomXBin = -1;
  int NomYBin = -1;

  bool isNC = false;

  const double* mode = &M3::Unity_D;
};
