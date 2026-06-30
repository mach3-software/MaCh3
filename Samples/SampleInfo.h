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

  /// In case experiment specific would like to have pointer to channel after using `GetOscChannel`, they can set pointer to this
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

/// @brief KS: Store info about MC sample
struct SampleInfo {
  /// Default constructor
  SampleInfo() = default;

  /// Destructor
  ~SampleInfo() {
    if(DataHist != nullptr) delete DataHist;
    if(MCHist   != nullptr) delete MCHist;
    if(W2Hist   != nullptr) delete W2Hist;
  }

  /// the strings associated with the variables used for the binning e.g. "RecoNeutrinoEnergy"
  std::vector<std::string> VarStr;

  /// @brief the name of this sample e.g."muon-like" used for printing
  std::string SampleTitle = "";
  /// @brief tag for sample used to easily set by which uncertainties should be affected
  std::string SampleName = "";

  /// @brief Keep track of the dimensions of the sample binning
  int nDimensions = M3::_BAD_INT_;

  /// names of mc files associated associated with this object
  std::vector<std::vector<std::string>> mc_files;
  /// names of spline files associated associated with this object
  std::vector<std::string> spline_files;

  /// Stores info about oscillation channel for a single sample
  std::vector<OscChannelInfo> OscChannels;

  /// histogram used for plotting storing data distribution
  TH1 *DataHist = nullptr;
  /// histogram used for plotting storing MC distribution
  TH1* MCHist = nullptr;
  /// histogram used for plotting storing W2 distribution
  TH1* W2Hist = nullptr;
};
