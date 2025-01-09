#pragma once
#include <vector>
#include <string>
#include "samplePDF/Structs.h"

/// @brief constructors are same for all three so put in here
struct FarDetectorCoreInfo {
  FarDetectorCoreInfo() : isNC{nullptr} {}
  FarDetectorCoreInfo(FarDetectorCoreInfo const &other) = delete;
  FarDetectorCoreInfo(FarDetectorCoreInfo &&other) = default;
  FarDetectorCoreInfo& operator=(FarDetectorCoreInfo const &other) = delete;
  FarDetectorCoreInfo& operator=(FarDetectorCoreInfo &&other) = delete;

  ~FarDetectorCoreInfo(){ delete [] isNC; }

  bool signal; ///< true if signue
  int nEvents; ///< how many MC events are there
  double ChannelIndex;
  std::string flavourName;

  std::vector<int*> Target; ///< target the interaction was on
  std::vector<const int*> nupdg;
  std::vector<const int*> nupdgUnosc;

  //THe x_var and y_vars that you're binning in
  std::vector<const double*> x_var;
  std::vector<const double*> y_var;
  std::vector<const double*> rw_etru;
  std::vector<const double*> rw_truecz;

  /// xsec bins
  std::vector< std::vector< int > > xsec_norms_bins;

  std::vector<int> nxsec_norm_pointers;
  std::vector<std::vector<const double*>> xsec_norm_pointers;

  std::vector<int> nxsec_spline_pointers;
  std::vector<std::vector<const M3::float_t*>> xsec_spline_pointers;

  std::vector<int> ntotal_weight_pointers;
  std::vector<std::vector<const M3::float_t*>> total_weight_pointers;
  std::vector<M3::float_t> total_w;

  std::vector<int> XBin;
  std::vector<int> YBin;
  std::vector<int> NomXBin;
  std::vector<int> NomYBin;

  bool *isNC;

  // histo pdf bins
  std::vector<double> rw_lower_xbinedge;       ///< lower to check if Eb has moved the erec bin
  std::vector<double> rw_lower_lower_xbinedge; ///< lower to check if Eb has moved the erec bin
  std::vector<double> rw_upper_xbinedge;       ///< upper to check if Eb has moved the erec bin
  std::vector<double> rw_upper_upper_xbinedge; ///< upper to check if Eb has moved the erec bin

  std::vector<const double*> mode;

  std::vector<const M3::float_t*> osc_w_pointer;
  std::vector<M3::float_t> xsec_w;
};
