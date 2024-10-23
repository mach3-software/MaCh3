#pragma once

/// @brief constructors are same for all three so put in here
struct fdmc_base {
  int nutype; // 2 = numu/signue | -2 = numub | 1 = nue | -1 = nueb           
  int oscnutype;    
  int nupdg;
  int nupdgUnosc;
  bool signal; // true if signue                                              
  int nEvents; // how many MC events are there
  double ChannelIndex;
  std::string flavourName;

  int **Target; // target the interaction was on

  int SampleDetID;

  //THe x_var and y_vars that you're binning in
  const double** x_var;
  const double** y_var;
  const double **rw_etru;
  const double **rw_truecz = NULL;

  /// xsec bins
  std::list< int > *xsec_norms_bins;

  /// DB Speedup bits
  double Unity;
  float Unity_F;
  int Unity_Int;
  double dummy_value = -999;

  int* nxsec_norm_pointers;
  const double*** xsec_norm_pointers;

  int* nxsec_spline_pointers;
  const double*** xsec_spline_pointers;

  int* ntotal_weight_pointers;
  const double*** total_weight_pointers;
  double* total_w;

  int* XBin;
  int* YBin;
  int* NomXBin;
  int* NomYBin;

  bool *isNC;

  // histo pdf bins
  double *rw_lower_xbinedge; // lower to check if Eb has moved the erec bin
  double *rw_lower_lower_xbinedge; // lower to check if Eb has moved the erec bin
  double *rw_upper_xbinedge; // upper to check if Eb has moved the erec bin
  double *rw_upper_upper_xbinedge; // upper to check if Eb has moved the erec bin

  double **mode;

  const _float_ **osc_w_pointer;
  double *xsec_w;
  splineFDBase *splineFile; 
};
