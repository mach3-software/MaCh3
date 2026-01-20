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
    if(dathist   != nullptr) delete dathist;
    if(dathist2d != nullptr) delete dathist2d;
    if(_hPDF1D   != nullptr) delete _hPDF1D;
    if(_hPDF2D   != nullptr) delete _hPDF2D;
  }

  /// the strings associated with the variables used for the binning e.g. "RecoNeutrinoEnergy"
  std::vector<std::string> VarStr;

  /// @brief the name of this sample e.g."muon-like"
  std::string SampleTitle = "";

  /// @brief Keep track of the dimensions of the sample binning
  int nDimensions = M3::_BAD_INT_;

  /// names of mc files associated associated with this object
  std::vector<std::string> mc_files;
  /// names of spline files associated associated with this object
  std::vector<std::string> spline_files;

  /// Stores info about oscillation channel for a single sample
  std::vector<OscChannelInfo> OscChannels;

  /// histogram used for plotting storing 1D data distribution
  TH1D *dathist = nullptr;
  /// histogram used for plotting storing 2D data distribution
  TH2D *dathist2d = nullptr;

  /// histogram used for plotting storing 1D MC distribution
  TH1D* _hPDF1D = nullptr;
  /// histogram used for plotting storing 2D MC distribution
  TH2D* _hPDF2D = nullptr;

  /// @brief Initialise histograms used for plotting
  void InitialiseHistograms() {
    TString histname1d = (VarStr[0]).c_str();
    TString histname2d = (VarStr[0] + "_" + VarStr[1]).c_str();
    TString histtitle = SampleTitle;

    //The binning here is arbitrary, now we get info from cfg so the
    //set1DBinning and set2Dbinning calls below will make the binning
    //to be what we actually want
    _hPDF1D   = new TH1D("h" + histname1d + SampleTitle, histtitle, 1, 0, 1);
    dathist   = new TH1D("d" + histname1d + SampleTitle, histtitle, 1, 0, 1);
    _hPDF2D   = new TH2D("h" + histname2d + SampleTitle, histtitle, 1, 0, 1, 1, 0, 1);
    dathist2d = new TH2D("d" + histname2d + SampleTitle, histtitle, 1, 0, 1, 1, 0, 1);

    _hPDF1D->SetDirectory(nullptr);
    dathist->SetDirectory(nullptr);
    _hPDF2D->SetDirectory(nullptr);
    dathist2d->SetDirectory(nullptr);

    // Set all titles so most of projections don't have empty titles...
    _hPDF1D->GetXaxis()->SetTitle(VarStr[0].c_str());
    _hPDF1D->GetYaxis()->SetTitle("Events");

    dathist->GetXaxis()->SetTitle(VarStr[0].c_str());
    dathist->GetYaxis()->SetTitle("Events");

    _hPDF2D->GetXaxis()->SetTitle(VarStr[0].c_str());
    _hPDF2D->GetYaxis()->SetTitle(VarStr[1].c_str());

    dathist2d->GetXaxis()->SetTitle(VarStr[0].c_str());
    dathist2d->GetYaxis()->SetTitle(VarStr[1].c_str());
  }
};

/// @brief constructors are same for all three so put in here
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski
struct EventInfo {
  /// @brief Default constructor.
  EventInfo(){}
  /// @brief Copy constructor (deleted to prevent copying).
  EventInfo(EventInfo const &other) = delete;
  /// @brief Move constructor (defaulted to allow moving).
  EventInfo(EventInfo &&other) = default;
  /// @brief Copy assignment operator (deleted).
  EventInfo& operator=(EventInfo const &other) = delete;
  /// @brief Move assignment operator (deleted).
  EventInfo& operator=(EventInfo &&other) = delete;
  /// @brief default destructor
  ~EventInfo(){}

  /// target the interaction was on
  const int* Target = 0;
  /// PDG of neutrino after oscillation
  const int* nupdg  = 0;
  /// PDG of neutrino before oscillation
  const int* nupdgUnosc = 0;

  /// Pointer to true Neutrino Energy
  const double* rw_etru = &M3::_BAD_DOUBLE_;
  /// Pointer to true cosine zenith
  const double* rw_truecz = &M3::_BAD_DOUBLE_;

  /// Pointers to normalisation weights which are being taken from Parameter Handler
  std::vector<const double*> xsec_norm_pointers;

  /// Pointers to weights like oscillation spline etc
  std::vector<const M3::float_t*> total_weight_pointers;

  /// The x_var and y_vars and beyond that you're binning in
  std::vector<const double*> KinVar;
  /// starting bins for each dimensions allowing to perform quick lookup
  std::vector<int> NomBin;

  /// Nominal sample to which event is associated
  int NominalSample = M3::_BAD_INT_;
  /// Is event NC or not
  bool isNC = false;

  /// Pointer to MaCh3 mode
  const double* mode = &M3::Unity_D;
};
