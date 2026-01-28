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

  /// histogram used for plotting storing data distribution
  TH1 *DataHist = nullptr;
  /// histogram used for plotting storing MC distribution
  TH1* MCHist = nullptr;
  /// histogram used for plotting storing W2 distribution
  TH1* W2Hist = nullptr;

  /// @brief Initialise histograms used for plotting
  void InitialiseHistograms() {
    std::string HistTitle = SampleTitle;

    //The binning here is arbitrary, now we get info from cfg so the
    //set1DBinning and set2Dbinning calls below will make the binning
    //to be what we actually want
    if(nDimensions == 1) {
      DataHist = new TH1D(("d" + HistTitle).c_str(), HistTitle.c_str(), 1, 0, 1);
      MCHist   = new TH1D(("h" + HistTitle).c_str(), HistTitle.c_str(), 1, 0, 1);
      W2Hist   = new TH1D(("w" + HistTitle).c_str(), HistTitle.c_str(), 1, 0, 1);

      // Set all titles so most of projections don't have empty titles...
      DataHist->GetXaxis()->SetTitle(VarStr[0].c_str());
      DataHist->GetYaxis()->SetTitle("Events");
      MCHist->GetXaxis()->SetTitle(VarStr[0].c_str());
      MCHist->GetYaxis()->SetTitle("Events");
      W2Hist->GetXaxis()->SetTitle(VarStr[0].c_str());
      W2Hist->GetYaxis()->SetTitle("Events");
    } else if(nDimensions == 2) {
      DataHist = new TH2D(("d" + HistTitle).c_str(), HistTitle.c_str(), 1, 0, 1, 1, 0, 1);
      MCHist   = new TH2D(("h" + HistTitle).c_str(), HistTitle.c_str(), 1, 0, 1, 1, 0, 1);
      W2Hist   = new TH2D(("w" + HistTitle).c_str(), HistTitle.c_str(), 1, 0, 1, 1, 0, 1);

      // Set all titles so most of projections don't have empty titles...
      DataHist->GetXaxis()->SetTitle(VarStr[0].c_str());
      DataHist->GetYaxis()->SetTitle(VarStr[1].c_str());
      MCHist->GetXaxis()->SetTitle(VarStr[0].c_str());
      MCHist->GetYaxis()->SetTitle(VarStr[1].c_str());
      W2Hist->GetXaxis()->SetTitle(VarStr[0].c_str());
      W2Hist->GetYaxis()->SetTitle(VarStr[1].c_str());
    } else {
      MACH3LOG_DEBUG("Not supported for Dim {}", nDimensions);
      return;
    }

    DataHist->SetDirectory(nullptr);
    MCHist->SetDirectory(nullptr);
    W2Hist->SetDirectory(nullptr);
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
  std::vector<const double*> norm_pointers;

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
