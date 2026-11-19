#pragma once

#include "Manager/Manager.h"
#include "Manager/MaCh3Modes.h"
#include "Splines/SplineStructs.h"
#include "Splines/SplineCommon.h"
#include "Samples/HistogramUtils.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT include
#include "TArrayF.h"
#include "TH3F.h"
_MaCh3_Safe_Include_End_ //}

/// @brief Struct to contain systematic properties
/// @author Henry Wallace
struct SystematicsStruct {
  /// name of this systematic
  std::string systName;

  std::string outputSystName;

  /// which interaction modes does this systematic affect?
  std::vector<int> intModes;

  /// pointers to the arrays we will read in from the weights file
  TArrayF* knotArray;

  /// which element # in the knot array corresponds to this systematics' nominal value?
  int nominalPosition;
};

/// @brief Base class for producing binned splines in MaCh3
class BinnedSplineMakerBase {
 public:
  /// @brief Constructor
  BinnedSplineMakerBase(const std::string& Config);
  /// @brief Destructor
  virtual ~BinnedSplineMakerBase();
  /// @brief The main function generating spline responses
  void MakeVariations();

 protected:
  /// @brief ETA function to generically convert a string from param handler to a kinematic type
  int ReturnKinematicParameterFromString(const std::string& KinematicStr) const;
  /// @brief ETA function to generically convert a kinematic type from param handler to a string
  std::string ReturnStringFromKinematicParameter(const int KinematicVariable) const;
  std::string GetSplineSuffix(const int spline_mode) const;
  int GetSplineModeFromGenerator(const int generatorMode) const;
  /// @brief The main initialiser, has to be called from inherited class :(
  void Initialise();
  /// @brief Hook allowing to load experiment spagyric variables trees etc.
  virtual void LoadVariables() = 0;
  /// @brief Hook allowing to load knots etc from experiment specific
  virtual void LoadKnots() = 0;
  /// @brief Way to get binning outside of config
  /// @warning will become deprecated
  virtual std::vector<double> GetETrueBinning() = 0;
  /// @brief simple getter for number of entries
  virtual Long64_t GetEntries() = 0;
  /// @brief Hook allowing to load info for a given entry
  virtual void GetEntry(Long64_t jentry) = 0;
  /// @brief Calculate weight for a given event using knot
  virtual double GetWeight(Long64_t jentry, size_t iSyst, int iKnot) = 0;
  /// @brief Calculate weight for a given event using value
  virtual double GetWeight(Long64_t jentry, size_t iSyst, double Value) = 0;
  /// @brief Return the value of an associated kinematic parameter for an event
  double ReturnKinematicParameter(const std::string& KinematicParameter, Long64_t iEvent) const {
    return ReturnKinematicParameter(ReturnKinematicParameterFromString(KinematicParameter), iEvent);
  }
  /// @brief Return the value of an associated kinematic parameter for an event
  virtual double ReturnKinematicParameter(const int binTypes, const Long64_t jentry) const = 0;

  /// pointer to config
  std::unique_ptr<Manager> Setting;
  std::vector<SystematicsStruct> systematicProperties;

  /// Mapping between string and kinematic enum
  const std::unordered_map<std::string, int>* KinematicParameters;
  /// Mapping between kinematic enum and string
  const std::unordered_map<int, std::string>* ReversedKinematicParameters;
 private:
  /// @brief Set binning and histograms for each spline
  void SetupBinning();
  void WriteSplines();
  void LoadSystematics();
  void AddSystematics();

  void FillSystematicProperties();
  /// @brief Fit a TGraph to each bin
  void FillGraphs();
  /// @brief Fit a spline to each TGraph
  void FillSplines();

  /// pointer to MaCh3 modes
  std::unique_ptr<MaCh3Modes> Modes;

  /// Histograms containing systematic variations at each knot.
  std::vector<std::vector<std::unique_ptr<TH3F>>> VariationHistogram;
  /// Names corresponding to the variation histograms.
  std::vector<std::string> VariationHistogramName;
  /// Cubic spline interpolators for each bin and systematic mode.
  std::vector<std::vector<std::vector<std::vector<std::unique_ptr<TSpline3>>>>> splines;
  /// Graphs storing systematic variation ratios at each knot for each bin.
  std::vector<std::vector<std::vector<std::vector<std::unique_ptr<TGraph>>>>> graphs;

  /// Enum for binning type to be used
  std::vector< int > BinningTypes;

  /// Suffix -> spline mode index
  std::map<std::string, int> fSuffixToSplineMode;

  /// MaCh3 mode index -> spline mode index
  std::map<int, int> fModeToSplineMode;

  /// Whether we use knot for weight calcuation or actual value
  bool UseKnotForWeight;
};
