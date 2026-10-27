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
  BinnedSplineMakerBase(const std::string& Config);
  virtual ~BinnedSplineMakerBase();
  void MakeVariations();

  /// @brief ETA function to generically convert a string from param handler to a kinematic type
  int ReturnKinematicParameterFromString(const std::string& KinematicStr) const;
  /// @brief ETA function to generically convert a kinematic type from param handler to a string
  std::string ReturnStringFromKinematicParameter(const int KinematicVariable) const;

 protected:
  void SetupBinning();
  void WriteSplines();
  void SetBinning(const std::vector<double>& var0Bins,
                  const std::vector<double>& var1Bins,
                  const std::vector<double>& var2Bins);
  void Initialise();
  void LoadSystematics();
  void AddSystematics();

  void FillSystematicProperties();
  void FillGraphs();
  void FillSplines();

  virtual void LoadVariables() = 0;
  virtual void LoadKnots() = 0;
  virtual std::vector<double> GetETrueBinning() = 0;
  double ReturnKinematicParameter(const std::string& KinematicParameter, Long64_t iEvent) const {
    return ReturnKinematicParameter(ReturnKinematicParameterFromString(KinematicParameter), iEvent);
  }
  virtual double ReturnKinematicParameter(const int binTypes, const Long64_t jentry) const = 0;
  virtual Long64_t GetEntries()  = 0;
  virtual void GetEntry(Long64_t jentry) = 0;
  virtual int GetMode(Long64_t jentry) = 0;
  virtual double GetWeight(Long64_t jentry, size_t iSyst, int iKnot) = 0;
  virtual double GetWeight(Long64_t jentry, size_t iSyst, double Value) = 0;

  std::unique_ptr<Manager> Setting;
  std::unique_ptr<MaCh3Modes> Modes;

  // Vectors of splines, graphs, and histograms
  std::vector< std::vector< std::unique_ptr<TH3F> > > dev_tmp;
  std::vector< std::string > dev_tmp_names;
  std::vector< std::vector< std::vector< std::vector< std::unique_ptr<TSpline3> > > > > splines;
  std::vector< std::vector< std::vector< std::vector< std::unique_ptr<TGraph> > > > >   graphs;

  std::vector< int >   binTypes_vec;
  std::vector<SystematicsStruct> systematicProperties;

  /// Suffix -> spline mode index
  std::map<std::string, int> fSuffixToSplineMode;

  /// MaCh3 mode index -> spline mode index
  std::map<int, int> fModeToSplineMode;

  bool UseKnotForWeight;

  /// Mapping between string and kinematic enum
  const std::unordered_map<std::string, int>* KinematicParameters;
  /// Mapping between kinematic enum and string
  const std::unordered_map<int, std::string>* ReversedKinematicParameters;
};
