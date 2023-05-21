#ifndef __OSCCLASS_CUDAPROB3_H__
#define __OSCCLASS_CUDAPROB3_H__

#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TStopwatch.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>

#include "manager/manager.h"

//Propagator includes
#ifdef USE_PROB3
  #include "BargerPropagator.h"
#else
  #include "beamcudapropagator.cuh"
  #include "atmoscudapropagator.cuh"
#endif

using FLOAT_T = double;
//using FLOAT_T = float;
using namespace cudaprob3; // namespace of the propagators

struct Corner{
  double x;
  double y;
};

struct Box{
  Corner BL;
  Corner UR;
};

class Oscillator {
 public:
  Oscillator(std::string ConfigName);

  int GetOscillogramNBins(int Switcher);

  void FillOscillogram(double* oscpar, double prodH, double Yp_Val=0.468);
  double ReturnProb(double NeutrinoEnergy, double Cosz, int InitialFlavour, int FinalFlavour);

  void RebinOscillogram(int Switcher, std::vector<double> NewBinning);

  //DB Setters
  void SetFillHistograms(bool fFillHistograms_=true) {fFillHistograms = fFillHistograms_; if (fFillHistograms && !HistogramsInitialised) {InitialiseHistograms();}}
  void SetOscillatorConfig(int ArrayConfig_=0);

  //DB Getters
  std::vector<double> ReturnSecondaryXBinEdges() {return SecondaryXBinEdges;}
  std::vector<double> ReturnSecondaryYBinEdges() {return SecondaryYBinEdges;}

  std::vector< std::vector< std::vector<TH2D*> > > ReturnPrimaryOscillogram() {return hPrimaryOscillogram;}
  std::vector< std::vector< std::vector<TH2D*> > > ReturnSecondaryOscillogram() {return hSecondaryOscillogram;}
  std::vector<TH2D*> ReturnOscillogramArray(int fPrimary);

  bool isLinear() {return IsLinear;}
  TH2D* ReturnPrimaryBinning() {return hPrimaryBinning;}

  const double* retPointer(int GenNeutrinoFlavour, int DetNeutrinoFlavour, double NeutrinoEnergy, double TrueCZ);

  std::vector<FLOAT_T> ReturnProductionHeightBinEdges() {return ProductionHeightBinEdges;}

  //DB Console output
  void PrintBinning();
  void PrintOscillatorConfig();

  void SaveOscillogramsToFile(TString FileName);

 private:
  //########################################################################################################################
  //Functions

  void isUsingGPU();

  void SetProductionHeightArray();
  void SetProductionHeightBinEdges();
  void DefineMiscValues();

  void ResizeArrays();
  void ResizeOscillogram(bool CoarseChanged, bool FineChanged);

  void DeletePropagator();
  void DeleteArrays();
  void DeleteOscillogram();

  void FillArrays();
  void FillArrays_ManyContrib_Area();
  void FillArrays_Standard();

  void InitOscillogram(TH2D* PrimaryHistTemplate,TH2D* SecondaryHistTemplate);
  void InitPropagator();
  void InitialiseHistograms();
  void DeleteHistograms();

  void Reset(int NeutrinoSignIndex, int InitialNeutrinoIndex, int FinalNeutrinoIndex);

  void FillPrimaryOscillogram(int NeutrinoSignIndex, int InitialNeutrinoIndex, int FinalNeutrinoIndex);
  void FillSecondaryHistograms(int NeutrinoSignIndex, int InitialNeutrinoIndex, int FinalNeutrinoIndex);
  void FillPrimaryHistograms(int NeutrinoSignIndex, int InitialNeutrinoIndex, int FinalNeutrinoIndex);

  bool isAlreadyCalculated(double* oscpar, double prodH, double Yp_Val);
  void SaveParams(double* oscpar, double prodH, double Yp_Val);
  void ResetSavedParams();

  void PrintBox(Box Box1);
  bool IsValidBox(Box Box1);
  double FractionOverlapped(Box Box1, Box Box2);

  inline int NeutrinoIndexToSign(int val) {
    return 2*val-1;
  }

  inline int NeutrinoIndexToFlavour(int val) {
    return val+1;
  }

  inline int NeutrinoSignToIndex(int val) {
    return (val<0) ? 0 : 1;
  }

  inline int NeutrinoFlavourToIndex(int val) {
    return (abs(val)-1);
  }

  void DefinePropagatorEnums();
  void CheckEarthDensityFile();

  std::vector<double> linspace(double Emin, double Emax, int nDiv);
  std::vector<double> logspace(double Emin, double Emax, int nDiv);

  std::vector<double> ReturnFineBinningFromCoarseBinnnig(int nFine, std::vector<double> CoarseBinning);

  //########################################################################################################################
  //Variables

  TString EarthDensityFile;
  TString ProductionHeightFileName;

  std::unique_ptr<AtmosCpuPropagator<FLOAT_T>> propagator;

  unsigned int nPrimaryHists;
  double** hPrimaryOscillogram_Arr;

  double* hPrimaryCounter_Arr;

  int nNeutrinoSigns;
  int nInitialFlavours;
  int nFinalFlavours;

  FLOAT_T* ProbList;

  std::vector< std::vector<int> > PrimaryBinContrib_Bin;
  std::vector< std::vector<double> > PrimaryBinContrib_Weight;

  std::vector<double> PrimaryXBinEdges;
  std::vector<double> PrimaryYBinEdges;

  std::vector<double> SecondaryXBinEdges;
  std::vector<double> SecondaryYBinEdges;

  std::vector<FLOAT_T> PrimaryXBinEvalPoints;
  std::vector<FLOAT_T> PrimaryYBinEvalPoints;

  std::vector<FLOAT_T> SecondaryXBinEvalPoints;
  std::vector<FLOAT_T> SecondaryYBinEvalPoints;

  std::vector<NeutrinoType> NeutrinoTypes;
  std::vector<TString> NeutrinoTypes_Names;

  std::vector< std::vector<TString> > OscChannels_Names;
  std::vector< std::vector<ProbType> > OscChannels;

  int nSecondaryBins;
  int nPrimaryBins;

  // Number of layers from propagator
  int nLayers;

  bool fFillHistograms;
  bool HistogramsInitialised;
  int ArrayConfig;

  int nPrimaryBinsX;
  int nPrimaryBinsY;

  int nSecondaryBinsX;
  int nSecondaryBinsY;

  double lCoarseEnergy;
  double hCoarseEnergy;
  int nCoarseEnergy;

  double lFineEnergy;
  double hFineEnergy;
  int nFineEnergy;

  double lCoarseCosz;
  double hCoarseCosz;
  int nCoarseCosz;

  double lFineCosz;
  double hFineCosz;
  int nFineCosz;

  std::vector< std::vector< std::vector<TH2D*> > > hPrimaryOscillogram;
  std::vector< std::vector< std::vector<TH2D*> > > hSecondaryOscillogram;
  TH2D* hPrimaryBinning;
  TH2D* hSecondaryBinning;

  double fprodH;
  double fYp_Val;
  double foscpar[6];
  double fBinning[4];

  int nOscpars;

  bool IsLinear;
  bool useFineBinsPerBin;
  bool RebinMode;
  int nMaxBin;

  std::string InputFileName;
  TString PrimaryHistKey;
  TString SecondaryHistKey;

  bool UseBinningTemplates;
  std::string PrimaryBinningTemplateName;
  std::string SecondaryBinningTemplateName;
  TString TemplateInputFileName;

  bool UseChemicalComposition;

  double fineCoarseRatioCosz;
  double fineCoarseRatioEnergy;

  bool UseProductionHeightAveraging;
  int nProductionHeightAveragingBins;
  double lProductionHeightRange;
  double hProductionHeightRange;

  std::vector<FLOAT_T> ProductionHeightBinEdges;
  std::vector<FLOAT_T> chemicalComposition_Nom;
};

#endif
