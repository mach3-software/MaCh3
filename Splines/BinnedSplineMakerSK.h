#pragma once
#include "BinnedSplineMakerBase.h"

// ROOT INCLUDES - Because ROOT sucks...
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TArrayF.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TMath.h"

struct ConstantsStruct{
  const double VNuc         = 27.0;       //MeV
  const double protonMass   = 938.272013; //MeV
  const double neutronMass  = 939.565346; //MeV
  const double deltaMass    = 1232.0;     //MeV
  const double electronMass = 0.510998;   //MeV
  const double muonMass     = 105.65836;  //MeV
  const double pionMass     = 139.57;     //MeV
};

// Now our enums :
enum iclassName{
  k1Re     = 11,
  k1Rmu    = 13,
  k1Re1de  = 14,
  kMRe1de  = 15,
  kMRMu1de = 31,
  kMRMu2de = 32
};

//ETA - adding in new function to calculate erec
inline Double_t CalculateErec(iclassName iclass, double LepMom, double costheta)
{
  struct ConstantsStruct constants;
  //Get the true neutrino PDG
  bool iselike = false;
  bool iscc1pi = false;

  //Now set some bools based on iclass
  switch(iclass) {
    case(k1Re):
      iselike = true;
      iscc1pi = false;
      break;
    case(k1Rmu):
      iselike = false;
      iscc1pi = false;
      break;
    case(k1Re1de):
      iselike = true;
      iscc1pi = true;
      break;
    case(kMRe1de):
      iselike = true;
      iscc1pi = true;
      break;
    case(kMRMu1de):
      iselike = false;
      iscc1pi = true;
      break;
    case(kMRMu2de):
      iselike = false;
      iscc1pi = true;
      break;
    default:
      std::cerr << "[ERR]: Unknown iclass so I don't know how to calculate erec :'( " << std::endl;
      std::cerr << __FILE__ << ":"<<__LINE__<<std::endl;
      throw;
      break;
  }

  Double_t erec;
  Double_t shiftedmom = LepMom;
  double mass;

  if(iselike) mass = constants.electronMass;
  else if(iscc1pi) mass = constants.pionMass; //ETA- there's a subtlety here. Muon momentum taken from a 2Ring PiPi fit so use pion mass to take this into account
  else        mass = constants.muonMass;

  double LeptonE = sqrt(mass*mass +shiftedmom*shiftedmom);
  erec = 0.;

  //!!calculate erec!!!
  if(!iscc1pi) {
    erec = (constants.neutronMass - constants.VNuc)*LeptonE - mass*mass*0.5;
    erec += constants.neutronMass*constants.VNuc - constants.VNuc*constants.VNuc*0.5;
    erec += (constants.protonMass*constants.protonMass - constants.neutronMass*constants.neutronMass)*0.5;
    erec /= (constants.neutronMass - constants.VNuc - LeptonE + shiftedmom*costheta);
    //End: Lifted from T2K-SK
  } else{
    erec = (constants.protonMass*LeptonE-mass*mass/2.);
    erec += (constants.deltaMass*constants.deltaMass - constants.protonMass*constants.protonMass)/2.;
    erec /= (constants.protonMass - LeptonE + shiftedmom*costheta);
    //End: Lifted from T2K-SK
  }

  return erec;
}

class BinnedSplineMakerSK : public BinnedSplineMakerBase {
 public:
  BinnedSplineMakerSK(const std::string& Config);
  virtual ~BinnedSplineMakerSK();

 protected:
   void LoadVariables() final;
   void LoadKnots() final;
   Int_t GetParticleID(int iEntry){
     fChain_sk->GetEntry(iEntry);
     return ipnu[0];
   }

  std::vector<double> GetETrueBinning() final;
  double GetVariableValue(int binTypes) final;

  Long64_t GetEntries() final {return fChain_sk->GetEntries();}
  void GetEntry(Long64_t jentry) final {fChain_sk->GetEntry(jentry);}
  int GetMode([[maybe_unused]] Long64_t jentry) final {return mode;}
  double GetEnu([[maybe_unused]] Long64_t jentry) final {return pnu[0];}
  double GetWeight(Long64_t jentry, size_t iSyst, int iKnot) final;
  int GetBinningType(const std::string& binName);
  std::vector<double> ExtractBinningFromTemplate(std::string template_file_name,
                                                               std::string  binning_name, bool atm);
  // Input file containing weights, properly formatted
  std::vector< TArrayF* > weightArray;

  // All of our non-mtuple variables
  // Pointer to the analyzed TTree or TChain
  TChain *fChain_sk;

  bool erec_mev;
  Int_t fCurrent_sk;
  bool is_atm;

  // Declaration of leaf types: fChain_sk
  Int_t           iclass;
  Int_t           Ibound;
  Double_t        wgtflx;
  Double_t        wgtosc;
  Double_t        E_rec;
  Float_t         fq1rdir[10][10][3];
  Float_t         fqmome;
  Float_t         fqmomm;
  //ETA adding in fiqun multi-ring variables needed for numuCC1pi sample(s)
  Float_t         fqmrdir[200][6][3];
  Float_t         fqmrmom[200][6];
  Float_t         fqmreloss[200][6];
  Int_t           ipp; //index of 2R pi+pi+ fit
  Int_t           iep; //index of 2R epi+ fit
  Int_t           mode;
  Int_t           ipnu[50];
  Float_t         pnu[50];

  // Atm Variables
  Float_t         genmom;
  Float_t         gencz;
  Float_t         solarweight;
  Float_t         oscweight3f;
  Double_t        FFF;

  //MC normalisation weight
  Double_t        wgtnvv;

  struct ConstantsStruct constants;

  enum binningType{kNoBinning, kErec, kTheta, kMom, kGenmom, kGencz};
};
