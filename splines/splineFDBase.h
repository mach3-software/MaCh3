#ifndef _splineFDBase_h_
#define _splineFDBase_h_

#ifndef __BAD_SPLINE__
#define __BAD_SPLINE__ 123456789
#endif

//ETA adding the option to use reduced TSpline3. Hopefully this aves us some memory!
// Do we use TF1 or TSpline3* for spline evaluations
#define USE_TSpline3_FD 1
#define USE_TSpline3_red_FD 2

// Can use:
//  TSpline3 (third order spline in ROOT)
//  TSpline3_red (reduced class third order spline in ROOT)

#define USE_SPLINE_FD USE_TSpline3_red_FD
// Set the __SPLINE_TYPE__ accordingly
#if USE_SPLINE_FD == USE_TSpline3_FD
#define __SPLINE_TYPE_FD__ TSpline3_FD
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
#define __SPLINE_TYPE_FD__ TSpline3_red_FD
#endif

#include "splineBase.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TKey.h"
#include "TROOT.h"
#include "TClass.h"
//ETA - need Structs.h because this is where TSpline3 and FastSplineEval are defined but maybe we should move them to a new header?
#include "../samplePDF/Structs.h"
//ETA - covariance xsec need to get the number of splines which apply to each det for example
#include "../covariance/covarianceXsec.h"

// Note: etrue-var1-var2 spline binning not yet implemented for 2017!

class splineFDBase : public splineBase
{
 public:
  splineFDBase(const char *spline, int nutype, int nevents, int DetID, covarianceXsec* xsec_cov = NULL); // constructor for etrue-var1 splines
  splineFDBase(const char *spline, int nutype, int nevents, double BinningOpt, int DetID, covarianceXsec* xsec_cov = NULL); // constructor for etrue-var1-var2 splines
  virtual ~splineFDBase();
  //TODO (ETA) - should we make these pure virutal functions? That way each experiment can read in their own spline files in a nice way?
  void SetupSplines();
  void SetupSplines(int BinningOpt);//~~~
  void SetSplineBinning();
  void SetSplineBinning(int BinningOpt);//~~~
  void GetSplineBins(int &nutype, bool &sig, double &enu, double &var1, unsigned int &enu_bin, unsigned int &var1_bin);
  void GetSplineBins(int &nutype, bool &sig, double &enu, double &var1, double &var2, unsigned int &enu_bin, unsigned int &bin_var1, unsigned int &bin_var2);//~~~
  //ETA it's nice to be able to know how many spline bins there are in samplePDFFDBase so adding these.
  //This way we can check that spline and sample binning is the same.
  int getNSplineBinsEnu() {return enu_spline->GetNbins();}
  int getNSplineBinsVar1() {return var1_spline->GetNbins();}
  int getNSplineBinsVar2() {return var2_spline->GetNbins();}

  //DB Function which interregates MaCh3Mode_to_SplineMode to determine which MaCh3Modes have splines and which modes piggy-back of the splines of other modes
  //TODO (ETA) - reimplement this in a generic way for all experiments.
  //This is to become pure virtual and defined in experiment specific implementation of splines
  void FindUniqueModes();
  //DB Function which removes any modes which piggy-back off other modes out of the vector given by covarianceXsec::GetSplineModeVecFromDetID() so they can be added in the XML but won't throw exception when trying to load up splines
  std::vector< std::vector<int> > StripDuplicatedModes(std::vector< std::vector<int> > InputVector);

  template <class T>
    double FastSplineEval(T* spline, const int SplineNumber);
  template <class T>
    bool SetSplineInfoArray(T* spline, int isyst);
  void FindSplineSegment();
  void SetSplineInfoArrays();
  void SetupSplineInfoArray(covarianceXsec * xsec);

  FastSplineInfo * SplineInfoArray;
  covarianceXsec * covxsec;
  int nSplineParams;
  std::vector<int> splineParsIndex;
  std::vector<std::string> splineParsNames;
  // ETA - new function to fill eventSplines 2D vec with splines for each event
  // pass it an array of modes, var1 bins and etrue bins from a samplePDF
  std::vector< std::vector<int> > getEventSplines(int &event, int mode, unsigned int &enu_bin, unsigned int &var1_bin);
  std::vector< std::vector<int> > getEventSplines(int &event, int mode, unsigned int &enu_bin, unsigned int &var1_bin, unsigned int &var2_bin);
  void calcWeights();

  const double* retPointer(int parambin,int modebin,int etruebin,int var1bin) {return &dev_1D_w[parambin][modebin][etruebin][var1bin];}
  const double* retPointer(int parambin,int modebin,int etruebin,int var1bin,int var2bin) {return &dev_2D_w[parambin][modebin][etruebin][var1bin][var2bin];}

  //TODO (ETA) - have a nice print function giving the modes to spline mode mapping etc.

 protected:
  int nutype; // 2 = numu/signue | -2 = numub | 1 = nue | -1 = nueb
  TFile *splinefile;
  
  // spline contents
  TAxis *enu_spline;
  TAxis *var1_spline;
  TAxis *var2_spline;
  
  int number_parms;
  int BinningOpt;
  int SampleDetID;

  //DB Variables related to determined which modes have splines and which piggy-back of other modes
  int nUniqueModes;
  std::vector<std::string> UniqueModeFarSplineNames;
  std::vector<int> DuplicatedFDModes;
  std::vector<int> MaCh3Mode_SplineMode_Map;

#if USE_SPLINE_FD == USE_TSpline3_FD

  std::vector< std::vector< std::vector< std::vector< TSpline3* > > > > dev_1D_vec;
  std::vector< std::vector< std::vector< std::vector< std::vector< TSpline3* > > > > > dev_2D_vec;

#elif USE_SPLINE_FD == USE_TSpline3_red_FD

  std::vector< std::vector< std::vector< std::vector< TSpline3_red* > > > > dev_1D_vec;
  std::vector< std::vector< std::vector< std::vector< std::vector< TSpline3_red* > > > > > dev_2D_vec;

#endif

  //Store weights for each eval spline
  std::vector< std::vector< std::vector< std::vector< double > > > > dev_1D_w;
  std::vector< std::vector< std::vector< std::vector< std::vector< double > > > > > dev_2D_w;

  //This basically just keeps a collection of one spline parameter
  //together with the name of the spline and is used to get the
  //splines out of the spline root file
  //There are two implementations depending on whether you're using
  //reduced TSpline3s or not.
#if USE_SPLINE_FD == USE_TSpline3_FD
  //First a struct to hold Enu-Var1 splines
  struct syst{
	std::string name;
	std::vector< std::vector< std::vector< TSpline3* > > > * spline;
  public:
    syst(std::string namein,    std::vector< std::vector< std::vector< TSpline3* > > >* splinein){
      name=namein;
      spline=splinein;
    };
  };

  //First a struct to hold Enu-Var1-Var-2 splines
  struct syst2D{
    std::string name;
    std::vector< std::vector< std::vector< std::vector< TSpline3* > > > >* spline;
	public:
    syst2D(std::string namein,    std::vector<std::vector< std::vector< std::vector< TSpline3* > > > >* splinein){
      name=namein;
      spline=splinein;
    };
  };

#elif USE_SPLINE_FD == USE_TSpline3_red_FD
  //First a struct to hold Enu-Var1 splines
  struct syst{
	std::string name;
	std::vector< std::vector< std::vector< TSpline3_red* > > > * spline;
  public:
    syst(std::string namein,    std::vector< std::vector< std::vector< TSpline3_red* > > >* splinein){
      name=namein;
      spline=splinein;
    };
  };

  //First a struct to hold Enu-Var1-Var-2 splines
  struct syst2D{
    std::string name;
    std::vector< std::vector< std::vector< std::vector< TSpline3_red* > > > >* spline;
	public:
    syst2D(std::string namein,    std::vector<std::vector< std::vector< std::vector< TSpline3_red* > > > >* splinein){
      name=namein;
      spline=splinein;
    };
  };
#endif
};
#endif
