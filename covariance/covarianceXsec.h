#ifndef _covarianceXsec_h_

#define _covarianceXsec_h_

// C++ includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <map>

// ROOT includes
#include <TDecompChol.h>
#include <TList.h>
#include <TStopwatch.h>
#include <TMatrixT.h>
#include <TVectorT.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TFile.h>
#include <TAxis.h>
#include <TSpline.h>

// MaCh3 includes
#include "covarianceBase.h"
#include "samplePDF/Structs.h"
#include "throwParms/ThrowParms.h"

class covarianceXsec : public covarianceBase {

  public:
  covarianceXsec(const char *name, const char *file, double threshold=-1,int firstpcapar=-999,int lastpcapar=-999);
    ~covarianceXsec();

    // Print information about the whole object once it is set
    void Print();

    void throwNominal(bool nomValues=true, int seed = 0);

    // General Getter functions not split by detector
    double GetLikelihood();
    const double GetParamUpperBound(int i) {return xsec_param_ub_a[i];}
    const double GetParamLowerBound(int i) {return xsec_param_lb_a[i];}
    const double GetParamPrior(int i)      {return xsec_param_prior_a[i];}
    const int  GetXSecParamID(int i, int j) const {return xsec_param_id_a[i][j];}
    const std::string & GetParameterName(int i) const {return xsec_param_names[i];}
    const int    GetNumParams()               {return nPars;}
    const char* GetParName(int i) const {return xsec_param_names[i].c_str();}

    const bool isParFlux(int i){
      return isFlux[i];
    }

    // Get functions for Near normalisation parameters
    const std::vector<XsecNorms4> GetNearNormPars() const{return NearNormParams;}
    const int                       GetNumNearNormParams() const  {return nNearNormParams;}

    // Get functions for Far normalisation parameters
    const std::vector<XsecNorms4> GetFarNormPars() const{return FarNormParams;}

    // Get functions for Near spline parameters
    const int                       GetNumNearSplineParams() const  {return nNearSplineParams;}
    const std::vector<std::string>& GetNearSplineParsNames() const  {return NearsplineParsNames;}
    const std::vector<std::string>& GetNearSplineFileParsNames() const  {return NearSplineFileParsNames;}
    const std::vector<int>&         GetNearSplineParsIndex() const  {return NearsplineParsIndex;}

    //DB Get spline parameters depending on given DetID
    const std::vector<std::string> GetSplineParsNamesFromDetID(int DetID);
    const std::vector<std::string> GetSplineFileParsNamesFromDetID(int DetID);
    const std::vector< std::vector<int> > GetSplineModeVecFromDetID(int DetID);
    const std::vector<int> GetSplineParsIndexFromDetID(int DetID);
    const int GetNumSplineParamsFromDetID(int DetID);

    //DB Get norm/func parameters depending on given DetID
    const std::vector<XsecNorms4> GetNormParsFromDetID(int DetID);
    const int GetNumFuncParamsFromDetID(int DetID);
    const std::vector<std::string> GetFuncParsNamesFromDetID(int DetID);
    const std::vector<int> GetFuncParsIndexFromDetID(int DetID);

    // Get functions for Far spline parameters
    const int                       GetNumFarSplineParams() const  {return nFarSplineParams;}
    const std::vector<std::string>& GetFarSplineParsNames() const  {return FarSplineParsNames;}
    const std::vector<std::string>& GetFarSplineFileParsNames() const  {return FarSplineFileParsNames;}
    //TVectorT<double> *GetFarSplineParsModes_TVec() {return FarSplineParsModes;}
    //const double GetFarSplineParsModes(int i) {return xsec_fd_spline_mode_a[i];}
    const std::vector<int>&         GetFarSplineParsIndex() const  {return FarSplineParsIndex;}
    const std::vector<int>&         GetFarSplineModeVec(int i) const {return FarSplineModes[i];}
	
    // Get functions for uniq spline parameters //!!decide what to do about these
    const int                       GetNumSplineParamsUniq() const  {return nSplineParamsUniq;}
    const std::vector<std::string>& GetSplineParsUniqNames() const  {return splineParsUniqNames;}
    const std::vector<int>&         GetSplineParsUniqIndex() const  {return splineParsUniqIndex;}

    // Get functions for shared spline parameters //!!decide what to do about these
    const int                       GetNumSplineParamsShare() const  {return nSplineParamsShare;}
    const std::vector<std::string>& GetSplineParsShareNames() const  {return splineParsShareNames;}
    const std::vector<int>&         GetSplineParsShareIndex() const  {return splineParsShareIndex;}
    const std::vector<int>&         GetSplineParsShareToUniq() const {return splineParsShareToUniq;}
    
    // Get functions for Near functional parameter (e.g. BeRPA)
    const int                       GetNumNearFuncParams() const     {return nNearFuncParams;}
    const std::vector<std::string>& GetNearFuncParsNames() const  {return NearfuncParsNames;}
    const std::vector<int>&         GetNearFuncParsIndex() const  {return NearfuncParsIndex;}

    // Get functions for Far functional parameter (e.g. BeRPA)
    const int                       GetNumFarFuncParams() const     {return nFarFuncParams;}
    const std::vector<std::string>& GetFarFuncParsNames() const  {return FarFuncParsNames;}
    const std::vector<int>&         GetFarFuncParsIndex() const  {return FarFuncParsIndex;}


    // Get nominal and prior for saving to the output file
    TVectorT<double> *GetNominal_TVec() {return xsec_param_nom;}
    TVectorT<double> *GetPrior_TVec()   {return xsec_param_prior;}


    // If we want to over-ride the default of running with a Gaussian prior on parameter i
    void setEvalLikelihood(int i, bool eL);
    void toggleFixParameter(int i);

    // What parameter Gets reweighted by what amount according to MCMC
    inline double calcReWeight(int bin){
	  if (bin >= 0 && bin < nPars) {
		return fParProp[bin];
	  } else {
		std::cerr << "Specified bin is <= 0 OR bin > npar!" << std::endl;
		std::cerr << "bin = " << bin << ", npar = " << nPars << std::endl;
		std::cerr << "This won't ruin much that this step in the MCMC, but does indicate something wrong in memory!" << std::endl;
		return 1.0;
	  }

	  return 1.0;  
	}; 

    //int *BoundaryHits;


  protected:
    // Helper functions to decide on what setup we're running
    void scanParameters();

    void initParams(double fScale);
    void setXsecParNames();

    // Vectors of the input root file
    TVectorT<double> *xsec_param_prior;
    TVectorT<double> *xsec_param_nom;
    TVectorT<double> *xsec_param_lb;
    TVectorT<double> *xsec_param_ub;
    TMatrixT<double> *xsec_param_id;

    //DB StepScaleReading
    TVectorT<double> *xsec_stepscale;
    std::vector<double> xsec_stepscale_vec;

    // TObjArrays from the input root file
    TObjArray* xsec_param_norm_modes;
    TObjArray* xsec_param_norm_horncurrents;
    TObjArray* xsec_param_norm_elem;
    TObjArray* xsec_param_norm_nupdg;
    TObjArray* xsec_param_norm_preoscnupdg;
    //TObjArray* xsec_param_norm_etru_bnd_low;
    //TObjArray* xsec_param_norm_etru_bnd_high;
    //TObjArray* xsec_param_norm_q2_true_bnd_low;
    //TObjArray* xsec_param_norm_q2_true_bnd_high;
	TObjArray* xsec_kinematic_type;
	TVectorT<double> *xsec_kinematic_ub;
	TVectorT<double> *xsec_kinematic_lb;

    // Here are some array equivalents (actually used in MCMC)
    int **xsec_param_id_a;
    // nominal values in MC
    double *xsec_param_nom_a;
    // lower bound
    double *xsec_param_lb_a;
    // upper bound
    double *xsec_param_ub_a;
    // priors from external data fit
    double *xsec_param_prior_a;

    // Contains the parameter names
    std::vector<std::string> xsec_param_names;
    
    //Contains the names of the Far spline objects in the spline files
    TObjArray* xsec_param_fd_spline_names;
    TObjArray* xsec_param_fd_spline_modes;

    std::vector<bool> isFlux;

    // Number of total parameters, just scanned from input root file
    int nPars;


    int nTotalNormParams;

    int nFaronlyNormParams;//Needed for consistency check
    int nLowEnergyAtmOnlyNormParams;//DB Needed for consistency check
    int nHighEnergyAtmOnlyNormParams;//DB Needed for consistency check

    std::vector<XsecNorms4> FarNormParams;
    int nFarNormParams;

    std::vector<XsecNorms4> NearNormParams;
    int nNearNormParams;

    // Number of Near spline parameters
    int nNearSplineParams;
    std::vector<std::string> NearsplineParsNames;
    std::vector<std::string> NearSplineFileParsNames;
    std::vector<int> NearsplineParsIndex;

    TObjArray* xsec_param_nd_spline_names;

    int nFarSplineParams;
    std::vector<std::string> FarSplineParsNames;
    std::vector<int> FarSplineParsIndex;
    //ETA - new object for storing spline name in File
    std::vector<std::string> FarSplineFileParsNames;
    //Mode which spline applies to
    //mode=12 (which is kMaCh3_nModes for NIWG2020 model) means it applies to all modes
    //TVectorT<double> *FarSplineParsModes;
    std::vector<std::vector<int> > FarSplineModes;
    //double *xsec_fd_spline_mode_a;

    // Number of spline parameters that aren't repeated
    int nSplineParamsUniq;
    std::vector<std::string> splineParsUniqNames;
    std::vector<int> splineParsUniqIndex;

    // Number of spline parameters that are shared
    int nSplineParamsShare;
    std::vector<std::string> splineParsShareNames;
    std::vector<int> splineParsShareIndex;
    std::vector<int> splineParsShareToUniq;

    // Number of functional parameter
    int nNearFuncParams;
    std::vector<std::string> NearfuncParsNames;
    std::vector<int> NearfuncParsIndex;

    int nFarFuncParams;
    std::vector<std::string> FarFuncParsNames;
    std::vector<int> FarFuncParsIndex;

};

#endif
