#pragma once

/// Run low or high memory versions of structs
/// N.B. for 64 bit systems sizeof(float) == sizeof(double) so not a huge effect
/// KS: Need more testing on FD
#ifdef _LOW_MEMORY_STRUCTS_
/// Custom floating point (float or double)
#define _float_ float
/// Custom integer (int or short int)
#define _int_ short int
/// Custom unsigned integer (unsigned short int or unsigned int)
#define _unsigned_int_ unsigned short int
#else
#define _float_ double
#define _int_ int
#define _unsigned_int_ unsigned int
#endif

/// KS: noexcept can help with performance but is terrible for debugging, this is meant to help easy way of of turning it on or off. In near future move this to struct or other central class.
//#define SafeException
#ifndef SafeException
#define _noexcept_ noexcept
#else
#define _noexcept_
#endif

/// KS: Using restrict limits the effects of pointer aliasing, aiding optimizations. While reading I found that there might be some compilers which don't have __restrict__. As always we use _restrict_ to more easily turn off restrict in the code
#ifndef DEBUG
#define _restrict_ __restrict__
#else
#define _restrict_
#endif

/// Number of overflow bins in TH2Poly,
#define _TH2PolyOverflowBins_ 9

/// Include some healthy defines for constructors
#define _BAD_DOUBLE_ -999.99
#define _BAD_INT_ -999

#define _DEFAULT_RETURN_VAL_ -999999.123456

// C++ includes
#include <sstream>
#include <fstream> 
#include <iostream>
#include <vector>
#include <iomanip>
#include <set>
#include <list>
#include <unordered_map>

// ROOT include
#include "TSpline.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TObjString.h"
#include "TH2Poly.h"
#include "TFile.h"

#ifdef MULTITHREAD
#include "omp.h"
#endif

// *******************
/// @brief Template to make vector out of an array of any length
template< typename T, size_t N >
std::vector<T> MakeVector( const T (&data)[N] ) {
// *******************
  return std::vector<T>(data, data+N);
}

// *******************
/// @brief KS: This is mad way of converting string to int. Why? To be able to use string with switch
constexpr unsigned int str2int(const char* str, int h = 0) {
// *******************
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

// *******************
/// @brief ETA - Normalisations for cross-section parameters
/// Carrier for whether you want to apply a systematic to an event or not
struct XsecNorms4 {
// *******************
    /// Name of parameters
    std::string name;
    /// Mode which parameter applies to
    std::vector<int> modes;
    /// Horn currents which parameter applies to
    std::vector<int> horncurrents;
    /// PDG which parameter applies to
    std::vector<int> pdgs;
    /// Preosc PDG which parameter applies to
    std::vector<int> preoscpdgs;
    /// Targets which parameter applies to
    std::vector<int> targets;
    /// Does this parameter have kinematic bounds
    bool hasKinBounds;
    /// Generic vector contain enum relating to a kinematic variable
    /// and lower and upper bounds. This can then be passed to IsEventSelected
    std::vector< std::vector<double> > Selection;

    /// Generic vector containing the string of kinematic type
    /// This then needs to be converted to a kinematic type enum
    /// within a samplePDF daughter class
    /// The bounds for each kinematic variable are given in Selection
    std::vector< std::string > KinematicVarStr;

    /// Parameter number of this normalisation in current systematic model
    int index;
};

/// Make an enum of the spline interpolation type
enum RespFuncType {
  kTSpline3_red,  //!< Uses TSpline3_red for interpolation
  kTF1_red,       //!< Uses TF1_red for interpolation
  kRespFuncTypes  //!< This only enumerates
};

/// Make an enum of the spline interpolation type
enum SplineInterpolation {
  kTSpline3,             //!< Default TSpline3 interpolation
  kLinear,               //!< Linear interpolation between knots
  kMonotonic,            //!< EM: DOES NOT make the entire spline monotonic, only the segments
  kAkima,                //!< EM: Akima spline iis allowed to be discontinuous in 2nd derivative and coefficients in any segment
  kLinearFunc,           //!< Liner interpolation using TF1 not spline
  kSplineInterpolations  //!< This only enumerates
};

// **************************************************
/// @brief Get function for TF1_red
/// @param i Interpolation type
inline std::string GetTF1(const SplineInterpolation i) {
  // **************************************************
  std::string Func = "";
  switch(i) {
    case SplineInterpolation::kLinearFunc:
      Func = "([1]+[0]*x)";
      break;
    default:
      std::cerr << "UNKNOWN SPECIFIED!" << std::endl;
      std::cerr << "You gave  " << static_cast<int>(i) << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
  return Func;
}

// **************************************************
/// @brief Convert a RespFuncType type to a SplineInterpolation
/// @param i Interpolation type
inline RespFuncType SplineInterpolation_ToRespFuncType(const SplineInterpolation i) {
// **************************************************
  RespFuncType Type = kRespFuncTypes;
  switch(i) {
    //  TSpline3 (third order spline in ROOT)
    case SplineInterpolation::kTSpline3:
      Type = RespFuncType::kTSpline3_red;
      break;
    case SplineInterpolation::kLinear:
      Type = RespFuncType::kTSpline3_red;
      break;
    case SplineInterpolation::kMonotonic:
      Type = RespFuncType::kTSpline3_red;
      break;
    //  (Experimental) Akima_Spline (crd order spline which is allowed to be discontinuous in 2nd deriv)
    case SplineInterpolation::kAkima:
      Type = RespFuncType::kTSpline3_red;
      break;
    case SplineInterpolation::kLinearFunc:
      Type = RespFuncType::kTF1_red;
      break;
    default:
      std::cerr << "UNKNOWN SPLINE INTERPOLATION SPECIFIED!" << std::endl;
      std::cerr << "You gave  " << static_cast<int>(i) << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
  return Type;
}

// **************************************************
/// @brief Convert a LLH type to a string
inline std::string SplineInterpolation_ToString(const SplineInterpolation i) {
// **************************************************
  std::string name = "";
  switch(i) {
    //  TSpline3 (third order spline in ROOT)
    case SplineInterpolation::kTSpline3:
      name = "TSpline3";
      break;
    case SplineInterpolation::kLinear:
      name = "Linear";
      break;
    case SplineInterpolation::kMonotonic:
      name = "Monotonic";
      break;
    //  (Experimental) Akima_Spline (crd order spline which is allowed to be discontinuous in 2nd deriv)
    case SplineInterpolation::kAkima:
      name = "Akima";
      break;
    case SplineInterpolation::kLinearFunc:
      name = "LinearFunc";
      break;
    default:
      std::cerr << "UNKNOWN SPLINE INTERPOLATION SPECIFIED!" << std::endl;
      std::cerr << "You gave  " << static_cast<int>(i) << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
  return name;
}

/// Make an enum of systematic type recognised by covariance class
/// @todo KS: Consider using enum class, it is generally recommended as safer. It will require many static_cast
enum SystType {
  kNorm,      //!< For normalisation parameters
  kSpline,    //!< For splined parameters (1D)
  kFunc,      //!< For functional parameters
  kSystTypes  //!< This only enumerates
};

// *******************
/// @brief KS: Struct holding info about Spline Systematics
struct XsecSplines1 {
// *******************
  /// Spline interpolation vector
  SplineInterpolation _SplineInterpolationType;

  /// Modes to which spline applies (valid only for binned splines)
  std::vector<int> _fSplineModes;

  /// EM: Cap spline knot lower value
  double _SplineKnotLowBound;
  /// EM: Cap spline knot higher value
  double _SplineKnotUpBound;
};

// **************************************************
/// @brief Convert a Syst type type to a string
inline std::string SystType_ToString(const SystType i) {
// **************************************************
  std::string name = "";
  switch(i) {
    case SystType::kNorm:
      name = "Norm";
      break;
    case SystType::kSpline:
      name = "Spline";
      break;
    case SystType::kFunc:
      name = "Functional";
      break;
    default:
      std::cerr << "UNKNOWN SYST TYPE SPECIFIED!" << std::endl;
      std::cerr << "You gave  " << static_cast<int>(i) << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
  return name;
}

// ***************************
// A handy namespace for variables extraction
namespace MaCh3Utils {
// ***************************

  // ***************************
  /// @brief Return mass for given PDG
  double GetMassFromPDG(int PDG);
  // ***************************

  extern std::unordered_map<int,int> KnownDetIDsMap;
  extern int nKnownDetIDs;

} // end MaCh3Utils namespace

// *****************
/// Enum to track the target material
enum TargetMat {
// *****************
  kTarget_H  = 1,   //!< Hydrogen
  kTarget_C  = 12,  //!< Carbon 12
  kTarget_N  = 14,
  kTarget_O  = 16,  //!< Oxygen 16
  kTarget_Al = 27,
  kTarget_Ar = 40,
  kTarget_Ti = 48,
  kTarget_Fe = 56,
  kTarget_Pb = 207
};

// *****************
/// @brief Converted the Target Mat to a string
inline std::string TargetMat_ToString(const TargetMat i) {
// *****************
  std::string name;

  switch(i) {
    case kTarget_H:
      name = "Hydrogen";
      break;
    case kTarget_C:
      name = "Carbon";
      break;
    case kTarget_N:
      name = "Nitrogen";
      break;
    case kTarget_O:
      name = "Oxygen";
      break;
    case kTarget_Al:
      name = "Aluminium";
      break;
    case kTarget_Ti:
      name = "Titanium";
      break;
    case kTarget_Fe:
      name = "Iron";
      break;
    case kTarget_Pb:
      name = "Lead";
      break;
    default:
      name = "TargetMat_Undefined";
      break;
  }

  return name;
}

// *****************
/// Enum to track the incoming neutrino species
enum NuPDG {
// *****************
  kNue = 12,
  kNumu = 14,
  kNutau = 16,
  kNueBar = -12,
  kNumuBar = -14,
  kNutauBar = -16
};

// *****************
/// Enum to track neutrino species for Prob3
enum ProbNu {
// *****************
  kProbNue = 1,
  kProbNumu = 2,
  kProbNutau = 3,
  kProbNueBar = -1,
  kProbNumuBar = -2,
  kProbNutauBar = -3
};

/// @brief ETA - Probs3++ doesn't use the PDG codes for the neutrino type so add in a small converter
inline int PDGToProbs(NuPDG pdg){

  int ReturnProbNu = -999;

  switch (pdg){
    case kNue:
	  ReturnProbNu = kProbNue;
	  break;
    case kNumu:
	  ReturnProbNu = kProbNumu;
	  break;
    case kNutau:
	  ReturnProbNu = kProbNutau;
	  break;
    case kNueBar:
	  ReturnProbNu = kProbNueBar;
	  break;
	case kNumuBar:
	  ReturnProbNu = kProbNumuBar;
	  break;
	case kNutauBar:
	  ReturnProbNu = kProbNutauBar;
	  break;
	default:
	  std::cout << "Unrecognised pdg for the neutrino so can't map this to an int for Prob3++" << std::endl;
	  break;
  }

  return ReturnProbNu;
}

inline int ProbsToPDG(ProbNu NuType){

  int ReturnNuPDG = -999;

  switch (NuType){
	case kProbNue:
	  ReturnNuPDG = static_cast<int>(kNue);
	  break;
	case kProbNumu:
	  ReturnNuPDG = static_cast<int>(kNumu);
	  break;
	case kProbNutau:
	  ReturnNuPDG = static_cast<int>(kNutau);
	  break;
	case kProbNueBar:
	  ReturnNuPDG = static_cast<int>(kNueBar);
	  break;
	case kProbNumuBar:
	  ReturnNuPDG = static_cast<int>(kNumuBar);
	  break;
	case kProbNutauBar:
	  ReturnNuPDG = static_cast<int>(kNutauBar);
	  break;
	default:
	  std::cout << "Unrecognised NuType for the neutrino so can't map this to a PDG code" << std::endl;
	  break;
  }

  return ReturnNuPDG;
}

/// Make an enum of the test statistic that we're using
enum TestStatistic {
  kPoisson,                //!< Standard Poisson likelihood
  kBarlowBeeston,          //!< Barlow-Beeston following Conway \cite Conway:2011in
  kIceCube,                //!< Based on \cite Arguelles:2019izp
  kPearson,                //!< Standard Pearson likelihood
  kDembinskiAbdelmottele,  //!< Based on \cite Dembinski:2022ios
  kNTestStatistics         //!< Number of test statistics
};

// **************************************************
/// @brief Convert a LLH type to a string
inline std::string TestStatistic_ToString(TestStatistic i) {
// **************************************************
  std::string name = "";

  switch(i) {
    case kPoisson:
    name = "Poisson";
    break;
    case kBarlowBeeston:
    name = "BarlowBeeston";
    break;
    case kIceCube:
    name = "IceCube";
    break;
    case kPearson:
    name = "Pearson";
    break;
    case kDembinskiAbdelmottele:
    name = "DembinskiAbdelmottele";
    break;
    default:
      std::cerr << "UNKNOWN LIKELHOOD SPECIFIED!" << std::endl;
      std::cerr << "You gave test-statistic " << i << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
  return name;
}

/// @brief WP: Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly* poly);

/// @brief WP: Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly* poly);

/// @brief WP: Poly Projectors
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins, bool computeErrors = false);
/// @brief WP: Poly Projectors
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins, bool computeErrors = false);

/// @brief KS: Convert TH2D to TH2Poly
TH2D* ConvertTH2PolyToTH2D(TH2Poly *poly, TH2D *TH2Dhist);
/// @brief KS: Convert TH2Poly to TH2D
TH2Poly* ConvertTH2DtoTH2Poly(TH2D *TH2Dhist);

/// @brief WP: Helper to Normalise histograms
TH2Poly* NormalisePoly(TH2Poly* Histogram);

/// @brief WP: Helper to scale th2poly analogous to th2d scale with option "width"
TH2Poly* PolyScaleWidth(TH2Poly *Histogram, double scale);
/// @brief WP: Helper to calc integral of th2poly analogous to th2d integra; with option "width"
double PolyIntegralWidth(TH2Poly *Histogram);

/// @brief KS: ROOT changes something with binning when moving from ROOT 5 to ROOT 6. If you open ROOT5 produced file with ROOT6 you will be missing 9 last bins
/// @brief However if you use ROOT6 and have ROOT6 file exactly the same code will work. Something have changed with how TH2Poly bins are stored in TFile
/// @param file ROOT file that we will make version checks
void CheckTH2PolyFileVersion(TFile *file);

/// @brief KS: Remove fitted TF1 from hist to make comparison easier
void RemoveFitter(TH1D* hist, const std::string& name);

/// @brief Helper to check if files exist or not
inline std::string file_exists(std::string filename) {
  std::ifstream infile(filename.c_str());
  if (!infile.good()) {
    std::cerr << "*** ERROR ***" << std::endl;
    std::cerr << "File " << filename << " does not exist" << std::endl;
    std::cerr << "Please try again" << std::endl;
    std::cerr << "*************" << std::endl;
    throw;
  }

  return filename;
}
/// @brief DB Get the Cernekov momentum threshold in MeV
double returnCherenkovThresholdMomentum(int PDG);

double CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2 = 0.0);
double CalculateEnu(double PLep, double cosTheta, double EB, bool neutrino);


enum CUDAProb_nu {
  e_e = 0,
  e_m = 1,
  e_t = 2,
  m_e = 3,
  m_m = 4,
  m_t = 5,
  t_e = 6,
  t_m = 7,
  t_t = 8
};
 

// ************************************************
/// @brief Get CUDAProb3 flavour from intital and final states
inline CUDAProb_nu GetCUDAProbFlavour(int nu_i, int nu_f) {
//*************************************************  
    
  switch (abs(nu_i)) {
  case 1:
    switch (abs(nu_f)) {
    case 1:
      return CUDAProb_nu::e_e;
      break;
    case 2:
      return CUDAProb_nu::e_m;
      break;
    case 3:
      return CUDAProb_nu::e_t;
      break;
	default:
	  std::cout << "Unknow flavour " << nu_f << std::endl;
	  throw;
    } 
  case 2:
    switch (abs(nu_f)) {
    case 1:
      return CUDAProb_nu::m_e;
      break;
    case 2:
      return CUDAProb_nu::m_m;
      break;
    case 3:
      return CUDAProb_nu::m_t;
      break;
	default:
	  std::cout << "Unknow flavour " << nu_f << std::endl;
	  throw;
    } 
  case 3:
    switch (abs(nu_f)) {
    case 1:
      return CUDAProb_nu::t_e;
      break;
    case 2:
      return CUDAProb_nu::t_m;
      break;
    case 3:
      return CUDAProb_nu::t_t;
      break;
	default:
	  std::cout << "Unknow flavour " << nu_f << std::endl;
	  throw;
    }
  default:
	std::cout << "Unknow flavour " << nu_i << std::endl;
	throw;
  }

}
