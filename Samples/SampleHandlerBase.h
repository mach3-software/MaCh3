#pragma once

//C++ includes
#include <assert.h>

//MaCh3 includes
#include "Samples/SampleStructs.h"
#include "Samples/HistogramUtils.h"
#include "Manager/Manager.h"
#include "Manager/MaCh3Modes.h"

_MaCh3_Safe_Include_Start_ //{
//ROOT includes
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
_MaCh3_Safe_Include_End_ //}

/// @brief Class responsible for handling implementation of samples used in analysis, reweighting and returning LLH
class SampleHandlerBase
{
 public:
   /// @brief The main constructor
   SampleHandlerBase();
  /// @brief destructor
  virtual ~SampleHandlerBase();

  /// @defgroup SampleHandlerSetters Sample Handler Setters
  /// Group of functions to set various parameters, names, and values.

  /// @defgroup SampleHandlerGetters Sample Handler Getters
  /// Group of functions to get various parameters, names, and values.

  /// @ingroup SampleHandlerGetters
  virtual inline M3::int_t GetNsamples(){ return nSamples; };
  /// @ingroup SampleHandlerGetters
  virtual inline std::string GetTitle()const {return "SampleHandler";};
  /// @ingroup SampleHandlerGetters
  virtual std::string GetSampleName(int Sample) const = 0;
  /// @ingroup SampleHandlerGetters
  virtual inline double GetSampleLikelihood(const int isample){(void) isample; return GetLikelihood();};
  /// @brief Allow to clean not used memory before fit starts
  virtual void CleanMemoryBeforeFit() = 0;
  /// @brief Store additional info in a chan
  virtual void SaveAdditionalInfo(TDirectory* Dir) {(void) Dir;};
  /// @brief Return pointer to MaCh3 modes
  /// @ingroup SampleHandlerGetters
  MaCh3Modes* GetMaCh3Modes() const { return Modes.get(); }
      
  virtual void Reweight()=0;
  /// @ingroup SampleHandlerGetters
  virtual double GetLikelihood() = 0;

  /// @ingroup SampleHandlerGetters
  unsigned int GetNEvents(){return nEvents;}
  /// @ingroup SampleHandlerGetters
  virtual int GetNMCSamples() { return nSamples; }
  /// @ingroup SampleHandlerGetters
  virtual int GetNOscChannels(){ return 1; }

  // WARNING KS: Needed for sigma var
  virtual void SetupBinning(const M3::int_t Selection, std::vector<double> &BinningX, std::vector<double> &BinningY){
    (void) Selection; (void) BinningX; (void) BinningY; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented");}
  virtual TH1* GetData(const int Selection) { (void) Selection; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }
  virtual TH2Poly* GetW2(const int Selection){ (void) Selection; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented");}
  virtual TH1* GetPDF(const int Selection){ (void) Selection; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented");}
  virtual inline TH1* GetPDFMode(const int Selection, const int Mode) {
    (void) Selection; (void) Mode; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }
  virtual inline std::string GetKinVarLabel(const int sample, const int Dimension) {
    (void) sample; (void) Dimension; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented");  };

  /// @brief Calculate test statistic for a single bin using Poisson
  /// @param data is data
  /// @param mc is mc
  /// @ingroup SampleHandlerGetters
  double GetTestStatLLH(double data, double mc) const;
  /// @brief Calculate test statistic for a single bin. Calculation depends on setting of fTestStatistic. Data and mc -> 0 cut-offs are defined in M3::_LOW_MC_BOUND_.
  /// @details Implemented fTestStatistic are kPoisson (with Stirling's approx.), kBarlowBeeston (arXiv:1103.0354), kDembinskiAbdelmotteleb (arXiv:2206.12346), kIceCube (arxiv:1901.04645), and kPearson.
  /// Test statistics require mc > 0, therefore low mc and data values are treated with cut-offs based on M3::_LOW_MC_BOUND_ = .00001 by default.
  /// For kPoisson, kBarlowBeeston, kDembinskiAbdelmotteleb, kPearson:
  /// data > _LOW_MC_BOUND_ & mc <= _LOW_MC_BOUND_: returns GetTestStatLLH(data, _LOW_MC_BOUND_, w2), with Poisson(data,_LOW_MC_BOUND_) limit for mc->0, w2->0.
  /// mc < data <= _LOW_MC_BOUND_: returns 0 (as if any data <= _LOW_MC_BOUND_ were effectively consistent with 0 data count), with a limit of 0 for mc->0.
  /// data = 0: returns mc (or mc/2. for kPearson), with a limit of 0 for mc->0.
  /// For kIceCube:
  /// mc < data returns the lower of IceCube(data,mc,w2) and Poisson(data,mc) penalties, with a Poisson(data,_LOW_MC_BOUND_) limit for mc->0, w2->0.
  /// @param data is data
  /// @param mc is mc
  /// @param w2 is is Sum(w_{i}^2) (sum of weights squared), which is sigma^2_{MC stats}
  /// @ingroup SampleHandlerGetters
  double GetTestStatLLH(const double data, const double mc, const double w2) const;
  /// @brief Set the test statistic to be used when calculating the binned likelihoods
  /// @param testStat The test statistic to use.
  /// @ingroup SampleHandlerGetters
  inline void SetTestStatistic(TestStatistic testStat){ fTestStatistic = testStat; }

protected:
  /// @brief CW: Redirect std::cout to silence some experiment specific libraries
  void QuietPlease();
  /// @brief CW: Redirect std::cout to silence some experiment specific libraries
  void NowTalk();

  /// @brief check if event is affected by following conditions, for example pdg, or modes etc
  template <typename T>
  bool MatchCondition(const std::vector<T>& allowedValues, const T& value) {
    if (allowedValues.empty()) {
      return true;  // Apply to all if no specific values are specified
    }
    return std::find(allowedValues.begin(), allowedValues.end(), value) != allowedValues.end();
  }

  /// Test statistic tells what kind of likelihood sample is using
  TestStatistic fTestStatistic;

  /// Keep the cout buffer
  std::streambuf *buf;
  /// Keep the cerr buffer
  std::streambuf *errbuf;

  /// Contains how many samples we've got
  M3::int_t nSamples;

  /// Number of MC events are there
  unsigned int nEvents;

  /// Holds information about used Generator and MaCh3 modes
  std::unique_ptr<MaCh3Modes> Modes;
};
