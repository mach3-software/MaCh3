#pragma once

//C++ includes
#include <assert.h>

//MaCh3 includes
#include "Samples/Structs.h"
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
   SampleHandlerBase();
  /// @brief destructor
  virtual ~SampleHandlerBase();

  virtual inline M3::int_t GetNsamples(){ return nSamples; };
  virtual inline std::string GetTitle()const {return "SampleHandler";};
  virtual std::string GetSampleName(int Sample) const = 0;
  virtual inline double GetSampleLikelihood(const int isample){(void) isample; return GetLikelihood();};

  /// @brief Return pointer to MaCh3 modes
  MaCh3Modes* GetMaCh3Modes() const { return Modes; }

  TH1D* Get1DHist();                                               
  TH2D* Get2DHist();
  TH1D* Get1DDataHist(){return dathist;}
  TH2D* Get2DDataHist(){return dathist2d;}
      
  virtual void Reweight()=0;
  virtual double GetLikelihood() = 0;

  virtual int GetNEventsInSample(int sample){ (void) sample; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }
  virtual int GetNMCSamples(){ throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }

  virtual void AddData(std::vector<double> &dat);
  virtual void AddData(std::vector< std::vector <double> > &dat);
  virtual void AddData(TH1D* binneddata);
  virtual void AddData(TH2D* binneddata);

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

  double GetTestStatLLH(double data, double mc) const;
  /// @brief Calculate test statistic for a single bin. Calculation depends on setting of fTestStatistic
  /// @param data is data
  /// @param mc is mc
  /// @param w2 is is Sum(w_{i}^2) (sum of weights squared), which is sigma^2_{MC stats}
  double GetTestStatLLH(const double data, const double mc, const double w2) const;
  /// @brief Set the test statistic to be used when calculating the binned likelihoods
  /// @param testStat The test statistic to use.
  inline void SetTestStatistic(TestStatistic testStat){ fTestStatistic = testStat; }

  virtual void Fill1DHist()=0;
  virtual void Fill2DHist()=0;

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
  /// KS: number of dimension for this sample
  int nDims;

  /// Holds information about used Generator and MaCh3 modes
  MaCh3Modes* Modes;

  TH1D *dathist; // tempstore for likelihood calc
  TH2D *dathist2d;

  // binned PDFs
  TH1D*_hPDF1D;
  TH2D*_hPDF2D;

  TRandom3* rnd;
};
