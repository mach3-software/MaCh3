#pragma once

//C++ includes
#include <assert.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
//ROOT includes
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#pragma GCC diagnostic pop

//MaCh3 includes
#include "samplePDF/Structs.h"
#include "samplePDF/HistogramUtils.h"
#include "manager/manager.h"

/// @brief Class responsible for handling implementation of samples used in analysis, reweighting and returning LLH
class samplePDFBase
{
 public:
   samplePDFBase();
  /// @brief destructor
  virtual ~samplePDFBase();

  virtual inline M3::int_t GetNsamples(){ return nSamples; };
  virtual inline std::string GetName()const {return "samplePDF";};
  virtual std::string GetSampleName(int Sample);
  virtual inline double getSampleLikelihood(const int isample){(void) isample; return GetLikelihood();};

  /// @brief Return pointer to MaCh3 modes
  MaCh3Modes* GetMaCh3Modes() const { return Modes; }

  TH1D* get1DHist();                                               
  TH2D* get2DHist();
  TH1D* get1DDataHist(){return dathist;}
  TH2D* get2DDataHist(){return dathist2d;}
      
  virtual void reweight()=0;
  virtual double GetLikelihood() = 0;

  virtual int getNEventsInSample(int sample){ (void) sample; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }
  virtual int getNMCSamples(){ throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }

  virtual void addData(std::vector<double> &dat);
  virtual void addData(std::vector< std::vector <double> > &dat);
  virtual void addData(TH1D* binneddata);
  virtual void addData(TH2D* binneddata);

  // WARNING KS: Needed for sigma var
  virtual void SetupBinning(const M3::int_t Selection, std::vector<double> &BinningX, std::vector<double> &BinningY){
    (void) Selection; (void) BinningX; (void) BinningY; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented");}
  virtual TH1* getData(const int Selection) { (void) Selection; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }
  virtual TH2Poly* getW2(const int Selection){ (void) Selection; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented");}
  virtual TH1* getPDF(const int Selection){ (void) Selection; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented");}
  virtual inline TH1* getPDFMode(const int Selection, const int Mode) {
    (void) Selection; (void) Mode; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }
  virtual inline std::string GetKinVarLabel(const int sample, const int Dimension) {
    (void) sample; (void) Dimension; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented");  };

  double getTestStatLLH(double data, double mc) const;
  /// @brief Calculate test statistic for a single bin. Calculation depends on setting of fTestStatistic
  /// @param data is data
  /// @param mc is mc
  /// @param w2 is is Sum(w_{i}^2) (sum of weights squared), which is sigma^2_{MC stats}
  double getTestStatLLH(const double data, const double mc, const double w2) const;
  /// @brief Set the test statistic to be used when calculating the binned likelihoods
  /// @param testStat The test statistic to use.
  inline void SetTestStatistic(TestStatistic testStat){ fTestStatistic = testStat; }

  virtual void fill1DHist()=0;
  virtual void fill2DHist()=0;

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
  /// Name of Sample
  std::vector<std::string> SampleName;

  //GetterForModes
  MaCh3Modes* Modes;

  TH1D *dathist; // tempstore for likelihood calc
  TH2D *dathist2d;

  // binned PDFs
  TH1D*_hPDF1D;
  TH2D*_hPDF2D;

  TRandom3* rnd;

  double pot;
};
