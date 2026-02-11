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
/// @ingroup CoreClasses
class SampleHandlerBase
{
 public:
   /// @brief The main constructor
   SampleHandlerBase();
  /// @brief destructor
  virtual ~SampleHandlerBase();

  /// @ingroup SampleHandlerGetters
  virtual inline M3::int_t GetNsamples(){ return nSamples; };
  /// @ingroup SampleHandlerGetters
  virtual std::string GetSampleTitle(const int Sample) const = 0;
  /// @ingroup SampleHandlerGetters
  virtual std::string GetName() const = 0;
  /// @ingroup SampleHandlerGetters
  virtual double GetSampleLikelihood(const int isample) const = 0;
  /// @brief Allow to clean not used memory before fit starts
  virtual void CleanMemoryBeforeFit() = 0;
  /// @brief Store additional info in a chan
  virtual void SaveAdditionalInfo(TDirectory* Dir) {(void) Dir;};
  /// @brief Return pointer to MaCh3 modes
  /// @ingroup SampleHandlerGetters
  MaCh3Modes* GetMaCh3Modes() const { return Modes.get(); }
      
  virtual void Reweight()=0;
  /// @ingroup SampleHandlerGetters
  virtual double GetLikelihood() const = 0;

  /// @brief Helper function to print rates for the samples with LLH
  /// @param DataOnly whether to print data only rates
  virtual void PrintRates(const bool DataOnly = false) = 0;

  /// @ingroup SampleHandlerGetters
  unsigned int GetNEvents() const {return nEvents;}
  /// @ingroup SampleHandlerGetters
  virtual int GetNOscChannels(const int iSample) const = 0;

  /// @brief Return Kinematic Variable name for specified sample and dimension for example "Reconstructed_Neutrino_Energy"
  /// @param iSample Sample index
  /// @param Dimension Dimension index
  /// @ingroup SampleHandlerGetters
  virtual std::string GetKinVarName(const int iSample, const int Dimension) const = 0;

  /// @brief Get Data histogram
  /// @ingroup SampleHandlerGetters
  virtual TH1* GetDataHist(const int Sample) = 0;
  /// @brief Get MC histogram
  /// @ingroup SampleHandlerGetters
  virtual TH1* GetMCHist(const int Sample) = 0;
  /// @brief Get W2 histogram
  /// @ingroup SampleHandlerGetters
  virtual TH1* GetW2Hist(const int Sample) = 0;

  /// @brief DB Function to differentiate 1D or 2D binning
  virtual int GetNDim(const int Sample) const = 0;
  virtual std::string GetFlavourName(const int iSample, const int iChannel) const = 0;

  /// @brief Return the binning used to draw a kinematic parameter
  virtual std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string &KinematicParameter) const = 0;

  virtual TH1* Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str,
                                            int kModeToFill = -1, int kChannelToFill = -1, int WeightStyle = 0, TAxis* Axis = nullptr) = 0;
  virtual TH1 *Get1DVarHist(const int iSample, const std::string &ProjectionVar,
                           const std::vector<KinematicCut> &EventSelectionVec = {}, int WeightStyle = 0,
                            TAxis *Axis = nullptr, const std::vector<KinematicCut> &SubEventSelectionVec = {}) = 0;
  virtual TH2* Get2DVarHist(const int iSample, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                            const std::vector< KinematicCut >& EventSelectionVec = {},
                            int WeightStyle = 0, TAxis* AxisX = nullptr, TAxis* AxisY = nullptr,
                            const std::vector< KinematicCut >& SubEventSelectionVec = {}) = 0;


  // WARNING KS: Needed for sigma var, but also remnants of T2K-ND280 code will be merged in SampleHandlerFD, stay tuned...
  virtual inline TH1* GetPDFMode(const int Selection, const int Mode) {
    (void) Selection; (void) Mode; throw MaCh3Exception(__FILE__ , __LINE__ , "Not implemented"); }
  ///////////

  /// @brief Calculate test statistic for a single bin using Poisson
  /// @param data is data
  /// @param mc is mc
  /// @ingroup SampleHandlerGetters
  double GetPoissonLLH(const double data, const double mc) const;

  /// @brief Calculate test statistic for a single bin. Calculation depends on setting of @p fTestStatistic. Data and mc -> 0 cut-offs are defined in M3::_LOW_MC_BOUND_.
  ///
  /// @details
  ///
  /// ### Poisson
  /// Standard Poisson log-likelihood (Stirling approximation) @cite BakerCousins1984
  /// \f[
  /// - \log \mathcal{L}_\mathrm{Poisson} = \sum_i N_i^\mathrm{MC} - N_i^\mathrm{data} +
  /// N_i^\mathrm{data} \ln \frac{N_i^\mathrm{data}}{N_i^\mathrm{MC}},
  /// \f]
  ///
  /// ### Pearson
  /// Standard Pearson likelihood @cite Pearson1900 (assumes Gaussian approximation of bin counts):
  /// \f[
  /// - \log \mathcal{L}_\mathrm{Pearson} = \sum_i \frac{(N_i^\mathrm{data} - N_i^\mathrm{MC})^2}{2 \, N_i^\mathrm{MC}}
  /// \f]
  ///
  /// ### Barlow-Beeston
  /// Based on @cite Barlow:1993dm and following Conway approximation (@cite Conway:2011in)
  /// The generation of MC is a stochastic process, so even identical settings can lead to different outputs
  /// (assuming that the seeds of the random number generator are different). This introduces uncertainty in
  /// MC distributions, especially in bins with low statistics.
  ///
  /// \f[
  /// - \log \mathcal{L}_\mathrm{BB} = - \log \mathcal{L}_\mathrm{Poisson} - \log \mathcal{L}_\mathrm{MC_{stat}}
  /// =  \sum_i \Biggl[ N_i^\mathrm{MC}(\vec{\theta}) - N_i^\mathrm{data} +
  /// N_i^\mathrm{data} \ln \frac{N_i^\mathrm{data}}{N_i^\mathrm{MC}(\vec{\theta})} +
  /// \frac{(\beta_i - 1)^2}{2 \sigma_{\beta_i}^2} \Biggr],
  /// \f]
  ///
  /// where \f$\beta_i\f$ is a scaling parameter between ideal ("true") and generated MC in a bin
  /// (\f$N^\mathrm{true}_{\mathrm{MC},i} = \beta_i N_i^\mathrm{MC}\f$), and
  /// \f$\sigma^2_{\beta_i} = \frac{\sum_i w_i^2}{N_i^\mathrm{MC}}\f$, with \f$\sum_i w_i^2\f$ being the sum of the
  /// squares of weights in bin \f$i\f$. Assuming \f$\beta_i\f$ follows a Gaussian, its mean can be found by solving
  /// the quadratic equation derived by Conway:
  ///
  /// \f[
  /// \beta_i^2 + (N_i^\mathrm{MC} \sigma_{\beta_i}^2 - 1)\beta_i - N_i^\mathrm{data} \sigma_{\beta_i}^2 = 0
  /// \f]
  ///
  /// ### Dembinski-Abdelmotteleb
  /// Alternative treatment of MC statistical uncertainty following Hans Dembinski and Ahmed Abdelmotteleb @cite Dembinski:2022ios
  ///
  /// This approach extends the Barlow-Beeston method. For each bin:
  /// \f[
  /// - \log \mathcal{L}_\mathrm{DA} = (N_i^{\mathrm{MC},\prime} - N_i^\mathrm{data} +
  /// N_i^\mathrm{data} \ln \frac{N_i^\mathrm{data}}{N_i^{\mathrm{MC},\prime}}) + k \beta - k + k \ln \frac{k}{k \beta}
  /// \f]
  /// where
  /// \f[
  /// k = \frac{(N_i^\mathrm{MC})^2}{\sum_i w_i^2}
  /// \f]
  /// and
  /// \f[
  /// \beta = \frac{N_i^\mathrm{data} + k}{N_i^\mathrm{MC} + k}, \quad
  /// N_i^{\mathrm{MC},\prime} = N_i^\mathrm{MC} \cdot \beta
  /// \f]
  ///
  /// ### IceCube
  /// Alternative likelihood definition described by the IceCube collaboration @cite Arguelles:2019izp
  /// \f[
  /// - \log \mathcal{L} = -  \sum_i \Biggl(
  ///     a_i \log(b_i) + \log[\Gamma(N_i^{\mathrm{data}}+a_i)]
  ///     - (N_i^{\mathrm{data}}+a_i)\log(b_i+1) - \log[\Gamma(a_i)]
  /// \Biggr),
  /// \f]
  /// where the auxiliary variables are
  /// \f[
  /// a_i = N^{\mathrm{gen}}_{\mathrm{MC},i} \, b_i + 1, \quad
  /// b_i = \frac{N^{\mathrm{gen}}_{\mathrm{MC},i}}{\sum_i w_i^2}.
  /// \f]
  ///
  /// ### Treatment of low data/mc
  /// Implemented fTestStatistic are @p kPoisson (with Stirling's approx.), @p kBarlowBeeston (arXiv:1103.0354), @p kDembinskiAbdelmotteleb (arXiv:2206.12346), @p kIceCube (arxiv:1901.04645), and @p kPearson.
  /// Test statistics require mc > 0, therefore low mc and data values are treated with cut-offs based on M3::_LOW_MC_BOUND_ = .00001 by default.
  /// For @p kPoisson, @p kBarlowBeeston, @p kDembinskiAbdelmotteleb, @p kPearson:
  /// data > _LOW_MC_BOUND_ & mc <= _LOW_MC_BOUND_: returns GetTestStatLLH(data, _LOW_MC_BOUND_, w2), with Poisson(data,_LOW_MC_BOUND_) limit for mc->0, w2->0.
  /// mc < data <= _LOW_MC_BOUND_: returns 0 (as if any data <= _LOW_MC_BOUND_ were effectively consistent with 0 data count), with a limit of 0 for mc->0.
  /// data = 0: returns mc (or mc/2. for @p kPearson), with a limit of 0 for mc->0.
  /// For @p kIceCube:
  /// mc < data returns the lower of IceCube(data,mc,w2) and Poisson(data,mc) penalties, with a Poisson(data,_LOW_MC_BOUND_) limit for mc->0, w2->0.
  /// @param data is data
  /// @param mc is mc
  /// @param w2 is \f$\sum_{i} w_{i}^2\f$ (sum of weights squared), which is \f$\sigma^2_{\text{MC stats}}\f$
  /// @ingroup SampleHandlerGetters
  double GetTestStatLLH(const double data, const double mc, const double w2) const;

  /// @brief Set the test statistic to be used when calculating the binned likelihoods
  /// @param testStat The test statistic to use.
  /// @ingroup SampleHandlerGetters
  void SetTestStatistic(TestStatistic testStat){ fTestStatistic = testStat; }
  /// @brief Get the test statistic used when calculating the binned likelihoods
  TestStatistic GetTestStatistic() const { return fTestStatistic; }

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
