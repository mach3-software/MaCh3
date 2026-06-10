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
///
/// @details It serves as basic interface for fit running, as well as other fucnitonalsities liek llh scan, sigma var or event posterior predictive disturibon.
/// Concrete implementations of this interface are responsible for defining
/// the specific structure of samples, event selections, and histogram filling.
///
/// @ingroup CoreClasses
///
/// @author Asher Kaboth
/// @author Richard Calland
class SampleHandlerInterface
{
 public:
   /// @brief The main constructor
   SampleHandlerInterface();
  /// @brief destructor
  virtual ~SampleHandlerInterface();

  /// @brief returns total number of samples
  virtual M3::int_t GetNSamples(){ return nSamples; };
  /// @brief Get fancy title for specified samples
  /// @param iSample Sample enumerator
  virtual std::string GetSampleTitle(const int iSample) const = 0;
  /// @brief Get name for Sample Handler
  virtual std::string GetName() const = 0;
  /// @brief Get likelihood (-logL) for a single sample
  /// @param iSample Sample enumerator
  virtual double GetSampleLikelihood(const int iSample) const = 0;
  /// @brief Allow to clean not used memory before fit starts
  virtual void CleanMemoryBeforeFit() = 0;
  /// @brief Store additional info in a chain
  /// @param Dir directory to which we save additional info
  virtual void SaveAdditionalInfo([[maybe_unused]] TDirectory* Dir) {};
  /// @brief Return pointer to MaCh3 modes
  MaCh3Modes* GetMaCh3Modes() const { return Modes.get(); }
  /// @brief main routine modifying MC prediction based on proposed parameter values
  virtual void Reweight()=0;
  /// @brief Return likelihood (-logL) for all samples
  virtual double GetLikelihood() const = 0;

  /// @brief Helper function to print rates for the samples with LLH
  /// @param DataOnly whether to print data only rates
  virtual void PrintRates(const bool DataOnly = false) = 0;

  /// @brief Return total number of events
  unsigned int GetNEvents() const {return nEvents;}
  /// @brief Get number of oscillation channels for a single sample
  /// @param iSample Sample enumerator
  virtual int GetNOscChannels(const int iSample) const = 0;

  /// @brief Return Kinematic Variable name for specified sample and dimension for example "Reconstructed_Neutrino_Energy"
  /// @param iSample Sample index
  /// @param Dimension Dimension index
  virtual std::string GetKinVarName(const int iSample, const int Dimension) const = 0;

  /// @brief Get Data histogram
  /// @param Sample Sample enumerator
  virtual const TH1* GetDataHist(const int Sample) = 0;
  /// @brief Get MC histogram
  /// @param Sample Sample enumerator
  virtual const TH1* GetMCHist(const int Sample) = 0;
  /// @brief Get W2 histogram
  /// @param Sample Sample enumerator
  virtual const TH1* GetW2Hist(const int Sample) = 0;

  /// @brief DB Get what dimensionality binning for given sample has
  /// @param Sample Number of sample
  virtual int GetNDim(const int Sample) const = 0;
  /// @brief Get the flavour name for a given sample and oscillation channel.
  /// @param iSample Index of the sample.
  /// @param iChannel Index of the oscillation channel within the sample.
  virtual std::string GetFlavourName(const int iSample, const int iChannel) const = 0;

  /// @brief Return the binning used to draw a kinematic parameter
  /// @param iSample Index of the sample.
  /// @param KinematicParameter name of variable
  virtual std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string &KinematicParameter) const = 0;

  /// @brief Build a 1D histogram for a given variable, optionally filtered by mode and channel.
  /// @param iSample Index of the sample.
  /// @param ProjectionVar_Str Name of the variable to project onto.
  /// @param kModeToFill Interaction mode to select (-1 means all modes).
  /// @param kChannelToFill Oscillation channel to select (-1 means all channels).
  /// @param WeightStyle Weighting scheme (e.g. 0 = nominal weights, 1 = unit weights).
  virtual std::unique_ptr<TH1> Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str,
                                                            const int kModeToFill = -1, const int kChannelToFill = -1,
                                                            const int WeightStyle = 0) = 0;
  /// @brief Build a 2D histogram for given variables, optionally filtered by mode and channel.
  /// @param iSample Index of the sample.
  /// @param ProjectionVar_StrX Name of the variable for the X axis.
  /// @param ProjectionVar_StrY Name of the variable for the Y axis.
  /// @param kModeToFill Interaction mode to select (-1 means all modes).
  /// @param kChannelToFill Oscillation channel to select (-1 means all channels).
  /// @param WeightStyle Weighting scheme (e.g. 0 = nominal weights, 1 = unit weights).
  virtual std::unique_ptr<TH2> Get2DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_StrX,
                                                            const std::string& ProjectionVar_StrY, int kModeToFill = -1,
                                                            const int kChannelToFill = -1, const int WeightStyle = 0) = 0;
  /// @brief Return 1D projection of MC into given 1D variable (doesn't have to be variable used in the fit)
  /// @param iSample Index of the sample.
  /// @param ProjectionVar name of variable
  /// @param EventSelectionVec Vector of additional cuts like cut on interaction mode
  /// @param WeightStyle Alow to modify weight for example if equal to 1 all weights are set to 1
  /// @param SubEventSelectionVec Vector of additional cuts for sub event (particle, ring etc.)
  virtual std::unique_ptr<TH1> Get1DVarHist(const int iSample, const std::string &ProjectionVar,
                                            const std::vector<KinematicCut> &EventSelectionVec = {}, int WeightStyle = 0,
                                            const std::vector<KinematicCut> &SubEventSelectionVec = {}) = 0;
  /// @brief Build a 2D projection of MC events into specified variables.
  /// @param iSample Index of the sample.
  /// @param ProjectionVarX Name of the variable for the X axis.
  /// @param ProjectionVarY Name of the variable for the Y axis.
  /// @param EventSelectionVec Vector of event-level selection cuts.
  /// @param WeightStyle Weighting scheme (e.g. 0 = nominal weights, 1 = unit weights).
  /// @param SubEventSelectionVec Vector of sub-event selection cuts.
  virtual std::unique_ptr<TH2> Get2DVarHist(const int iSample, const std::string& ProjectionVarX,
                                            const std::string& ProjectionVarY,
                                            const std::vector< KinematicCut >& EventSelectionVec = {},
                                            const int WeightStyle = 0, const std::vector< KinematicCut >& SubEventSelectionVec = {}) = 0;

  /// @brief Calculate test statistic for a single bin using Poisson
  /// @param data is data
  /// @param mc is mc
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
  double GetTestStatLLH(const double data, const double mc, const double w2) const;

  /// @brief Set the test statistic to be used when calculating the binned likelihoods
  /// @param testStat The test statistic to use.
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
