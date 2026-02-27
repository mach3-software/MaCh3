#pragma once

#include "Fitters/FitterBase.h"

class ParameterHandlerGeneric;

// ***************************
/// @brief KS: Summary of sample info to be used by
struct PredictiveSample {
// ***************************
  // Name of sample
  std::string Name;
  /// Pointer to SampleHandler
  const SampleHandlerBase* SamHandler;
  /// Local SampleId in SampleHandler
  int LocalId;
  /// Sample Dimension
  int Dimenstion;
};

/// @brief Implementation of Prior/Posterior Predictive and Bayesian p-Value calculations following the approach described in @cite gelman1996posterior.
/// @details For more information, visit the <a href="PosteriorPredictive.html">Posterior Predictive page</a>.
///
/// @author Asher Kaboth
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Yashwanth S Prabhu
/// @author Kamil Skwarczynski
/// @author Patrick Dunne
/// @author Clarence Wret

/// @todo add BIC, DIC, WAIC
/// @todo add Rate $p$-value
/// @todo add plots by mode
/// @todo add rate plate and rate error reduction
/// @todo unify code with SampleSummary
class PredictiveThrower : public FitterBase {
 public:
   /// @brief Constructor
   /// @param fitMan A pointer to a manager object, which will handle all settings.
  PredictiveThrower(Manager * const fitMan);
  /// @brief Destructor
  virtual ~PredictiveThrower();

  /// @brief Produce toys by throwing from MCMC
  void ProduceToys();

  /// @brief Main routine responsible for producing posterior predictive distributions and $p$-value
  void RunPredictiveAnalysis();

  /// @brief This is not used in this class
  void RunMCMC() override {
    MACH3LOG_ERROR("{} is not supported in {}", __func__, GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  };

 private:
  /// @brief This set some params to prior value this way you can evaluate errors from subset of errors
  void SetParamters();

  /// @brief Setup useful variables etc before stating toy generation
  void SetupToyGeneration();

  /// @brief Load existing toys
  bool LoadToys();
  /// @brief Save histograms for a single MCMC Throw/Toy
  void WriteToy(TDirectory* ToyDirectory, TDirectory* Toy_1DDirectory, TDirectory* Toy_2DDirectory, const int iToy);
  /// @brief Setup sample information
  void SetupSampleInformation();

  /// @brief Get Fancy parameters stored in mcmc chains for passed ParameterHandler
  std::vector<std::string> GetStoredFancyName(ParameterHandlerBase* Systematics) const;

  /// @brief Produce posterior predictive distribution
  std::vector<std::unique_ptr<TH1>> MakePredictive(const std::vector<std::vector<std::unique_ptr<TH1>>>& Toys,
                                                   const std::vector<TDirectory*>& Director,
                                                   const std::string& suffix,
                                                   const bool DebugHistograms);
  /// @brief Load 1D projections and later produce violin plots for each
  void Study1DProjections(const std::vector<TDirectory*>& SampleDirectories) const;
  /// @brief Produce Violin style spectra
  void ProduceSpectra(const std::vector<std::vector<std::vector<std::unique_ptr<TH1D>>>>& Toys,
                      const std::vector<TDirectory*>& Director,
                      const std::string suffix) const;

  /// @brief Make Poisson fluctuation of TH1D hist
  void MakeFluctuatedHistogram(TH1* FluctHist, TH1* PolyHist, const int Dim);
  /// @brief Calculate Posterior Predictive $p$-value
  void PosteriorPredictivepValue(const std::vector<std::unique_ptr<TH1>>& PostPred_mc,
                                 const std::vector<TDirectory*>& SampleDir);


  /// @brief Helper functions to calculate likelihoods using TH1
  /// @param Data histogram with data distribution for a single sample
  /// @param MC histogram with MC distribution for a single sample
  /// @param W2 histogram with W2 distribution for a single sample
  double GetLLH(const std::unique_ptr<TH1>& DatHist,
                                   const std::unique_ptr<TH1>& MCHist,
                                   const std::unique_ptr<TH1>& W2Hist,
                                   const SampleHandlerBase* SampleHandler);

  /// @brief Produce Chi2 plot for a single sample based on which $p$-value is calculated
  void MakeChi2Plots(const std::vector<std::vector<double>>& Chi2_x,
                     const std::string& Chi2_x_title,
                     const std::vector<std::vector<double>>& Chi2_y,
                     const std::string& Chi2_y_title,
                     const std::vector<TDirectory*>& SampleDir,
                     const std::string Title);


  /// @brief Construct a human-readable label describing a specific analysis bin.
  /// @param hist Histogram providing the binning definition.
  /// @param uniform Flag indicating whether the histogram uses regular axis
  ///        binning (TH1/TH2) or irregular polygonal binning (e.g. TH2Poly).
  /// @param Dim Dimensionality of the original distribution.
  /// @param bins Vector of per-dimension bin indices in analysis coordinates.
  std::string GetBinName(TH1* hist,
                         const bool uniform,
                         const int Dim,
                         const std::vector<int>& bins) const;
  /// @brief Create per-bin posterior histograms for a given sample.
  ///
  /// For each analysis bin of the input histogram, this function allocates a new
  /// 1D histogram intended to accumulate the distribution of predicted event
  /// counts (e.g. across throws, toys, or posterior evaluations).
  ///
  /// The number of output histograms therefore equals the number of physical bins:
  /// - TH1  → N histograms
  /// - TH2  → Nx × Ny histograms
  /// - TH2Poly → one histogram per polygon bin
  ///
  /// @param hist Input histogram defining the bin structure for this sample.
  /// @param SampleId Index identifying the sample in SampleInfo.
  /// @param Dim Dimensionality of the original distribution.
  /// @param suffix String appended to histogram names (e.g. to distinguish stages).
  std::vector<std::unique_ptr<TH1D>> PerBinHistogram(TH1* hist,
                                                     const int SampleId,
                                                     const int Dim,
                                                     const std::string& suffix) const;

  /// @brief Evaluate prior/post predictive distribution for beta parameters (used for evaluating impact MC statistical uncertainty)
  void StudyBetaParameters(TDirectory* PredictiveDir);
  /// @brief Make the 1D Event Rate Hist
  void MakeCutEventRate(TH1D *Histogram, const double DataRate) const;
  /// @brief Produce distribution of number of events for each sample
  void RateAnalysis(const std::vector<std::vector<std::unique_ptr<TH1>>>& Toys,
                                       const std::vector<TDirectory*>& SampleDirectories) const;

  /// KS: Use Full LLH or only sample contribution based on discussion with Asher we almost always only want the sample likelihood
  bool FullLLH;
  /// KS: Count total number of model parameters which can be used for stuff like BIC
  int NModelParams;
  /// Whether it is Prior or Posterior predictive
  bool Is_PriorPredictive;

  /// Number of toys we are generating analysing
  int TotalNumberOfSamples;

  /// Handy struct for all sample info
  std::vector<PredictiveSample> SampleInfo;

  /// Number of toys we are generating analysing
  int Ntoys;
  /// KS: Names of parameter groups that will not be varied
  std::vector<std::string> ParameterGroupsNotVaried;
  /// KS: Index of parameters that will be varied
  std::unordered_set<int> ParameterOnlyToVary;

  /// Pointer to El Generico
  ParameterHandlerGeneric* ModelSystematic;

  /// Vector of Data histograms
  std::vector<std::unique_ptr<TH1>> Data_Hist;
  /// Vector of MC histograms
  std::vector<std::unique_ptr<TH1>> MC_Nom_Hist;
  /// Vector of W2 histograms
  std::vector<std::unique_ptr<TH1>> W2_Nom_Hist;

  /// Vector of MC histograms per sample and toy experiment.
  /// Indexed as [sample][toy].
  std::vector<std::vector<std::unique_ptr<TH1>>> MC_Hist_Toy;
  /// Vector of W² histograms per sample and toy experiment.
  /// Indexed as [sample][toy]
  std::vector<std::vector<std::unique_ptr<TH1>>> W2_Hist_Toy;

  /// Reweighting factors applied for each toy, by default 1
  std::vector<double> ReweightWeight;
  /// Penalty term values for each toy by default 0
  std::vector<double> PenaltyTerm;

  /// KS: We have two methods for Poissonian fluctuation
  bool StandardFluctuation;
};

