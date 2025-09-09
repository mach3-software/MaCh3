#pragma once

#include "Fitters/FitterBase.h"

class ParameterHandlerGeneric;

/// @brief Implementation of Prior/Posterior Predictive and Bayesian p-Value calculations following the approach described in @cite gelman1996posterior.
/// @details For more information, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/10.-Posterior-Predictive,-p%E2%80%90value-etc.).
/// @author Asher Kaboth
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Yashwanth S Prabhu
/// @author Kamil Skwarczynski
/// @author Patrick Dunne
/// @author Clarence Wret

/// @todo add BIC, DIC, WAIC
/// @todo add ability yo make projection for Get1DDiscVar
/// @todo add ability for TH2D
/// @todo speed improvements
/// @todo add Rate $p$-value
/// @todo unify code with SampleSummary
class PredictiveThrower : public FitterBase {
 public:
   /// @brief Constructor
   /// @param fitMan A pointer to a manager object, which will handle all settings.
  PredictiveThrower(manager * const fitMan);
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

  /// @brief Get name of class
  inline std::string GetName() const override {return "PredictiveThrower";};
 private:
  /// @brief This set some params to prior value this way you can evaluate errors from subset of errors
  void SetParamters();

  /// @brief Setup useful variables etc before stating toy generation
  void SetupToyGeneration();

  /// @brief Load existing toys
  bool LoadToys();

  /// @brief Setup sample information
  void SetupSampleInformation();

  /// @brief Produce posterior predictive distribution
  std::unique_ptr<TH1D> MakePredictive(const std::vector<std::unique_ptr<TH1D>>& Toys,
                                       const std::string& Sample_Name,
                                       const std::string& suffix,
                                       const bool DebugHistograms);

  /// @brief Produce Violin style spectra
  std::vector<std::unique_ptr<TH2D>> ProduceSpectra(const std::vector<std::vector<std::unique_ptr<TH1D>>>& Toys,
                                                    const std::string suffix);

  /// @brief Calculate Posterior Predictive $p$-value
  void PosteriorPredictivepValue(const std::vector<std::unique_ptr<TH1D>>& PostPred_mc,
                                 //const std::vector<std::unique_ptr<TH1D>>& PostPred_w2,
                                 const std::vector<TDirectory*>& SampleDir);


  /// @brief Helper functions to calculate likelihoods using TH1D
  /// @param Data histogram with data distribution for a single sample
  /// @param MC histogram with MC distribution for a single sample
  /// @param W2 histogram with W2 distribution for a single sample
  double GetLLH(const std::unique_ptr<TH1D>& DatHist,
                                   const std::unique_ptr<TH1D>& MCHist,
                                   const std::unique_ptr<TH1D>& W2Hist,
                                   SampleHandlerBase* SampleHandler);

  /// @brief Produce Chi2 plot for a single sample based on which $p$-value is calculated
  void MakeChi2Plots(const std::vector<std::vector<double>>& Chi2_x,
                     const std::string& Chi2_x_title,
                     const std::vector<std::vector<double>>& Chi2_y,
                     const std::string& Chi2_y_title,
                     const std::vector<TDirectory*>& SampleDir,
                     const std::string Title);

  /// KS: Use Full LLH or only sample contribution based on discussion with Asher we almost always only want the sample likelihood
  bool FullLLH;
  /// KS: Count total number of model parameters which can be used for stuff like BIC
  int NModelParams;
  /// Whether it is Prior or Posterior predictive
  bool Is_PriorPredictive;

  /// Number of toys we are generating analysing
  int TotalNumberOfSamples;
  /// Name of a single sample
  std::vector<std::string> SampleNames;
  /// Maps if sample with given SampleHandler, useful if we have more than one sample in single object
  std::vector<int> SampleObjectMap;

  /// Number of toys we are generating analysing
  int Ntoys;
  /// KS: Names of parameter groups that will not be varied
  std::vector<std::string> ParameterGroupsNotVaried;
  /// KS: Index of parameters groups that will be varied
  std::unordered_set<int> ParameterOnlyToVary;

  /// Pointer to El Generico
  ParameterHandlerGeneric* ModelSystematic;

  /// Vector of Data histograms
  std::vector<std::unique_ptr<TH1D>> Data_Hist;
  /// Vector of MC histograms
  std::vector<std::unique_ptr<TH1D>> MC_Nom_Hist;
  /// Vector of W2 histograms
  std::vector<std::unique_ptr<TH1D>> W2_Nom_Hist;

  /// Vector of MC histograms per sample and toy experiment.
  /// Indexed as [sample][toy].
  std::vector<std::vector<std::unique_ptr<TH1D>>> MC_Hist_Toy;
  /// Vector of W² histograms per sample and toy experiment.
  /// Indexed as [sample][toy]
  std::vector<std::vector<std::unique_ptr<TH1D>>> W2_Hist_Toy;

  /// Reweighting factors applied for each toy, by default 1
  std::vector<double> ReweightWeight;
  /// Penalty term values for each toy by default 0
  std::vector<double> PenaltyTerm;
};

