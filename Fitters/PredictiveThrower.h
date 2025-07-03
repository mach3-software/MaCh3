#pragma once

#include "Fitters/FitterBase.h"

class ParameterHandlerGeneric;

/// @brief Implementation of Prior/Posterior Predictive and Bayesian p-Value calculations following the approach described in @cite gelman1996posterior.
/// @author Asher Kaboth
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Yashwanth S Prabhu
/// @author Kamil Skwarczynski
/// @author Patrick Dunne
/// @author Clarence Wret
class PredictiveThrower : public FitterBase {
 public:
   /// @brief Constructor
   /// @param fitMan A pointer to a manager object, which will handle all settings.
  PredictiveThrower(manager * const fitMan);
  /// @brief Destructor
  virtual ~PredictiveThrower();

  /// @brief Produce toys by throwing from MCMC
  void ProduceToys();

  void RunPredictiveAnalysis();

  /// @brief This is not used in this class
  void RunMCMC() override {
    MACH3LOG_ERROR("{} is not supported in {}", __func__, GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  };

  /// @brief Get name of class
  inline std::string GetName()const {return "PredictiveThrower";};
 private:
  /// @brief This set some params to prior value this way you can evaluate errors from subset of errors
  void SetParamters();

  /// @brief Setup useful variables etc before stating toy generation
  void SetupToyGeneration();

  /// TODO!!!!
  std::unique_ptr<TH1D> MakeSpectra(const std::unique_ptr<TH2D>& Spectra,
                                    const std::string& Sample_Name,
                                    const std::string& suffix);
  /// TODO!!!!
  std::vector<std::unique_ptr<TH2D>> ProduceSpectra(const std::vector<std::vector<std::unique_ptr<TH1D>>>& Toys,
                                                    const std::string suffix);

  /// TODO!!!!
  void PosteriorPredictivepValue(const std::vector<std::unique_ptr<TH1D>>& PostPred_mc,
                                 //const std::vector<std::unique_ptr<TH1D>>& PostPred_w2,
                                 const std::vector<TDirectory*>& SampleDir);


  /// TODO!!!!
  double GetLLH(const std::unique_ptr<TH1D>& DatHist,
                                   const std::unique_ptr<TH1D>& MCHist,
                                   const std::unique_ptr<TH1D>& W2Hist,
                                   SampleHandlerBase* SampleHandler);

  /// TODO!!!!
  void MakeChi2Plots(const std::vector<std::vector<double>>& Chi2_x,
                     const std::string& Chi2_x_title,
                     const std::vector<std::vector<double>>& Chi2_y,
                     const std::string& Chi2_y_title,
                     const std::vector<TDirectory*>& SampleDir,
                     const std::string Tittle);

  /// KS: Use Full LLH or only sample contribution based on discussion with Asher we almost always only want the sample likelihood
  bool FullLLH;

  /// KS: Count total number of model parameters which can be used for stuff like BIC
  int NModelParams;

  /// Whether it is Prior or Posterior predictive
  bool Is_PriorPredictive;

  /// Number of toys we are generating analysing
  int Ntoys;
  /// Number of toys we are generating analysing
  int TotalNumberOfSamples;
  std::vector<std::string> SampleNames;
  /// KS: Names of parameter groups that will not be varied
  std::vector<std::string> ParameterGroupsNotVaried;

  /// Pointer to El Generico
  ParameterHandlerGeneric* ModelSystematic;


  std::unique_ptr<TH1D> MakePredictive(const std::vector<std::unique_ptr<TH1D>>& Toys,
                                                          const std::string& Sample_Name,
                                                          const std::string& suffix);

  std::vector<std::unique_ptr<TH1D>> Data_Hist;


  /// TODO!!!!
  //[samples][toy]
  std::vector<std::vector<std::unique_ptr<TH1D>>> MC_Hist_Toy;
  /// TODO!!!!
  std::vector<std::vector<std::unique_ptr<TH1D>>> W2_Hist_Toy;

  /// TODO!!!!
  std::vector<double> ReweightWeight;
  /// TODO!!!!
  std::vector<double> PenaltyTerm;
};

