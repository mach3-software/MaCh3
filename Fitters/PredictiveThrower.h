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

  /// KS: Use Full LLH or only sample contribution based on discussion with Asher we almost always only want the sample likelihood
  bool FullLLH;

  //KS: Count total number of model parameters which can be used for stuff like BIC
  int NModelParams;

  ParameterHandlerGeneric* ModelSystematic;


};

