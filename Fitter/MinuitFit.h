#pragma once

//MaCh3 includes
#include "Fitter/LikelihoodFit.h"

// ROOT includes
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

/// @brief Implementation of Minuit fitting algorithm
/// @cite James:2004xla
/// @author Kamil Skwarczynski
class MinuitFit : public LikelihoodFit {
 public:
  /// @brief Constructor
  MinuitFit(manager * const fitMan);
  /// @brief Destructor
  virtual ~MinuitFit();

  /// @brief Actual implementation of Minuit Fit algorithm
  void runMCMC() override;

  /// @brief Get name of class
  inline std::string GetName()const {return "MinuitFit";};

 private:
  /// Pointer to minimizer, which most often is Minuit
  ROOT::Math::Minimizer* minuit;
};

