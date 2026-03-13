#pragma once

//MaCh3 includes
#include "Fitters/LikelihoodFit.h"

// ROOT includes
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

/// @brief Implementation of Minuit fitting algorithm
/// @cite James:2004xla
/// @author Kamil Skwarczynski
///
/// @ingroup FittingAlgorithms
class MinuitFit : public LikelihoodFit {
 public:
  /// @brief Constructor
  MinuitFit(Manager * const fitMan);
  /// @brief Destructor
  virtual ~MinuitFit();

  /// @brief Actual implementation of Minuit Fit algorithm
  void RunMCMC() override;
 private:
  /// Pointer to minimizer, which most often is Minuit
  std::unique_ptr<ROOT::Math::Minimizer> minuit;
};

