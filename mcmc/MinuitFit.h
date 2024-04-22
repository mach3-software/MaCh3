#pragma once

#include "LikelihoodFit.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

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

