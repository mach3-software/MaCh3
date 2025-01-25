#pragma once

// MaCh3 includes
#include "covariance/covarianceBase.h"

/// @brief Class responsible for handling of neutrino oscillation  parameters.
/// @author Richard Calland
/// @author Asher Kaboth
class covarianceOsc : public covarianceBase
{
 public:
  /// @brief Constructor
  covarianceOsc(const std::vector<std::string>& FileNames, std::string name = "osc_cov", double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
  /// @brief Destructor
  virtual ~covarianceOsc();
  /// @brief Propose MCMC step, including mass flipping
  void proposeStep() override;
  /// @brief Sets whether to flip delta M23.
  void setFlipDeltaM23(bool flip){flipdelM = flip;}
  /// @brief Get pointers to Osc params from detId
  std::vector<const double*> GetOscParsFromDetID(const std::string& DetID);
  /// @brief KS: Print all useful information's after initialization
  void Print();

 protected:
    /// @brief HW :: This method is a tad hacky but modular arithmetic gives me a headache.
    /// @author Henry Wallace
    void CircularPrior(const int i, const double LowBound, const double UpBound);
    /// Do we flip DeltaM23 or not
    bool flipdelM;
    /// There is special treatment for delta CP, therefore store enum keeping track when to apply special treatment
    int kDeltaCP;
    /// There is special treatment for DeltaM23, therefore store enum keeping track when to apply special treatment
    int kDeltaM23;
    /// Enum for special treatment of sin(theta23).
    int kSinTheta23;
};
