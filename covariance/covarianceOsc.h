#pragma once

// MaCh3 includes
#include "covariance/covarianceBase.h"

/// @brief Class responsible for handling of neutrino oscillation  parameters.
class covarianceOsc : public covarianceBase
{
 public:
  /// @brief Constructor
  covarianceOsc(const std::vector<std::string>& FileNames, const char *name = "osc_cov", double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
  /// @brief Destructor
  virtual ~covarianceOsc();
  /// @brief Checks if parameter didn't go outside of physical bounds
  inline int CheckBounds() override;
  /// @brief Retrieves the proposed parameters.
  double *getPropPars();
  /// @brief Propose MCMC step, including mass flipping
  void proposeStep() override;
  /// @brief Sets whether to flip delta M23.
  void setFlipDeltaM23(bool flip){flipdelM = flip;}

  /// @brief Get baseline
  inline double GetPathLength() { return _fPropVal[kBaseline];}
  /// @brief Get density
  inline double GetDensity() { return _fPropVal[kDensity];}

  /// @brief KS: Print all useful information's after initialization
  void Print();

 protected:

    /// Do we flip DeltaM23 or not
    bool flipdelM;

    /// There is special treatment for delta CP, therefore store enum keeping track when to apply special treatment
    int kDeltaCP;
    /// There is special treatment for DeltaM23, therefore store enum keeping track when to apply special treatment
    int kDeltaM23;
    /// Enum for special treatment of sin(theta23).
    int kSinTheta23;
    /// Enum for special treatment of of baseline
    int kBaseline;
    /// Enum for special treatment of of density
    int kDensity;
};
