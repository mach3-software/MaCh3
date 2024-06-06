#pragma once

// MaCh3 includes
#include "covariance/covarianceBase.h"

/// @brief Class responsible for handling of neutrino oscillation  parameters.
class covarianceOsc : public covarianceBase
{
 public:
  /// @brief Constructor
  covarianceOsc(const char* name, const char *file);
  /// @brief Destructor
  virtual ~covarianceOsc();
  /// @brief Same implementation as in base class but accounts for applying reactor prior
  /// @warning May become deprecated and use only base class instead
  double GetLikelihood() override;
  /// @brief Checks if parameter didn't go outside of physical bounds
  inline int CheckBounds() override;
  /// @brief Retrieves the proposed parameters.
  double *getPropPars();
  /// @brief Propose MCMC step, including mass fliping
  void proposeStep() override;
  /// @brief Sets whether to flip delta M23.
  void setFlipDeltaM23(bool flip){flipdelM = flip;}
  /// @brief Sets whether to flip beta.
  void setFlipBeta(bool flip){flipBeta = flip;}
  /// @brief Sets whether to use reactor prior.
  /// @warning May become deprecated
  void useReactorPrior(bool reactor){reactorPrior = reactor;};
  /// @brief Sets extra branches for processing.
  void setExtraBranches(TTree &tree);

  /// @brief Get baseline
  inline double GetPathLength() { return L;}
  /// @brief Get density
  inline double GetDensity() { return density;}
  /// @brief Check if we are performing beta study or not
  inline bool GetPerformBetaStudy(){return PerformBetaStudy;}
  /// @brief KS: Print all useful information's after initialization
  void Print();
  /// @brief KS: Currently prob3++/probgpu requires particular order so we need to check this is the case
  void CheckOrderOfParams();

 protected:
    /// Value of baseline
    double L;
    /// Value of density
    double density;
    /// Do we flip DeltaM23 or not
    bool flipdelM;
    /// Bool whether we apply additional reactor prior or not
    bool reactorPrior;
    /// Do we flip Beta or not
    bool flipBeta;
    /// Value of all osc parameters
    double *oscpars1;
    /// Bool whether we are performing beta study or not
    bool PerformBetaStudy;

    /// There is special treatment for delta CP, therefore store enum keeping track when to apply special treatment
    int kDeltaCP;
    /// There is special treatment for DeltaM23, therefore store enum keeping track when to apply special treatment
    int kDeltaM23;
    /// Enum for special treatment of sin(theta23).
    int kSinTheta23;
    /// Enum for special treatment of beta.
    int kBeta;
};

