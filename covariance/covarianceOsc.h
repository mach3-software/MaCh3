#pragma once

// MaCh3 includes
#include "covariance/covarianceBase.h"

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

  inline double GetPathLength() { return L;}
  inline double GetDensity() { return density;}
  inline bool GetPerformBetaStudy(){return PerformBetaStudy;}
  /// KS: Print all useful information's after initialization
  void Print();
  /// KS: Currently prob3++/probgpu requires particular order so we need to check this is the case
  void CheckOrderOfParams();

 protected:
    double L;
    double density;
    bool flipdelM;
    bool reactorPrior;
    bool flipBeta;
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

