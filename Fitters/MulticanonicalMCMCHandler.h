#pragma once
#include <ostream>
#include <string>
#include <vector>

#include "Parameters/ParameterHandlerBase.h"
#include "Manager/Manager.h"

_MaCh3_Safe_Include_Start_ //{
#include "TSpline.h"
_MaCh3_Safe_Include_End_ //}


namespace M3 {
  enum class BiasFunction {
    kGaussian,
    kVonMises,
    kGeneralisedGaussian
  };


  /// Normalisation constant for the n=2 umbrella sampling Gaussian.
  ///
  /// The factor comes from:
  ///   2n / Gamma(1/(2n))
  /// with n = 2:
  ///   4 / Gamma(1/4) = 0.906402477055
  ///
  constexpr double UmbrellaGaussianNormFactor = 0.906402477055;
}

/// @brief Helper class for configuring and evaluating multicanonical umbrella weights.
/// @details this method is used to bias the MCMC sampling of chains to particular regions of the parameter space
/// The method relies on a bias function defined via a yaml configuration file in experimental specific code
/// With a bias function defined chains run as normal and are combined after the fact using the information
/// stored in the chains saved yaml configuration file.
///
/// @author David Riley
/// @cite Matthews2018umbrella
class MulticanonicalMCMCHandler {
public:
  /// @brief Constructor
  MulticanonicalMCMCHandler();
  /// @brief Destructor
  virtual ~MulticanonicalMCMCHandler();

#ifdef DEBUG
  /// @brief Set the optional debug output stream.
  void setDebugStream(std::ostream* os, bool enabled);
#endif

  /// @brief Locate the systematic object which contains the parameters used by the handler. Stores the parameter indices.
  void FindOscCovParams(const std::vector<ParameterHandlerBase*>& systematics);

  /// @brief Adjust the parameter-of-interest's step scale according to the width of the umbrella bias function.
  /// @details This is used to adjust the step scale of the parameter to keep the size of the step scale
  /// relative to the new biased space the chain sees the same. Current assumes that your are sampling delta cp (ie cirular on 2pi)
  /// @todo grab the step scale from the objejct directly and update based on the radtio to the umbrella width, removing assumption that its dcp 
  void AdjustUmbrellaStepScale(const std::vector<ParameterHandlerBase*>& systematics);

  /// @brief Read multicanonical configuration from the yaml via manager and initialise handler state.
  void InitializeMulticanonicalHandlerConfig(Manager* fitMan, std::vector<ParameterHandlerBase*>& systematics);

  /// @brief Initialise the starting values used by the multicanonical parameter handling.
  /// @details This moves the starting point of the chain to the centre of the umbrella bias function. Only affects the parameter of interest. 
  void InitializeMulticanonicalParams(std::vector<ParameterHandlerBase*>& systematics);

  /// @brief Compute the multicanonical penalty for the configured bias mode.
  /// @details This wraps the various bias function implementations and returns the appropriate penalty
  /// @return Returns the log likelihood penalty for the selected bias function at the proposed parameter values.
  double GetMulticanonicalWeight(double deltacp, double delm23_value);

  /// @brief Compute the multicanonical penalty using a spline.
  double GetMulticanonicalWeightSpline(double deltacp, double delm23_value);

  /// @brief Compute a Gaussian multicanonical penalty.
  double GetMulticanonicalWeightGaussian(double deltacp);

  /// @brief Compute a triple-Gaussian multicanonical penalty.
  /// @details This was mostsly for testing the first implementation, probably not used anymore.
  double GetMulticanonicalWeightTripleGaussian(double deltacp);

  /// @brief Compute a von Mises multicanonical penalty.
  /// @details This is the circular analogue of a gaussian, meaning it handles the wrapping of the parameter space at 2pi automatically
  /// this avoids the need to calculate the weight multiple times for each parameter to ensure it recieves a bias even when it jumps the boundary
  /// the normalisation grows extremely quickly, widths of less that 0.05 should not be used until properly tested.
  double GetMulticanonicalWeightVonMises(double deltacp);

  /// @brief Compute a generalised-Gaussian multicanonical penalty.
  /// @details the generalised gaussian is like a gaussian with an extra factor of n on the exponent of the gaussian. This allows for a stronger 
  /// constraint in the tails of the distribution, more like a tophat function. The normalisation is a little ugly and uses gamma functions, so for now it is hardcoded to n=2.
  /// This is the current recommended configuration for umbrella sampling.
  /// @todo implememnt the wrapping with circular distance function to avoid the need to calculate the weight multiple times for each parameter to ensure it recieves a bias even when it jumps the boundary
  /// @todo implement the normalisation for n != 2, and allow n to be set in the yaml configuration file
  double GetMulticanonicalWeightGenGaussian(double deltacp);

  /// @brief Compute the circular distance between two angles.
  /// @details This is used to calculate the distance between two angles in a circular space, such as delta_cp. It returns the shortest distance between the two angles, taking into account the wrapping at 2pi.
  double circularDistance(double x, double mean);

  /// @brief Wraps the generalised gaussian function for a given x, mean, and width. Required to handle the wrappipng of the parameter space at 2pi.
  double generalisedGaussian2(double x, double mean, double width);

  /// @brief Index of the oscillation-covariance systematic in the current fit.
  int oscCovVar;
  /// @brief Parameter index used for the multicanonical delta_cp weight.
  int multicanonicalVar;
  /// @brief Parameter index used for the multicanonical delm2_23 weight.
  int multicanonicalVar_dm23;

  /// @brief Selected bias function for multicanonical weights.
  M3::BiasFunction umbrellaBiasFunction;

  /// @brief Configured bias function name for logging.
  std::string umbrellaBiasFunctionName;

  /// @brief Toggle for spline-based multicanonical weights.
  bool multicanonicalSpline;

protected:
  /// @brief Global scale factor applied to the multicanonical penalty. 1 is full strength, 0 is no penalty.
  double multicanonicalBeta;

  /// @brief delta_cp value used during proposal evaluation.
  double delta_cp_value;
  /// @brief delm2_23 value used during proposal evaluation.
  double delm23_value;

  /// @brief Spline for the IO branch, if spline mode is enabled.
  TSpline3* dcp_spline_IO;
  /// @brief Spline for the NO branch, if spline mode is enabled.
  TSpline3* dcp_spline_NO;

  /// @brief Umbrella centre used for the current configuration.
  double umbrellaMean;
  /// @brief Umbrella width used for the current configuration.
  double umbrellaWidth;
  /// @brief Number of total umbrellas used for the current configuration. 
  /// @details Currently only used when auto-overlap mode is enabled. This is not fully tested.
  int umbrellaNumber;
  /// @brief Toggle for deriving umbrella widths from a desired # of sigma overlaps between umbrellas.
  bool umbrellaOverlapMode;
  /// @brief Requested overlap for evenly spaced umbrellas.
  double umbrellaSigmaOverlap;
  /// @brief Toggle for rescaling the step size based on the umbrella width.
  bool umbrellaAdjustStepScale;
  /// @brief Additional scale factor applied when rescaling the step size. This is for fine tuning
  double umbrellaStepScaleFactor;
  /// @brief Optional flip-window control. 
  /// @details briefly tested feature to mirrow windows about their centres. Not tested but may be useful in the future. Currently disabled and not used.
  bool flipWindow;

  /// @brief Von Mises kappa parameter. Analogue of sigma for a gaussian
  double vonMises_kappa;
  /// @brief Cached I0(kappa) value for the von Mises form.
  double vonMises_I0_kappa;

#ifdef DEBUG
  std::ostream* debugStream = nullptr;
  bool debugEnabled = false;
#endif
};
