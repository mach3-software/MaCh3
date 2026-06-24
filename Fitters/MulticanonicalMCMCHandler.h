#pragma once

#include "TSpline.h"
#include "Parameters/ParameterHandlerBase.h"
#include "Manager/Manager.h"
#include <ostream>
#include <string>
#include <vector>

class MulticanonicalMCMCHandler {
public:
  enum class BiasFunction { Gaussian, VonMises, GeneralisedGaussian };

  /// @brief Constructor
  MulticanonicalMCMCHandler();
  /// @brief Destructor
  virtual ~MulticanonicalMCMCHandler();

#ifdef DEBUG
  void setDebugStream(std::ostream* os, bool enabled);
#endif

  void FindOscCovParams(const std::vector<ParameterHandlerBase*>& systematics);
  void AdjustUmbrellaStepScale(const std::vector<ParameterHandlerBase*>& systematics);

  void InitializeMulticanonicalHandlerConfig(Manager* fitMan, std::vector<ParameterHandlerBase*>& systematics);

  void InitializeMulticanonicalParams(std::vector<ParameterHandlerBase*>& systematics);

  /// getters for the multicanonical weights these will take in delta values and return a llh penalty based on the chosen multicanonical method
  double GetMulticanonicalWeight(double deltacp, double delm23_value);

  double GetMulticanonicalWeightSpline(double deltacp, double delm23_value);

  double GetMulticanonicalWeightGaussian(double deltacp);

  double GetMulticanonicalWeightTripleGaussian(double deltacp);

  double GetMulticanonicalWeightVonMises(double deltacp);

  double GetMulticanonicalWeightGenGaussian(double deltacp);

  double circularDistance(double x, double mean);

  /// bias function implementations
  double generalisedGaussian2(double x, double mean, double width);

  /// osc_cov systematic variable we wish to apply multicanonical to
  int oscCovVar;
  /// multi-canonical par number
  int multicanonicalVar;
  /// multi-canonical par number
  int multicanonicalVar_dm23;

  /// selected bias function for non-spline multicanonical weights
  BiasFunction umbrellaBiasFunction;

  /// configured bias function name for logging and compatibility
  std::string umbrellaBiasFunctionName;

  /// multi-canonical spline toggle on/off
  bool multicanonicalSpline;

protected:
  /// multi-canonical beta
  double multicanonicalBeta;
  /// osc_cov systematic variable we wish to apply multicanonical to

  /// delta_cp parameter value
  double delta_cp_value;
  /// dm23 parameter value
  double delm23_value;

  /// multi-canonical spline object
  TSpline3* dcp_spline_IO;
  TSpline3* dcp_spline_NO;

  /// umbrella mean
  double umbrellaMean;
  /// umbrella width
  double umbrellaWidth;
  /// umbrella number
  int umbrellaNumber;
  /// Toggle for setting umbrella widths based on umbrella overlap
  bool umbrellaOverlapMode;
  /// the desired overlap of evenly placed umbrellas
  double umbrellaSigmaOverlap;
  /// umbrella auto adjust step scale mode
  bool umbrellaAdjustStepScale;
  /// umbrella relative step scale
  double umbrellaStepScaleFactor;
  /// flip window toggle
  bool flipWindow;

  /// von Mises kappa parameter
  double vonMises_kappa;
  /// von Mises I0(kappa) precomputed value
  double vonMises_I0_kappa;

#ifdef DEBUG
  std::ostream* debugStream = nullptr;
  bool debugEnabled = false;
#endif
};
