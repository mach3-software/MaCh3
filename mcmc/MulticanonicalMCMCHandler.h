#pragma once

#include "manager/manager.h"
#include "covariance/covarianceBase.h"
#include <vector>
#include "TSpline.h"

class MulticanonicalMCMCHandler  {
  public:
    /// @brief Constructor
    MulticanonicalMCMCHandler();
    /// @brief Destructor
    virtual ~MulticanonicalMCMCHandler();

    void InitializeMulticanonicalHandlerConfig(manager* fitMan, std::vector<covarianceBase*>& systematics);

    void InitializeMulticanonicalParams(std::vector<covarianceBase*>& systematics);

    double GetMulticanonicalWeightSpline(double deltacp, double delm23_value);

    double GetMulticanonicalWeightSeparate(double deltacp);

    double GetMulticanonicalWeightGaussian(double deltacp);

    double GetMulticanonicalWeightVonMises(double deltacp);

    /// osc_cov systematic variable we wish to apply multicanonical to
    int oscCovVar;
    /// multi-canonical par number
    int multicanonicalVar;
    /// multi-canonical par number
    int multicanonicalVar_dm23;


    /// multi-canonical separate toggle on/off
    bool multicanonicalSeparate;

    /// multi-canonical spline toggle on/off
    bool multicanonicalSpline;
  protected:
    /// multi-canonical beta
    double multicanonicalBeta;
    /// multi-canonical sigma
    double multicanonicalSigma;
    /// osc_cov systematic variable we wish to apply multicanonical to

    /// delta_cp parameter value
    double delta_cp_value;
    /// dm23 parameter value
    double delm23_value;

     /// multi-canonical spline object
    TSpline3 *dcp_spline_IO;
    TSpline3 *dcp_spline_NO;
    
    /// multi-canonical separate mean
    double multicanonicalSeparateMean;
    /// multi-canonical separate sigma
    double multicanonicalSeparateSigma;
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

    /// von Mises mode toggle
    bool vonMises_mode;
    /// von Mises kappa parameter
    double vonMises_kappa;
    /// von Mises I0(kappa) precomputed value
    double vonMises_I0_kappa;

};
