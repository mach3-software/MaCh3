#pragma once

// C++ includes
#include <cstdlib>
#include <iomanip>

// MaCh3  includes
#include "Splines/SplineStructs.h"
#include "Splines/SplineCommon.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT include
#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"
#include "TString.h"
#include "TIterator.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
_MaCh3_Safe_Include_End_ //}

/// @brief Base class for calculating weight from spline
class SplineBase {
  public:
    /// @brief Constructor
    SplineBase();
    /// @brief Destructor
    virtual ~SplineBase();

    /// @brief  CW: This Eval should be used when using two separate x,{y,a,b,c,d} arrays to store the weights; probably the best one here! Same thing but pass parameter spline segments instead of variations
    virtual void Evaluate() = 0;

    /// @brief Get class name
    virtual inline std::string GetName()const {return "SplineBase";};

    /// @brief Get number of spline parameters
    short int GetNParams()const {return nParams;};

  protected:
    /// @brief CW:Code used in step by step reweighting, Find Spline Segment for each param
    void FindSplineSegment();
    /// @brief CPU based code which eval weight for each spline
    virtual void CalcSplineWeights() = 0;
    /// @brief Calc total event weight
    virtual void ModifyWeights() = 0;

    /// @brief CW: Gets the polynomial coefficients for TF1
    /// @param spl pointer to TF1_red that will be checked
    /// @param nPoints number of knots
    /// @param coeffs Array holding coefficients for each knot
    void getTF1Coeff(TF1_red* &spl, int &nPoints, float *&coeffs);

    /// Array of FastSplineInfo structs: keeps information on each xsec spline for fast evaluation
    /// Method identical to TSpline3::Eval(double) but faster because less operations
    std::vector<FastSplineInfo> SplineInfoArray;
    /// Store currently found segment they are not in FastSplineInfo as in case of GPU we need to copy paste it to GPU
    /// @warning this is being used sometimes by GPU, therefore must be raw pointer!
    short int *SplineSegments;
    /// Store parameter values they are not in FastSplineInfo as in case of GPU we need to copy paste it to GPU
    float *ParamValues;
    /// Number of parameters that have splines
    short int nParams;
};
