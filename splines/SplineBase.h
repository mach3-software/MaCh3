#pragma once

// C++ includes
#include <cstdlib>
#include <iomanip>

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

#ifdef MULTITHREAD
#include "omp.h"
#endif

// MaCh3  includes
#include "manager/MaCh3Logger.h"
#include "splines/SplineStructs.h"
#include "splines/SplineCommon.h"

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

  protected:

    /// @brief CW:Code used in step by step reweighting, Find Spline Segment for each param
    virtual void FindSplineSegment() = 0;
    /// @brief CPU based code which eval weight for each spline
    virtual void CalcSplineWeights() = 0;
    /// @brief Calc total event weight
    virtual void ModifyWeights() = 0;
};
