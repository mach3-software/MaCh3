#pragma once

// C++ includes
#include <math.h>
#include <map>

// ROOT includes
#include "TList.h"

// MaCh3 includes
#include "covariance/covarianceBase.h"

/// @brief Class responsible for handling of systematic error parameters with different types defined in the config. Like spline , normalisation parameters etc.
class covarianceXsec : public covarianceBase {

  public:
    /// @brief Constructor
    covarianceXsec(std::vector<std::string> FileNames, double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
    /// @brief Destructor
    ~covarianceXsec();

    /// @brief ETA - trying out the yaml parsing
    inline void InitXsecFromConfig();
    /// @brief Initialise Norm params
    inline void SetupNormPars();
    /// @brief Print information about the whole object once it is set
    inline void Print();

    // General Getter functions not split by detector
    /// @brief ETA - just return the int of the DetID, this can be removed to do a string comp at some point.
    int GetParDetID(const int i) const { return _fDetID[i];};
    /// @brief ETA - just return a string of "spline", "norm" or "functional"
    const char*  GetParamType(const int i) const {return _fParamType[i].c_str();}

    /// @brief Get interpolation type vector
    const std::vector<SplineInterpolation>& GetSplineInterpolation() const{return _fSplineInterpolationType;}
    /// @brief Get interpolation type for a given parameter
    SplineInterpolation GetParSplineInterpolation(int i) {return _fSplineInterpolationType.at(i);}

    /// @brief EM: value at which we cap spline knot weight
    double GetParSplineKnotUpperBound(const int i) {return fSplineKnotUpBound[i];}
    /// @brief EM: value at which we cap spline knot weight
    double GetParSplineKnotLowerBound(const int i) {return fSplineKnotLowBound[i];}

    /// @brief DB Get spline parameters depending on given DetID
    const std::vector<std::string> GetSplineParsNamesFromDetID(int DetID);
    /// @brief DB Get spline parameters depending on given DetID
    const std::vector<std::string> GetSplineFileParsNamesFromDetID(int DetID);

    /// ETA - what does this even do?
    const std::vector<std::string> GetFDSplineFileParsNamesFromDetID(int DetID);
    const std::vector<std::string> GetNDSplineFileParsNamesFromDetID(int DetID);
    const std::vector< std::vector<int> > GetSplineModeVecFromDetID(int DetID);
    const std::vector<int> GetSplineParsIndexFromDetID(int DetID);
    int GetNumSplineParamsFromDetID(int DetID);

    /// @brief DB Get norm/func parameters depending on given DetID
    const std::vector<XsecNorms4> GetNormParsFromDetID(int DetID);

    int GetNumFuncParamsFromDetID(int DetID);
    const std::vector<std::string> GetFuncParsNamesFromDetID(int DetID);
    const std::vector<int> GetFuncParsIndexFromDetID(int DetID);

    /// KS: For most covariances nominal and fparInit (prior) are the same, however for Xsec those can be different
    /// For example Sigma Var are done around nominal in ND280, no idea why though...
    std::vector<double> getNominalArray() override
    {
      std::vector<double> nominal;
      for (int i = 0; i < size; i++)
      {
        nominal.push_back(_fPreFitValue.at(i));
      }
      return nominal;
    }
    /// @brief Get nominal for a given param
    inline double getNominal(const int i) override { return _fPreFitValue.at(i); };
   
    /// @brief Is parameter a flux param or not. This might become deprecated in future
    bool IsParFlux(const int i){ return isFlux[i]; }

    /// @brief KS Function to set to nominal flux parameters
    void setXsecOnlyParameters();
    /// @brief KS Function to set to nominal flux  parameters
    void setFluxOnlyParameters();
    
  protected:
    /// @brief Initialise CovarianceXsec
    void initParams(const double fScale);
    /// Is parameter flux or not, This might become deprecated in future
    std::vector<bool> isFlux;

    /// Tells to which samples object param should be applied
    std::vector<int> _fDetID;
    //std::vector<std::string> _fDetString;
    /// Type of parameter like norm, spline etc.
    std::vector<std::string> _fParamType;

    //Some "usual" variables. Don't think we really need the ND/FD split
    std::vector<std::vector<int>> _fNormModes;
    std::vector<std::vector<int>> _fTargetNuclei;
    std::vector<std::vector<int>> _fNeutrinoFlavour;
    std::vector<std::vector<int>> _fNeutrinoFlavourUnosc;

    //Variables related to spline systematics
    std::vector<std::string> _fNDSplineNames;
    std::vector<std::string> _fFDSplineNames;
    std::vector<std::vector<int>> _fFDSplineModes;
    /// Spline interpolation vector
    std::vector<SplineInterpolation> _fSplineInterpolationType;

    /// EM: Cap spline knot lower value
    std::vector<double> fSplineKnotLowBound;
    /// EM: Cap spline knot higher value
    std::vector<double> fSplineKnotUpBound;

    /// Information to be able to apply generic cuts
    std::vector<std::vector<std::string>> _fKinematicPars;
    /// Information to be able to apply generic cuts
    std::vector<std::vector<std::vector<double>>> _fKinematicBounds;

    /// Vector containing info for normalisation systematics
    std::vector<XsecNorms4> NormParams;
};
