#pragma once

// C++ includes
#include <math.h>
#include <map>

// ROOT includes
#include "TList.h"

// MaCh3 includes
#include "covariance/covarianceBase.h"

/// @brief Class responsible for handling of systematic error parameters with different types defined in the config. Like spline, normalisation parameters etc.
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/02.-Implementation-of-Systematic).
class covarianceXsec : public covarianceBase {
  public:
    /// @brief Constructor
    /// @param FileNames A vector of strings representing the YAML files used for initialisation of matrix
    /// @param name Matrix name
    /// @param threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
    /// @param FirstPCAdpar First PCA parameter that will be decomposed.
    /// @param LastPCAdpar First PCA parameter that will be decomposed.
    covarianceXsec(const std::vector<std::string>& FileNames, const char *name = "xsec_cov", double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
    /// @brief Destructor
    ~covarianceXsec();

    /// @brief Print information about the whole object once it is set
    inline void Print();

    /// @brief KS: Check if matrix is correctly initialised
    void CheckCorrectInitialisation();

    // General Getter functions not split by detector
    /// @brief ETA - just return the int of the DetID, this can be removed to do a string comp at some point.
    /// @param i parameter index
    inline int GetParDetID(const int i) const { return _fDetID[i];};
    /// @brief ETA - just return a string of "spline", "norm" or "functional"
    /// @param i parameter index
    inline std::string GetParamTypeString(const int i) const { return SystType_ToString(_fParamType[i]); }
    /// @brief Returns enum describing our param type
    /// @param i parameter index
    inline SystType GetParamType(const int i) const {return _fParamType[i];}

    /// @brief Get interpolation type for a given parameter
    /// @param i spline parameter index, not confuse with global index
    inline SplineInterpolation GetParSplineInterpolation(const int i) {return SplineParams.at(i).SplineInterpolationType;}

    /// @brief EM: value at which we cap spline knot weight
    /// @param i spline parameter index, not confuse with global index
    inline double GetParSplineKnotUpperBound(const int i) {return SplineParams.at(i).SplineKnotUpBound;}
    /// @brief EM: value at which we cap spline knot weight
    /// @param i spline parameter index, not confuse with global index
    inline double GetParSplineKnotLowerBound(const int i) {return SplineParams.at(i).SplineKnotLowBound;}

    /// @brief DB Grab the number of parameters for the relevant DetID
    /// @param Type Type of syst, for example kNorm, kSpline etc
    int GetNumParamsFromDetID(const int DetID, const SystType Type);
    /// @brief DB Grab the parameter names for the relevant DetID
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<std::string> GetParsNamesFromDetID(const int DetID, const SystType Type);
    /// @brief DB Grab the parameter indices for the relevant DetID
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<int> GetParsIndexFromDetID(const int DetID, const SystType Type);

    /// @brief DB Get spline parameters depending on given DetID
    const std::vector<std::string> GetSplineParsNamesFromDetID(const int DetID);
    /// @brief DB Get spline parameters depending on given DetID
    const std::vector<std::string> GetSplineFileParsNamesFromDetID(const int DetID);

    /// ETA - what does this even do?
    const std::vector<std::string> GetFDSplineFileParsNamesFromDetID(const int DetID);
    const std::vector<std::string> GetNDSplineFileParsNamesFromDetID(const int DetID);
    /// @brief DB Grab the Spline Modes for the relevant DetID
    const std::vector< std::vector<int> > GetSplineModeVecFromDetID(const int DetID);
    /// @brief DB Grab the Spline Indices for the relevant DetID
    const std::vector<int> GetSplineParsIndexFromDetID(const int DetID){return GetParsIndexFromDetID(DetID, kSpline);}

    /// @brief DB Grab the Number of splines for the relevant DetID
    int GetNumSplineParamsFromDetID(const int DetID){return GetNumParamsFromDetID(DetID, kSpline);}

    /// @brief DB Get norm/func parameters depending on given DetID
    const std::vector<XsecNorms4> GetNormParsFromDetID(const int DetID);

    /// @brief DB Grab the number of Normalisation parameters for the relevant DetID
    int GetNumFuncParamsFromDetID(const int DetID){return GetNumParamsFromDetID(DetID, kFunc);}
    /// @brief DB Grab the Functional parameter names for the relevant DetID
    const std::vector<std::string> GetFuncParsNamesFromDetID(const int DetID){return GetParsNamesFromDetID(DetID, kFunc);}
    /// @brief DB Grab the Functional parameter indices for the relevant DetID
    const std::vector<int> GetFuncParsIndexFromDetID(const int DetID){return GetParsIndexFromDetID(DetID, kFunc);}

    /// @brief KS: For most covariances nominal and fparInit (prior) are the same, however for Xsec those can be different
    /// For example Sigma Var are done around nominal in ND280, no idea why though...
    std::vector<double> getNominalArray() override
    {
      std::vector<double> nominal(_fNumPar);
      for (int i = 0; i < _fNumPar; i++) {
        nominal[i] = _fPreFitValue.at(i);
      }
      return nominal;
    }
    /// @brief Get nominal for a given param
    /// @param i parameter index
    inline double getNominal(const int i) override { return _fPreFitValue.at(i); };
    /// @brief Is parameter a flux param or not. This might become deprecated in future
    /// @param i parameter index
    /// @warning Will become deprecated
    inline bool IsParFlux(const int i){ return isFlux[i]; }
    /// @brief KS Function to set to nominal flux parameters
    /// @warning Will become deprecated
    void setXsecOnlyParameters();
    /// @brief KS Function to set to nominal flux  parameters
    /// @warning Will become deprecated
    void setFluxOnlyParameters();
    
  protected:
    /// @brief Initialise CovarianceXsec
    void initParams();
    /// @brief ETA - trying out the yaml parsing
    inline void InitXsecFromConfig();
    /// @brief Get Norm params
    /// @param param Yaml node describing param
    /// @param Index Global parameter index
    inline XsecNorms4 GetXsecNorm(const YAML::Node& param, const int Index);
    /// @brief Get Spline params
    /// @param param Yaml node describing param
    inline XsecSplines1 GetXsecSpline(const YAML::Node& param);

    /// Is parameter flux or not, This might become deprecated in future
    /// @warning Will become deprecated
    std::vector<bool> isFlux;

    /// Tells to which samples object param should be applied
    std::vector<int> _fDetID;
    /// Type of parameter like norm, spline etc.
    std::vector<SystType> _fParamType;

    //Variables related to spline systematics
    std::vector<std::string> _fNDSplineNames;
    std::vector<std::string> _fFDSplineNames;
    std::vector<std::vector<int>> _fFDSplineModes;

    /// Vector containing info for normalisation systematics
    std::vector<XsecSplines1> SplineParams;

    /// Vector containing info for normalisation systematics
    std::vector<XsecNorms4> NormParams;
};
