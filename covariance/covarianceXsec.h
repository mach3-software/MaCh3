#pragma once

// MaCh3 includes
#include "covariance/covarianceBase.h"
#include "samplePDF/Structs.h"

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
    inline SplineInterpolation GetParSplineInterpolation(const int i) {return SplineParams.at(i)._SplineInterpolationType;}

    //DB Get spline parameters depending on given DetID
    const std::vector<int> GetGlobalSystIndexFromDetID(const int DetID, const SystType Type);
    /// @brief EM: value at which we cap spline knot weight
    /// @param i spline parameter index, not confuse with global index
    inline double GetParSplineKnotUpperBound(const int i) {return SplineParams.at(i)._SplineKnotUpBound;}
    /// @brief EM: value at which we cap spline knot weight
    /// @param i spline parameter index, not confuse with global index
    inline double GetParSplineKnotLowerBound(const int i) {return SplineParams.at(i)._SplineKnotLowBound;}

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

    /// @brief DB Grab the Spline Modes for the relevant DetID
    const std::vector< std::vector<int> > GetSplineModeVecFromDetID(const int DetID);
    /// @brief DB Grab the Spline Indices for the relevant DetID
    const std::vector<int> GetSplineParsIndexFromDetID(const int DetID){return GetParsIndexFromDetID(DetID, SystType::kSpline);}
    /// @brief ETA Grab the index of the spline relative to the _fSplineNames vector.
    const std::vector<int> GetSplineSystIndexFromDetID(const int DetID){return GetSystIndexFromDetID(DetID, SystType::kSpline);};
    /// @brief Grab the index of the syst relative to global numbering.
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<int> GetSystIndexFromDetID(const int DetID, const SystType Type);

    /// @brief DB Grab the Number of splines for the relevant DetID
    int GetNumSplineParamsFromDetID(const int DetID){return GetNumParamsFromDetID(DetID, SystType::kSpline);}

    /// @brief DB Get norm/func parameters depending on given DetID
    const std::vector<XsecNorms4> GetNormParsFromDetID(const int DetID);

    /// @brief DB Grab the number of Normalisation parameters for the relevant DetID
    int GetNumFuncParamsFromDetID(const int DetID){return GetNumParamsFromDetID(DetID, SystType::kFunc);}
    /// @brief DB Grab the Functional parameter names for the relevant DetID
    const std::vector<std::string> GetFuncParsNamesFromDetID(const int DetID){return GetParsNamesFromDetID(DetID, SystType::kFunc);}
    /// @brief DB Grab the Functional parameter indices for the relevant DetID
    const std::vector<int> GetFuncParsIndexFromDetID(const int DetID){return GetParsIndexFromDetID(DetID, SystType::kFunc);}

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

    /// @brief Checks if parameter belongs to a given group
    /// @param i parameter index
    /// @param Group name of group, like Xsec or Flux
    /// @return bool telling whether param is part of group
    bool IsParFromGroup(const int i, const std::string& Group);

    /// @brief KS Function to set to prior parameters of a given group
    /// @param Group name of group, like Xsec or Flux
    void SetGroupOnlyParameters(const std::string& Group);

    /// @brief Dump Matrix to ROOT file, useful when we need to pass matrix info to another fitting group
    /// @param Name Name of TFile to which we save stuff
    /// @warning This is mostly used for backward compatibility
    void DumpMatrixToFile(const std::string& Name);
  protected:
    /// @brief Print information about the whole object once it is set
    void Print();
    /// @brief Prints general information about the covarianceXsec object.
    void PrintGlobablInfo();
    /// @brief Prints normalization parameters.
    void PrintNormParams();
    /// @brief Prints spline parameters.
    void PrintSplineParams();
    /// @brief Prints functional parameters.
    void PrintFunctionalParams();
    /// @brief Prints groups of parameters.
    void PrintParameterGroups();

    /// @brief KS: Check if matrix is correctly initialised
    void CheckCorrectInitialisation();

    /// @brief Iterates over parameters and applies a filter and action function.
    ///
    /// This template function provides a way to iterate over parameters associated
    /// with a specific Detector ID (DetID). It applies a filter function to determine
    /// which parameters to process and an action function to define what to do
    /// with the selected parameters.
    ///
    /// @tparam FilterFunc The type of the filter function used to determine
    /// which parameters to include.
    /// @tparam ActionFunc The type of the action function applied to each selected
    /// parameter.
    /// @param DetID The Detector ID used to filter parameters.
    template <typename FilterFunc, typename ActionFunc>
    void IterateOverParams(const int DetID, FilterFunc filter, ActionFunc action);

    /// @brief Initializes the systematic parameters from the configuration file.
    /// This function loads parameters like normalizations and splines from the provided YAML file.
    /// @note This is used internally during the object's initialization process.
    void initParams();
    /// @brief Parses the YAML configuration to set up cross-section parameters.
    /// The YAML file defines the types of systematic errors, interpolation types, and bounds for splines.
    inline void InitXsecFromConfig();
    /// @brief Get Norm params
    /// @param param Yaml node describing param
    /// @param Index Global parameter index
    inline XsecNorms4 GetXsecNorm(const YAML::Node& param, const int Index);
    /// @brief Get Spline params
    /// @param param Yaml node describing param
    inline XsecSplines1 GetXsecSpline(const YAML::Node& param);

    /// Tells to which samples object param should be applied
    std::vector<int> _fDetID;
    /// Type of parameter like norm, spline etc.
    std::vector<SystType> _fParamType;

    /// Name of spline in TTree (TBranch),
    std::vector<std::string> _fSplineNames;

    /// KS: Allow to group parameters for example to affect only cross-section or only flux etc.
    std::vector<std::string> _ParameterGroup;

    /// Map between number of given parameter type with global parameter numbering. For example 2nd norm param may be 10-th global param
    std::vector<std::map<int, int>> _fSystToGlobalSystIndexMap;

    /// Vector containing info for normalisation systematics
    std::vector<XsecSplines1> SplineParams;

    /// Vector containing info for normalisation systematics
    std::vector<XsecNorms4> NormParams;
};
