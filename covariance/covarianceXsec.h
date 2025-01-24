#pragma once

// MaCh3 includes
#include "covariance/covarianceBase.h"
#include "samplePDF/Structs.h"

/// @brief Class responsible for handling of systematic error parameters with different types defined in the config. Like spline, normalisation parameters etc.
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/02.-Implementation-of-Systematic).
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski
class covarianceXsec : public covarianceBase {
  public:
    /// @brief Constructor
    /// @param FileNames A vector of strings representing the YAML files used for initialisation of matrix
    /// @param name Matrix name
    /// @param threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
    /// @param FirstPCAdpar First PCA parameter that will be decomposed.
    /// @param LastPCAdpar First PCA parameter that will be decomposed.
    covarianceXsec(const std::vector<std::string>& FileNames, std::string name = "xsec_cov", double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
    /// @brief Destructor
    ~covarianceXsec();

    // General Getter functions not split by detector
    /// @brief ETA - just return the int of the DetID, this can be removed to do a string comp at some point.
    /// @param i parameter index
    inline std::vector<std::string> GetParDetID(const int i) const { return _fDetID[i];};
    /// @brief ETA - just return a string of "spline", "norm" or "functional"
    /// @param i parameter index
    inline std::string GetParamTypeString(const int i) const { return SystType_ToString(_fParamType[i]); }
    /// @brief Returns enum describing our param type
    /// @param i parameter index
    inline SystType GetParamType(const int i) const {return _fParamType[i];}

    /// @brief Get interpolation type for a given parameter
    /// @param i spline parameter index, not confuse with global index
    inline SplineInterpolation GetParSplineInterpolation(const int i) {return SplineParams.at(i)._SplineInterpolationType;}
    /// @brief Get the interpolation types for splines affecting a particular DetID
    const std::vector<SplineInterpolation> GetSplineInterpolationFromDetID(const std::string& DetID);
    /// @brief Get the name of the spline associated with the spline at index i
    /// @param i spline parameter index, not to be confused with global index
    std::string GetParSplineName(const int i) {return _fSplineNames[i];}

    /// @brief DB Get spline parameters depending on given DetID
    const std::vector<int> GetGlobalSystIndexFromDetID(const std::string& DetID, const SystType Type);
    /// @brief EM: value at which we cap spline knot weight
    /// @param i spline parameter index, not confuse with global index
    inline double GetParSplineKnotUpperBound(const int i) {return SplineParams.at(i)._SplineKnotUpBound;}
    /// @brief EM: value at which we cap spline knot weight
    /// @param i spline parameter index, not confuse with global index
    inline double GetParSplineKnotLowerBound(const int i) {return SplineParams.at(i)._SplineKnotLowBound;}

    /// @brief DB Grab the number of parameters for the relevant DetID
    /// @param Type Type of syst, for example kNorm, kSpline etc
    int GetNumParamsFromDetID(const std::string& DetID, const SystType Type);
    /// @brief DB Grab the parameter names for the relevant DetID
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<std::string> GetParsNamesFromDetID(const std::string& DetID, const SystType Type);
    /// @brief DB Grab the parameter indices for the relevant DetID
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<int> GetParsIndexFromDetID(const std::string& DetID, const SystType Type);

    /// @brief DB Get spline parameters depending on given DetID
    const std::vector<std::string> GetSplineParsNamesFromDetID(const std::string& DetID);
    /// @brief DB Get spline parameters depending on given DetID
    const std::vector<std::string> GetSplineFileParsNamesFromDetID(const std::string& DetID);

    /// @brief DB Grab the Spline Modes for the relevant DetID
    const std::vector< std::vector<int> > GetSplineModeVecFromDetID(const std::string& DetID);
    /// @brief Grab the index of the syst relative to global numbering.
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<int> GetSystIndexFromDetID(const std::string& DetID, const SystType Type);
    /// @brief DB Get norm/func parameters depending on given DetID
    const std::vector<XsecNorms4> GetNormParsFromDetID(const std::string& DetID);

    /// @brief KS: For most covariances prior and fparInit (prior) are the same, however for Xsec those can be different
    std::vector<double> getNominalArray() override
    {
      std::vector<double> prior(_fNumPar);
      for (int i = 0; i < _fNumPar; i++) {
        prior[i] = _fPreFitValue.at(i);
      }
      return prior;
    }
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
    void IterateOverParams(const std::string& DetID, FilterFunc filter, ActionFunc action);

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
