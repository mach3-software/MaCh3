#pragma once

// MaCh3 includes
#include "Parameters/ParameterHandlerBase.h"
#include "Samples/SampleStructs.h"

/// @brief Class responsible for handling of systematic error parameters with different types defined in the config. Like spline, normalisation parameters etc.
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/02.-Implementation-of-Systematic).
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski
class ParameterHandlerGeneric : public ParameterHandlerBase {
  public:
    /// @brief Constructor
    /// @param FileNames A vector of strings representing the YAML files used for initialisation of matrix
    /// @param name Matrix name
    /// @param threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
    /// @param FirstPCAdpar First PCA parameter that will be decomposed.
    /// @param LastPCAdpar First PCA parameter that will be decomposed.
    ParameterHandlerGeneric(const std::vector<std::string>& FileNames, std::string name = "xsec_cov", double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
    /// @brief Destructor
    ~ParameterHandlerGeneric();

    // General Getter functions not split by detector
    /// @brief ETA - just return the int of the SampleName, this can be removed to do a string comp at some point.
    /// @param i parameter index
    inline std::vector<std::string> GetParSampleID(const int i) const { return _fSampleNames[i];};
    /// @brief ETA - just return a string of "spline", "norm" or "functional"
    /// @param i parameter index
    inline std::string GetParamTypeString(const int i) const { return SystType_ToString(_fParamType[i]); }
    /// @brief Returns enum describing our param type
    /// @param i parameter index
    inline SystType GetParamType(const int i) const {return _fParamType[i];}

    /// @brief Get interpolation type for a given parameter
    /// @param i spline parameter index, not confuse with global index
    inline SplineInterpolation GetParSplineInterpolation(const int i) const {return SplineParams.at(i)._SplineInterpolationType;}
    /// @brief Get the interpolation types for splines affecting a particular SampleName
    const std::vector<SplineInterpolation> GetSplineInterpolationFromSampleName(const std::string& SampleName);
    /// @brief Get the name of the spline associated with the spline at index i
    /// @param i spline parameter index, not to be confused with global index
    std::string GetParSplineName(const int i) const {return _fSplineNames[i];}

    /// @brief DB Get spline parameters depending on given SampleName
    const std::vector<int> GetGlobalSystIndexFromSampleName(const std::string& SampleName, const SystType Type);
    /// @brief EM: value at which we cap spline knot weight
    /// @param i spline parameter index, not confuse with global index
    inline double GetParSplineKnotUpperBound(const int i) const {return SplineParams.at(i)._SplineKnotUpBound;}
    /// @brief EM: value at which we cap spline knot weight
    /// @param i spline parameter index, not confuse with global index
    inline double GetParSplineKnotLowerBound(const int i) const {return SplineParams.at(i)._SplineKnotLowBound;}

    /// @brief DB Grab the number of parameters for the relevant SampleName
    /// @param SampleName property of SampleHandler class based on which we select whether to apply uncertainties or not
    /// @param Type Type of syst, for example kNorm, kSpline etc
    int GetNumParamsFromSampleName(const std::string& SampleName, const SystType Type);
    /// @brief DB Grab the parameter names for the relevant SampleName
    /// @param SampleName property of SampleHandler class based on which we select whether to apply uncertainties or not
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<std::string> GetParsNamesFromSampleName(const std::string& SampleName, const SystType Type);
    /// @brief DB Grab the parameter indices for the relevant SampleName
    /// @param SampleName property of SampleHandler class based on which we select whether to apply uncertainties or not
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<int> GetParsIndexFromSampleName(const std::string& SampleName, const SystType Type);

    /// @brief DB Get spline parameters depending on given SampleName
    const std::vector<std::string> GetSplineParsNamesFromSampleName(const std::string& SampleName);
    /// @brief DB Get spline parameters depending on given SampleName
    const std::vector<std::string> GetSplineFileParsNamesFromSampleName(const std::string& SampleName);

    /// @brief DB Grab the Spline Modes for the relevant SampleName
    const std::vector< std::vector<int> > GetSplineModeVecFromSampleName(const std::string& SampleName);
    /// @brief Grab the index of the syst relative to global numbering.
    /// @param SampleName property of SampleHandler class based on which we select whether to apply uncertainties or not
    /// @param Type Type of syst, for example kNorm, kSpline etc
    const std::vector<int> GetSystIndexFromSampleName(const std::string& SampleName, const SystType Type) const;
    /// @brief DB Get norm/func parameters depending on given SampleName
    const std::vector<NormParameter> GetNormParsFromSampleName(const std::string& SampleName) const;
    /// @brief HH Get functional parameters for the relevant SampleName
    const std::vector<FunctionalParameter> GetFunctionalParametersFromSampleName(const std::string& SampleName) const;
    /// @brief KS: Grab the Spline parameters for the relevant SampleName
    const std::vector<SplineParameter> GetSplineParsFromSampleName(const std::string& SampleName) const;

    /// @brief Checks if parameter belongs to a given group
    /// @param i parameter index
    /// @param Group name of group, like Xsec or Flux
    /// @return bool telling whether param is part of group
    bool IsParFromGroup(const int i, const std::string& Group) const;

    /// @brief KS: Check how many parameters are associated with given group
    int GetNumParFromGroup(const std::string& Group) const;
    /// @brief KS: Get names of all unique parameter groups
    std::vector<std::string> GetUniqueParameterGroups();

    /// @brief KS Function to set to prior parameters of a given group or values from vector
    /// @param Group name of group, like Xsec or Flux
    /// @param Pars Values which will overwrite proposed step
    /// @note this mimic functionality of @ParameterHandlerBase::SetParameters
    void SetGroupOnlyParameters(const std::string& Group, const std::vector<double>& Pars = {});
    /// @brief KS Function to set to prior parameters of a given groups or values from vector
    /// @param Groups vector of group names, like Xsec or Flux
    void SetGroupOnlyParameters(const std::vector<std::string>& Groups);
    
    /// @brief TN Method to set parameters within a group to be fixed to their prior values
    /// @param Group name of the parameter group (Xsec, Flux, Osc, etc.)
    void SetFixGroupOnlyParameters(const std::string& Group);
    /// @brief TN Method to set parameters of certain groups to be fixed to their prior values
    /// @param Groups vector of group names (e.g. {"Xsec", "Flux"})
    void SetFixGroupOnlyParameters(const std::vector<std::string>& Groups);

    /// @brief TN Method to set parameters within a group to be treated as free
    /// @param Group name of the parameter group (Xsec, Flux, Osc, etc.)
    void SetFreeGroupOnlyParameters(const std::string& Group);
    /// @brief TN Method to set parameters of certain groups to be treated as free
    /// @param Groups vector of group names (e.g. {"Xsec", "Flux"})
    void SetFreeGroupOnlyParameters(const std::vector<std::string>& Groups);

    /// @brief TN Method to toggle fix/free parameters within a group
    /// @param Group name of the parameter group (Xsec, Flux, Osc, etc.)
    void ToggleFixGroupOnlyParameters(const std::string& Group);   
    /// @brief TN Method to toggle fix/free parameters within given groups
    /// @param Group vector of group names (e.g. {"Xsec", "Flux"})
    void ToggleFixGroupOnlyParameters(const std::vector<std::string>& Groups);

    /// @brief Dump Matrix to ROOT file, useful when we need to pass matrix info to another fitting group
    /// @param Name Name of TFile to which we save stuff
    /// @warning This is mostly used for backward compatibility
    void DumpMatrixToFile(const std::string& Name);

    /// @brief Get pointers to Osc params from Sample name
    std::vector<const double*> GetOscParsFromSampleName(const std::string& SampleName);

  protected:
    /// @brief Print information about the whole object once it is set
    void Print();
    /// @brief Prints general information about the ParameterHandler object.
    void PrintGlobablInfo();
    /// @brief Prints normalization parameters.
    void PrintNormParams();
    /// @brief Prints spline parameters.
    void PrintSplineParams();
    /// @brief Prints functional parameters.
    void PrintFunctionalParams();
    /// @brief Prints oscillation parameters.
    void PrintOscillationParams();
    /// @brief Prints groups of parameters.
    void PrintParameterGroups();

    /// @brief KS: Check if matrix is correctly initialised
    void CheckCorrectInitialisation();

    /// @brief Iterates over parameters and applies a filter and action function.
    ///
    /// This template function provides a way to iterate over parameters associated
    /// with a specific Sample ID (SampleName). It applies a filter function to determine
    /// which parameters to process and an action function to define what to do
    /// with the selected parameters.
    ///
    /// @tparam FilterFunc The type of the filter function used to determine
    /// which parameters to include.
    /// @tparam ActionFunc The type of the action function applied to each selected
    /// parameter.
    /// @param SampleName The Sample ID used to filter parameters.
    template <typename FilterFunc, typename ActionFunc>
    void IterateOverParams(const std::string& SampleName, FilterFunc filter, ActionFunc action);

    /// @brief Initializes the systematic parameters from the configuration file.
    /// This function loads parameters like normalizations and splines from the provided YAML file.
    /// @note This is used internally during the object's initialization process.
    void InitParams();

    /// @brief Parses the YAML configuration to set up cross-section parameters.
    /// The YAML file defines the types of systematic errors, interpolation types, and bounds for splines.
    inline void InitParametersTypeFromConfig();

    /// @brief Get Norm params
    /// @param param Yaml node describing param
    /// @param Index Global parameter index
    inline NormParameter GetNormParameter(const YAML::Node& param, const int Index);

    /// @brief Get Osc params
    /// @param param Yaml node describing param
    /// @param Index Global parameter index
    inline OscillationParameter GetOscillationParameters(const YAML::Node& param, const int Index);

    /// @brief Get Func params
    /// @param param Yaml node describing param
    /// @param Index Global parameter index
    inline FunctionalParameter GetFunctionalParameters(const YAML::Node& param, const int Index);
    /// @brief Get Spline params
    /// @param param Yaml node describing param
    /// @param Index Global parameter index
    inline SplineParameter GetSplineParameter(const YAML::Node& param, const int Index);
    /// @brief Fill base parameters
    /// @param param Yaml node describing param
    /// @param Index Global parameter index
    /// @param Parameter Object storing info
    inline void GetBaseParameter(const YAML::Node& param, const int Index, TypeParameterBase& Parameter);

    /// @brief Retrieve parameters that apply to a given sample name.
    /// @tparam ParamT Type of parameter (e.g., FunctionalParameter, NormParameter).
    /// @param indexMap Map from local to global parameter indices.
    /// @param params Vector of all parameters of the given type.
    /// @param SampleName The name of the sample to filter applicable parameters for.
    /// @return Vector of parameters of type ParamT that apply to the specified sample.
    template<typename ParamT>
    std::vector<ParamT> GetTypeParamsFromSampleName(const std::map<int, int>& indexMap, const std::vector<ParamT>& params, const std::string& SampleName) const;

    /// Type of parameter like norm, spline etc.
    std::vector<SystType> _fParamType;

    /// Name of spline in TTree (TBranch),
    std::vector<std::string> _fSplineNames;

    /// KS: Allow to group parameters for example to affect only cross-section or only flux etc.
    std::vector<std::string> _ParameterGroup;

    /// Map between number of given parameter type with global parameter numbering. For example 2nd norm param may be 10-th global param
    std::vector<std::map<int, int>> _fSystToGlobalSystIndexMap;

    /// Vector containing info for normalisation systematics
    std::vector<SplineParameter> SplineParams;

    /// Vector containing info for normalisation systematics
    std::vector<NormParameter> NormParams;

    /// Vector containing info for functional systematics
    std::vector<FunctionalParameter> FuncParams;

    /// Vector containing info for functional systematics
    std::vector<OscillationParameter> OscParams;
};
