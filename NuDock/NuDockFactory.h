#pragma once

/// @file NuDockFactory.h
/// @brief Factory utilities for creating and configuring NuDock communication objects.
///
/// Provides helper functions to initialise NuDock instances from MaCh3 configuration, 
/// as well as conversion utilities between MaCh3 and NuDock oscillation parameter conventions.
/// @author Hank Hua

#include <Manager/Manager.h>
#include "Parameters/ParameterHandlerGeneric.h"
_MaCh3_Safe_Include_Start_ //{
#include "nudock.hpp"
_MaCh3_Safe_Include_End_ //}

/// @brief SampleNames tag marking a parameter as belonging to the NuDock interface.
constexpr const char* kNuDockSampleTag = "NuDock";

/// @brief Initialise a NuDock communication object from a manager configuration.
///
/// @param man Pointer to the MaCh3 manager holding the NuDock YAML config block.
/// @param[in,out] nudock_ptr Unique pointer that will own the newly created NuDock object.
/// @throw MaCh3Exception if no NuDock configuration is found, or if an
///        unsupported communication type / verbosity level is specified.
void InitialiseNuDockObj(Manager *man, std::unique_ptr<NuDock> &nudock_ptr);

/// @brief Convert an oscillation parameter value from MaCh3 convention to NuDock convention.
///
/// @param param_name The MaCh3-side parameter name (e.g. "sin2th_12").
/// @param[in,out] param_value The parameter value to convert in-place.
void FormatOscParsForNuDock(const std::string &param_name, double &param_value);

/// @brief Convert an oscillation parameter value from NuDock convention to MaCh3 convention.
///
/// @param param_name The MaCh3-side parameter name (e.g. "sin2th_12").
/// @param[in,out] param_value The parameter value to convert in-place.
void FormatOscParsForMaCh3(const std::string &param_name, double &param_value);

/// @brief Collect the parameters that make up the NuDock oscillation-parameter set.
///
/// @param parHandler Parameter handler to query.
/// @return Indices of the tagged oscillation parameters.
std::vector<int> GetNuDockOscParIndices(const ParameterHandlerGeneric *parHandler);

/// @brief Throw an error if any standard oscillation parameters are not tagged for NuDock.
///
/// @param parHandler Parameter handler to check.
/// @throw MaCh3Exception if a legacy oscillation parameter is present but untagged.
void AssertNuDockOscParamsTagged(const ParameterHandlerGeneric *parHandler);

/// @brief Map from NuDock oscillation parameter names to MaCh3 names.
///
/// Keys are NuDock-style names (e.g. "Theta12"), values are MaCh3-style names (e.g. "sin2th_12").
extern const std::unordered_map<std::string, std::string> NuDockOscNameMap;

/// @brief Map from MaCh3 oscillation parameter names to NuDock names.
///
/// Keys are MaCh3-style names (e.g. "sin2th_12"), values are NuDock-style names (e.g. "Theta12").
extern const std::unordered_map<std::string, std::string> NuDockOscNameMap_r;
