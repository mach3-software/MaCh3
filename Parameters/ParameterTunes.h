#pragma once

// MaCh3 Includes
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerUtils.h"

/// @brief Class storing tune values, for example allowing flexibly to retrieve parameter values for for PostND tunes
/// @author Kamil Skwarczynski
class ParameterTunes{
 public:
  /// @brief Constructor
  ParameterTunes(){};
  ParameterTunes(const YAML::Node& Settings);
  /// @brief
  virtual ~ParameterTunes();

  /// @brief Return vector of tune vales for each parameter
  std::vector<double> GetTune(const int TuneNumber) const;
  /// @brief Return vector of tune vales for each parameter
  std::vector<double> GetTune(const std::string& TuneName) const;
  /// @brief Simply print all tunes and associated values
  void PrintTunes() const;
 private:
  /// Name of each Tune
  std::vector<std::string> TuneNames;
  /// Values for each Tune and Parameter
  std::vector<std::vector<double>> TuneValues;
  /// Map between tune name and value
  std::unordered_map<std::string, int> TuneMap;
};
