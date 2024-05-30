#pragma once

#include <cmath>
#include <map>
#include <string>

//MaCh3 includes
#include "manager/YamlHelper.h"
#include "manager/Monitor.h"
#include "manager/MaCh3Exception.h"

/// Enumerator of MaCh3Mode
using MaCh3Modes_t = int;

/// @brief KS: Class containing information for a single MaCh3Mode
struct MaCh3ModeInfo {
  /// Mode name
  std::string Name;
  /// Mode fancy name
  std::string FancyName;
  /// Mapping between mode and generator integers
  std::vector<int> GeneratorMaping;
  /// @brief KS: Checks MaCh3 modes is associated with a given generator mode
  inline bool IsMode(const int GenMode) {
    bool exists = std::find(GeneratorMaping.begin(), GeneratorMaping.end(), GenMode) != GeneratorMaping.end();
    return exists;
  }
};

/// @brief KS: Class describing MaCh3 modes used in the analysis, it is being initialised from config
class MaCh3Modes {
 public:
  /// @brief KS: Initialise MaCh3 modes
  MaCh3Modes(std::string const &filename);
  /// @brief KS: Empty destructor
  virtual ~MaCh3Modes(){};

  /// @brief KS: Print info about initialised modes
  void Print();

  /// @brief KS: Get number of modes, keep in mind actual number is +1 greater due to unknown category
  inline int GetNModes(){return NModes;}
  /// @brief KS: Get mode number based on name, if mode not known you will get UNKNOWN_BAD
  MaCh3Modes_t GetMode(std::string name);
  /// @brief KS: Get normal name of mode, if mode not known you will get UNKNOWN_BAD
  std::string GetMaCh3ModeName(const int Index);
  /// @brief KS: Get fancy name of mode, if mode not known you will get UNKNOWN_BAD
  std::string GetMaCh3ModeFancyName(const int Index);
  /// @brief KS: Get MaCh3 mode from generator mode
  MaCh3Modes_t GetModeFromGenerator(const int Index);
  /// @brief Get class name
  inline std::string GetName()const {return "MaCh3Modes";};

 private:
  /// @brief KS: Make sure we don't have two modes with the same name
  inline MaCh3Modes_t EnsureModeNameRegistered(std::string const &name);

  /// @brief KS: Add new mode
  inline void DeclareNewMode(std::string const &name,
                      std::string const &fancyname,
                      std::vector<int> const &GenMap);

  /// @brief KS: Fill ModeMap
  inline void PrepareMap();

  /// KS: Handy map which helps find mode number based on string
  std::map<std::string, MaCh3Modes_t> Mode;
  /// KS: Main map storing info about used modes
  std::map<MaCh3Modes_t, MaCh3ModeInfo> fMode;
  /// KS: Handy map helping us find MaCh3 mode based on Generator mode value
  std::vector<int> ModeMap;

  /// KS: Name of loaded modes
  std::string Title;
  /// KS: Name of generator like NEUT, NuWro etc. this is to make stuff fancy
  std::string Generator;
  /// KS: Number of modes, keep in mind actual number is +1 greater due to unknown category
  int NModes;
};
