#pragma once

#include <cmath>
#include <map>
#include <string>

//MaCh3 includes
#include "manager/YamlHelper.h"
#include "manager/Monitor.h"
#include "manager/MaCh3Exception.h"

using MaCh3Modes_t = int;

/// @brief Class containing information for a single MaCh3Mode
struct ModeInfo {
  std::string Name;

  std::string FancyName;

  std::vector<int> GeneratorMaping;

  inline bool IsMode(const int GenMode) {
    bool exists = std::find(GeneratorMaping.begin(), GeneratorMaping.end(), GenMode) != GeneratorMaping.end();
    return exists;
  }
};

/// @brief Class describing MaCh3 modes used in the analysis, it is being initalised from config
class MaCh3Modes {
 public:

  MaCh3Modes(std::string const &filename);
  virtual ~MaCh3Modes(){};

  inline int GetNModes(){return NModes;}
  MaCh3Modes_t GetMode(std::string name);
  std::string GetMaCh3ModeName(const int Index);
  std::string GetMaCh3ModeFancyName(const int Index);
  MaCh3Modes_t GetModeFromGenerator(const int Index);

  void Print();

 private:
  inline MaCh3Modes_t EnsureModeNameRegistered(std::string const &name);

  inline void DeclareNewMode(std::string const &name,
                      std::string const &fancyname,
                      std::vector<int> const &GenMap);

  inline void PrepareMap();

  std::map<std::string, MaCh3Modes_t> Mode;
  std::map<MaCh3Modes_t, ModeInfo> fMode;
  std::vector<int> ModeMap;

  YAML::Node config;

  std::string Title;
  std::string Generator;
  int NModes;
};
