#pragma once

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

// ROOT include
#include "TTree.h"
#include "TBranch.h"
#include "TMacro.h"

#include "manager/MaCh3Logger.h"

#include "samplePDF/Structs.h"
#include "manager/YamlHelper.h"
#include "manager/Monitor.h"

//KS: Joy of forward declaration https://gieseanw.wordpress.com/2018/02/25/the-joys-of-forward-declarations-results-from-the-real-world/
class TFile;

class manager {

public:
  manager(std::string const &);
  virtual ~manager();

  /// @brief Get likelihood type defined in the config
  inline int GetMCStatLLH(){return mc_stat_llh;}

  /// @brief Return name of config
  inline std::string GetFileName(){return FileName;}

  /// @brief Return config
  inline YAML::Node const &raw(){return config;}

  /// @brief Add manager useful information's to TFile, in most cases to Fitter
  void SaveSettings(TFile* const OutputFile);

  /// @brief Print currently used config
  void Print();

private:
  YAML::Node config;
  std::string FileName;
  int mc_stat_llh;

};
