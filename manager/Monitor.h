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
#include "TChain.h"
#include "TStopwatch.h"

#include "manager/MaCh3Logger.h"

#include "samplePDF/Structs.h"
#include "manager/YamlHelper.h"


namespace MaCh3Utils {
  /// @brief KS: Prints welcome message with MaCh3 logo
  void MaCh3Welcome();
  /// @brief KS: Find out more about operational system
  void GetOSInfo();
  /// @brief KS: Check what CPU you are using
  void GetCPUInfo();
  /// @brief KS: Check what GPU you are using
  void GetGPUInfo();
  /// @brief KS: Convoluted code to grab output from terminal to string
  std::string TerminalToString(const char* cmd);
  /// @brief KS: Check what CPU you are using
  void EstimateDataTransferRate(TChain* chain, const int entry);
  /// @brief KS: Find out about Disk usage
  void GetDiskUsage();
  /// @brief KS: Simply print progress bar
  void PrintProgressBar(const int Done, const int All);
  /// @brief KS: Get version of MaCh3
  std::string GetMaCh3Version();
  /// @brief CW: Get info like RAM
  int getValue(std::string Type);
  /// @brief CW: Get memory, which is probably silly
  int parseLine(const std::string& line);
  /// @brief KS: Print Yaml config using logger
  void PrintConfig(const YAML::Node& node);
}
