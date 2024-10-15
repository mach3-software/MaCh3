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

// MaCh3 includes
#include "manager/MaCh3Logger.h"
#include "samplePDF/Structs.h"
#include "manager/MaCh3Exception.h"
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
  /// @param cmd The terminal command to execute.
  /// @return The output of the terminal command as a string.
  std::string TerminalToString(std::string cmd);
  /// @brief KS: Check what CPU you are using
  void EstimateDataTransferRate(TChain* chain, const int entry);
  /// @brief KS: Find out about Disk usage
  void GetDiskUsage();
  /// @brief KS: Simply print progress bar
  /// @param Done The number of tasks completed.
  /// @param All The total number of tasks.
  /// @details This function prints a progress bar to the terminal, indicating the percentage of tasks completed.
  void PrintProgressBar(const int Done, const int All);
  /// @brief KS: Get version of MaCh3
  /// @return The current MaCh3 version as a string.
  /// @details This function fetches and returns the version of the MaCh3 software being used.
  std::string GetMaCh3Version();
  /// @brief CW: Get info like RAM
  /// @param Type The type of system information to retrieve (e.g., RAM, CPU usage).
  /// @return The requested system information as an integer.
  /// @details This function fetches system information like RAM usage or other hardware details based on the specified type.
  int getValue(const std::string& Type);
  /// @brief CW: Get memory, which is probably silly
  /// @param line The line of text to parse.
  /// @return The extracted memory value as an integer.
  /// @details This function is used to parse a line of text and extract memory-related information.
  int parseLine(const std::string& line);
  /// @brief KS: Print Yaml config using logger
  /// @param node yaml config node
  void PrintConfig(const YAML::Node& node);
  /// @brief KS: Almost all MaCh3 executables have the same usage, prepare simple printer
  /// @param argc The number of command-line arguments.
  /// @param argv The array of command-line arguments.
  /// @details This function prints a simple usage guide for MaCh3 executables, typically called when incorrect arguments are passed.
  void MaCh3Usage(int argc, char **argv);
}
