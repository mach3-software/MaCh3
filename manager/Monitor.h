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

  /**
   * @brief KS: Prints welcome message with MaCh3 logo
   */
  void MaCh3Welcome();
  //KS: Check what CPU you are using
  void GetCPUInfo();
  //KS: Convoluted code to grab output from terminal to string
  std::string TerminalToString(const char* cmd);
  //KS: Check what CPU you are using
  void EstimateDataTransferRate(TChain* chain, const int entry);
  //KS: Simply print progress bar
  void PrintProgressBar(const int Done, const int All);
  //CW: Get info like RAM
  int getValue(std::string Type);
  int parseLine(const std::string& line);
}
