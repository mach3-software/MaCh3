#pragma once

// MaCh3 includes
#include "Manager/Manager.h"
#include "Manager/Random.h"

/// @brief Helper to check if files exist or not
inline std::string file_exists(std::string filename) {
  std::ifstream infile(filename.c_str());
  if (!infile.good()) {
    MACH3LOG_ERROR("File {} does not exist", filename);
    MACH3LOG_ERROR("Please try again");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return filename;
}


namespace M3 {
/// @brief Opens a ROOT file with the given name and mode.
///
/// This function wraps ROOT’s `TFile` constructor and checks whether the file was opened
/// successfully and give some useful debugging information
///
/// @param Name The name or path of the file to open.
/// @param Type The file open mode (e.g., "READ", "RECREATE", "UPDATE").
/// @param File The name of the file where the exception occurred.
/// @param Line The line number where the exception occurred.
/// @return Pointer to the opened ROOT file.
TFile* Open(const std::string& Name, const std::string& Type, const std::string& File, const int Line);
} //end M3
