#pragma once

// C++ Includes
#include <iostream>
#include <string>
#include <stdexcept>

// MaCh3 Includes
#include "manager/MaCh3Logger.h"

/// @brief Custom exception class for MaCh3 errors.
class MaCh3Exception : public std::exception {
public:
  /// @brief Constructs a MaCh3Exception object with the specified error message.
  /// @param File The name of the file where the exception occurred.
  /// @param Line The line number where the exception occurred.
  /// @param Message The error message describing the exception (optional).
  explicit MaCh3Exception(std::string File, int Line, std::string Message = "")
  {
    size_t lastSlashPos = File.find_last_of('/');
    std::string fileName = (lastSlashPos != std::string::npos) ? File.substr(lastSlashPos + 1) : File;

    errorMessage = ((Message.empty()) ? "Terminating MaCh3" : Message);
    // KS: Set logger format where we only have have "information type", line would be confusing
    #ifndef USE_FPGA
      spdlog::set_pattern("[%^%l%$] %v");
    #endif
    MACH3LOG_ERROR("Find me here: {}::{}", fileName, Line);
  }

  /// @brief Returns the error message associated with this exception.
  /// @return A pointer to the error message string.
  const char* what() const noexcept override {
    return errorMessage.c_str();
  }

private:
  /// The error message associated with this exception.
  std::string errorMessage;
};
