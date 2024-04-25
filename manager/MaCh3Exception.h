#pragma once

#include <iostream>
#include <string>
#include <stdexcept>

#include "manager/MaCh3Logger.h"

/// @brief Custom exception class for MaCh3 errors.
class MaCh3Except : public std::exception {
public:
  /// @brief Constructs a MaCh3Except object with the specified error message.
  /// @param message The error message to associate with the exception.
  explicit MaCh3Except(const char* message) : errorMessage(message) {}

  /// @brief Returns the error message associated with this exception.
  /// @return A pointer to the error message string.
  const char* what() const noexcept override {
    return errorMessage;
  }

private:
  /// The error message associated with this exception.
  const char* errorMessage;
};

/// Custom wrapper around throw, to ensure we always get line and file where error appeared
inline void MaCh3Exception(std::string File, const int Line, std::string Message = "")
{
  size_t lastSlashPos = File.find_last_of('/');
  std::string fileName = (lastSlashPos != std::string::npos) ? File.substr(lastSlashPos + 1) : File;

  std::string ErrorMessage = (Message.empty()) ? "Terminating MaCh3" : Message;

  MACH3LOG_ERROR("Find me here:  {}::{}", fileName, Line);
  throw MaCh3Except(ErrorMessage.c_str());
}

