#pragma once

// spdlog Includes
#ifndef USE_FPGA
  #include "spdlog/spdlog.h"
#endif

// C++ Includes
#include <iostream>
#include <sstream>
#include <functional>
#include <string>

//KS: Based on this https://github.com/gabime/spdlog/blob/a2b4262090fd3f005c2315dcb5be2f0f1774a005/include/spdlog/spdlog.h#L284

#ifndef USE_FPGA
  #define MACH3LOG_TRACE SPDLOG_TRACE
  #define MACH3LOG_DEBUG SPDLOG_DEBUG
  #define MACH3LOG_INFO SPDLOG_INFO
  #define MACH3LOG_WARN SPDLOG_WARN
  #define MACH3LOG_ERROR SPDLOG_ERROR
  #define MACH3LOG_CRITICAL SPDLOG_CRITICAL
  #define MACH3LOG_OFF SPDLOG_OFF
#else
  #define MACH3LOG_TRACE
  #define MACH3LOG_DEBUG
  #define MACH3LOG_INFO
  #define MACH3LOG_WARN
  #define MACH3LOG_ERROR
  #define MACH3LOG_CRITICAL
  #define MACH3LOG_OFF
#endif

/// @brief Set messaging format of the logger
inline void SetMaCh3LoggerFormat()
{
  #ifndef USE_FPGA
    //KS: %H for hour, %M for minute, %S for second, [%s:%#] for class and line
    //For documentation see https://github.com/gabime/spdlog/wiki/3.-Custom-formatting
    #ifdef DEBUG
    //spdlog::set_pattern("[%H:%M:%S][%s:%#][%^%l%$] %v");
    spdlog::set_pattern("[%s:%#][%^%l%$] %v");
    #else
    //spdlog::set_pattern("[%H:%M:%S][%s][%^%l%$] %v");
    spdlog::set_pattern("[%s][%^%l%$] %v");
    #endif
  #endif
}

/// @brief KS: This is bit convoluted but this is to allow redirecting cout and errors from external library into MaCh3 logger format
/// @tparam Func The type of the function to be called, which outputs to stdout and stderr.
/// @tparam LogFunc The type of the logging function, typically a lambda that formats and logs messages.
/// @tparam Args The types of the arguments to be passed to `func`.
/// @param LibName The name of the library or component from which the log message originates.
/// @param logFunction The lambda function responsible for formatting and logging messages.
///                    It should accept a single const std::string& parameter.
/// @param func The function to be called, whose output needs to be captured and logged.
/// @param args The arguments to be passed to `func`.
/// @code
/// void ExampleFunc(int x, const std::string& str) {
///     std::cout << "Output from exampleFunc: " << x << ", " << str << "\n";
///     std::cerr << "Error from exampleFunc: " << x << ", " << str << "\n";
/// }
///
/// int main() {
///     LoggerPrint("ExampleLib", [](const std::string& message) { MACH3LOG_INFO("{}", message); },
///                   static_cast<void(*)(int, const std::string&)>(ExampleFunc), 666, "Number of the BEAST");
///     return 0;
/// }
/// @endcode
/// @note second argument is lambda fucniton whih convers mach3 so template works, it is bit faff...
///  This approach allows seamless integration of `func` into an existing logging mechanism without modifying `func` itself
template <typename Func, typename LogFunc, typename... Args>
void LoggerPrint(const std::string& LibName, LogFunc logFunction, Func&& func, Args&&... args)
{
  // Create a stringstream to capture the output
  std::stringstream sss;
  // Save the original stream buffers
  std::streambuf* coutBuf = std::cout.rdbuf();
  std::streambuf* cerrBuf = std::cerr.rdbuf();

  // Redirect std::cout and std::cerr to the stringstream buffer
  std::cout.rdbuf(sss.rdbuf());
  std::cerr.rdbuf(sss.rdbuf());

  // Call the provided function
  func(std::forward<Args>(args)...);

  // Restore the original stream buffers
  std::cout.rdbuf(coutBuf);
  std::cerr.rdbuf(cerrBuf);

  std::string line;
  while (std::getline(sss, line))
  {
    #ifndef USE_FPGA
      auto formatted_message = fmt::format("[{}] {}", LibName, line);
      logFunction(formatted_message);
    #else
      logFunction(line);
    #endif
  }
}
