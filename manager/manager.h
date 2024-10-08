#pragma once

// MaCh3 Includes
#include "manager/MaCh3Logger.h"
#include "samplePDF/Structs.h"
#include "manager/YamlHelper.h"
#include "manager/Monitor.h"
#include "manager/MaCh3Exception.h"
#include "manager/MaCh3Modes.h"

//KS: Joy of forward declaration https://gieseanw.wordpress.com/2018/02/25/the-joys-of-forward-declarations-results-from-the-real-world/
class TFile;

/// @brief The manager class is responsible for managing configurations and settings.
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/01.-Manager-and-config-handling).
class manager {
public:
  /// @brief Constructs a manager object with the specified file name.
  /// @param filename The name of the configuration file.
  explicit manager(std::string const &filename);
  /// @brief Destroys the manager object.
  virtual ~manager();

  /// @brief Add manager useful information's to TFile, in most cases to Fitter
  /// @param OutputFile The ROOT TFile to which the information will be added.
  void SaveSettings(TFile* const OutputFile);

  /// @brief Print currently used config
  void Print();

  /// @brief Get likelihood type defined in the config
  inline int GetMCStatLLH(){return mc_stat_llh;}
  /// @brief Return name of config
  inline std::string GetFileName(){return FileName;}
  /// @brief Return config
  inline YAML::Node const &raw(){return config;}
  /// @brief Return pointer to MaCh3 modes
  MaCh3Modes* GetMaCh3Modes() const { return Modes; }
  /// @brief Get class name
  inline std::string GetName()const {return "Manager";};

  /// @brief Overrides the configuration settings based on provided arguments.
  ///
  /// This function allows you to set configuration options for the manager.
  /// It accepts either two or three string arguments:
  /// - For two arguments, the first argument is a key, and the second is the value.
  /// - For three arguments, the first two arguments are keys, and the third is the value.
  ///
  /// @param args The arguments to override the configuration.
  ///             - When two arguments are provided, they represent the key and value, respectively.
  ///             - When three arguments are provided, they represent two keys and a value.
  ///
  /// @note Example usage:
  /// @code
  /// FitManager->OverrideConfig("General", "OutputFile", "Wooimbouttamakeanameformyselfere.root");
  /// @endcode
  template <typename... Args>
  void OverrideConfig(Args... args) {
    static_assert(sizeof...(args) == 2 || sizeof...(args) == 3,
                  "OverrideConfig accepts either 2 or 3 arguments.");

    auto args_tuple = std::make_tuple(args...); // Create a tuple from the parameter pack

    if constexpr (sizeof...(args) == 2) {
      std::string blarb1 = std::get<0>(args_tuple); // First argument
      std::string result = std::get<1>(args_tuple); // Second argument

      config[blarb1] = result;
    }
    else if constexpr (sizeof...(args) == 3) {
      std::string blarb1 = std::get<0>(args_tuple); // First argument
      std::string blarb2 = std::get<1>(args_tuple); // Second argument
      std::string result = std::get<2>(args_tuple); // Third argument

      config[blarb1][blarb2] = result;
    }
  }
private:
  /// The YAML node containing the configuration data.
  YAML::Node config;
  /// The name of the configuration file.
  std::string FileName;
  /// The likelihood type defined in the configuration.
  int mc_stat_llh;
  /// MaCh3 Modes
  MaCh3Modes* Modes;
};
