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
/// @author Ed Atkin
/// @author Kamil Skwarczynski
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
  int GetMCStatLLH();
  /// @brief Return name of config
  inline std::string GetFileName(){return FileName;}
  /// @brief Return config
  inline YAML::Node const &raw(){return config;}
  /// @brief Return pointer to MaCh3 modes
  MaCh3Modes* GetMaCh3Modes() const { return Modes; }
  /// @brief Get class name
  inline std::string GetName()const {return "Manager";};

  /// @brief Overrides the configuration settings based on provided arguments.
  /// @note Example usage:
  /// @code
  /// FitManager->OverrideSettings("General", "OutputFile", "Wooimbouttamakeanameformyselfere.root");
  /// @endcode
  template <typename... Args>
  void OverrideSettings(Args&&... args) {
    OverrideConfig(config, std::forward<Args>(args)...);
  }

private:
  /// The YAML node containing the configuration data.
  YAML::Node config;
  /// The name of the configuration file.
  std::string FileName;
  /// MaCh3 Modes
  MaCh3Modes* Modes;
};
