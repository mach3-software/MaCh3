#pragma once

// MaCh3 Includes
#include "Samples/SampleStructs.h"
#include "Manager/YamlHelper.h"
#include "Manager/Monitor.h"
#include "Manager/MaCh3Exception.h"

#ifdef MPIENABLED
#include <mpi.h>
#endif

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
  /// @brief Constructs a manager object with the specified YAML
  /// @note This is primarily used when initializing from a previous chain, allowing the creation of a manager instance based embedded YAML in that chain.
  /// @param ConfigNode Actual YAML config
  manager(const YAML::Node ConfigNode);

  /// @brief Destroys the manager object.
  virtual ~manager();

  /// @brief Add manager useful information's to TFile, in most cases to Fitter
  /// @param OutputFile The ROOT TFile to which the information will be added.
  void SaveSettings(TFile* const OutputFile) const;

  /// @brief Print currently used config
  void Print() const;

  /// @brief Get likelihood type defined in the config
  int GetMCStatLLH() const;
  /// @brief Return name of config
  inline std::string GetFileName() const {return FileName;}
  /// @brief Return config
  inline YAML::Node const &raw() const {return config;}
  /// @brief Get class name
  inline std::string GetName() const {return "Manager";};

  /// @brief Overrides the configuration settings based on provided arguments.
  /// @note Example usage:
  /// @code
  /// FitManager->OverrideSettings("General", "OutputFile", "Wooimbouttamakeanameformyselfere.root");
  /// @endcode
  template <typename... Args>
  void OverrideSettings(Args&&... args) {
    OverrideConfig(config, std::forward<Args>(args)...);
  }

  private :
    /// @brief Common inialiser for both constructors
    void Initialise();

  void InitialiseMPI();

  /// The YAML node containing the configuration data.
  YAML::Node config;
  /// The name of the configuration file.
  std::string FileName;
};
