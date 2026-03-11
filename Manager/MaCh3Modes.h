#pragma once

// C++ Includes
#include <cmath>
#include <map>
#include <string>

//MaCh3 includes
#include "Manager/YamlHelper.h"
#include "Manager/Monitor.h"
#include "Manager/MaCh3Exception.h"

/// Enumerator of MaCh3Mode
using MaCh3Modes_t = int;

/// @brief KS: Class containing information for a single MaCh3Mode
struct MaCh3ModeInfo {
  /// Mode name
  std::string Name;
  /// Mode fancy name
  std::string FancyName;
  /// Mode color for plotting purposes
  int PlotColor;
  /// Mapping between mode and generator integers
  std::vector<int> GeneratorMaping;
  /// IsNC check
  bool IsNC;
  /// Spline suffix
  std::string SplineSuffix;
  /// @brief KS: Checks MaCh3 modes is associated with a given generator mode
  inline bool IsMode(const int GenMode) const {
    bool exists = std::find(GeneratorMaping.begin(), GeneratorMaping.end(), GenMode) != GeneratorMaping.end();
    return exists;
  }
};

/// @brief KS: Class describing MaCh3 modes used in the analysis, it is being initialised from config
///
/// The class provides a translation layer between generator-specific interaction
/// mode definitions and the analysis-level categories used inside MaCh3.
///
/// Different neutrino generators use different numbering and naming conventions.
/// Additionally, an analysis may wish to group or redefine these modes.
/// `MaCh3Modes` allows the user to define an experiment-specific mapping via
/// a YAML configuration file.
///
/// The class builds:
///  - A list of MaCh3 analysis modes
///  - A mapping from generator mode IDs → MaCh3 categories
///  - A protected fallback category (`UNKNOWN_BAD`) for unmapped values
///
/// @section MaCh3Modes_Config Configuration Example
/// Example YAML configuration used to initialise the modes:
///
/// @code{.yaml}
/// #########################################
/// Title: "T2K Interaction modes"
///
/// GeneratorName: "NEUT"
///
/// MaCh3Modes: [
///   "CCQE",
///   "CC1pipm",
///   "2p2h",
/// ]
///
/// CCQE:
///   Name: "CCQE"
///   GeneratorMaping: [1]      # 1 = CCQE
///
/// CC1pipm:
///   Name: "CC 1#pi^{#pm}"
///   GeneratorMaping: [11, 13] # 11 = CC1pi+1p, 13 = CC1pi+1n
///
/// 2p2h:
///   Name: "2p2h"
///   GeneratorMaping: [2]      # 2 = (Nieves) MEC / 2p2h
/// @endcode
///
/// @section MaCh3Modes_Output Example Initialisation Dump
/// After loading the configuration, MaCh3 produces a mapping table like:
///
/// @code{.text}
/// Printing MaCh3 Modes Called: T2K Interaction modes
/// ===============================================================
/// #   | Name     | FancyName           | NEUT Modes
/// ---------------------------------------------------------------
/// 0   | CCQE     | CCQE                | 1
/// 1   | CC1pipm  | CC 1#pi^{#pm}       | 11 13
/// 2   | 2p2h     | 2p2h                | 2
/// ===============================================================
///
/// NEUT Modes | Name
/// ---------------------------
/// 0          | UNKNOWN_BAD
/// 1          | CCQE
/// 2          | 2p2h
/// 3–10       | UNKNOWN_BAD
/// 11         | CC1pipm
/// 12         | UNKNOWN_BAD
/// 13         | CC1pipm
/// ===============================================================
/// @endcode
///
/// @section MaCh3Modes_Unknown Unknown Category Handling
/// Generators such as NEUT define many interaction codes that may not be used
/// in a given analysis. Any generator mode not explicitly mapped is assigned
/// to `UNKNOWN_BAD`.
///
/// This behaviour:
///  - Prevents segfaults from invalid lookups
///  - Allows sparse or partial mode definitions
///  - Provides a safe catch-all for unexpected generator values
///
/// @section MaCh3Modes_Usage Usage Examples
///
/// Access a specific mode:
/// @code{.cpp}
/// MaCh3Modes_t kMaCh3_CCQE = Modes->GetMode("CCQE");
/// @endcode
///
/// Loop over all defined modes:
/// @code{.cpp}
/// for (int j = 0; j < Modes->GetNModes() + 1; ++j) {
///   Modes->GetMaCh3ModeName(j);
/// }
/// @endcode
///
/// @note If you attempt to access a mode outside the defined range, the class
/// will return the `UNKNOWN_BAD` category instead of throwing.
///
///
/// @author Kamil Skwarczynski
/// @author Daniel Barrow
class MaCh3Modes {
 public:
  /// @brief KS: Initialise MaCh3 modes using path to config
  MaCh3Modes(std::string const &filename);
  /// @brief KS: Initialise MaCh3 modes using config
  MaCh3Modes(const YAML::Node& config);
  /// @brief KS: Empty destructor
  virtual ~MaCh3Modes(){};

  /// @brief KS: Print info about initialised modes
  void Print() const;

  /// @brief KS: Get number of modes, keep in mind actual number is +1 greater due to unknown category
  inline int GetNModes() const {return NModes;}
  /// @brief KS: Get mode number based on name, if mode not known you will get UNKNOWN_BAD
  /// @name of MaCh3 mode
  MaCh3Modes_t GetMode(const std::string& name) const;
  /// @brief KS: Get normal name of mode, if mode not known you will get UNKNOWN_BAD
  std::string GetMaCh3ModeName(const int Index) const;
  /// @brief KS: Get normal name of mode, if mode not known you will get UNKNOWN_BAD
  int GetMaCh3ModePlotColor(const int Index) const;
  /// @brief KS: Get fancy name of mode, if mode not known you will get UNKNOWN_BAD
  std::string GetMaCh3ModeFancyName(const int Index) const;
  /// @brief DB: Get IsNC (a check whether the given MaCh3 corresponds to a Neutral Current mode)
  bool IsMaCh3ModeNC(const int Index) const;
  /// @brief DB: Get binned spline mode suffic from MaCh3 Mode
  std::string GetSplineSuffixFromMaCh3Mode(const int Index);
  /// @brief KS: Get MaCh3 mode from generator mode
  MaCh3Modes_t GetModeFromGenerator(const int Index) const;
  /// @brief Get class name
  inline std::string GetName() const {return "MaCh3Modes";};
  /// @brief Return count of CC modes
  inline int GetNCCModes() const {return nCCModes;};

 private:
  /// @brief KS: Initialise MaCh3 modes based on provided config
  /// @param config YAML-based config used for class initialisation
  inline void Initialise(const YAML::Node& config);

  /// @brief KS: Make sure we don't have two modes with the same name
  /// @name of MaCh3 mode
  inline MaCh3Modes_t EnsureModeNameRegistered(std::string const &name);

  /// @brief KS: Add new mode
  inline void DeclareNewMode(std::string const &name,
                             std::string const &fancyname,
                             int PlotColor,
                             std::vector<int> const &GenMap,
                             const bool IsNC,
                             const std::string& SplineSuffix);

  /// @brief KS: Fill ModeMap
  inline void PrepareMap();

  /// KS: Handy map which helps find mode number based on string
  std::map<std::string, MaCh3Modes_t> Mode;
  /// KS: Main map storing info about used modes
  std::map<MaCh3Modes_t, MaCh3ModeInfo> fMode;
  /// KS: Handy map helping us find MaCh3 mode based on Generator mode value
  std::vector<int> ModeMap;

  /// KS: Name of loaded modes
  std::string Title;
  /// KS: Name of generator like NEUT, NuWro etc. this is to make stuff fancy
  std::string Generator;
  /// KS: Number of modes, keep in mind actual number is +1 greater due to unknown category
  int NModes;
  /// DB: Number of CC modes
  int nCCModes;
};
