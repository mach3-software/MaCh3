#pragma once

// MaCh3 Includes
#include "Manager/MaCh3Logger.h"
#include "Manager/YamlHelper.h"
#include "Manager/MaCh3Exception.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT Things
#include "TColor.h"
#include "TH1.h"
#include "TStyle.h"
_MaCh3_Safe_Include_End_ //}


namespace MaCh3Plotting {
/// @author Ewan Miller

/// @brief EW: provides centralized styling utilities for plots, including name prettification and style application
/// @ingroup MaCh3Plotting
class StyleManager {
public:
  /// @brief Constructor
  /// @param configName The style config to read from
  StyleManager(std::string configName);
  
  // NO COPYING!
  StyleManager( const StyleManager& ) = delete;
  StyleManager( StyleManager&& ) = default;

  ~StyleManager(){
    MACH3LOG_DEBUG("##### Deleting StyleManager Instance #####");
  }

  /// @brief Convert hideous and vulgar internal parameter name into a beautiful presentable name
  /// @details The pretty parameter names should be specified in the style config file
  /// @param origName The "internal" name used to uniquely identify the parameter inside the
  /// plotting code
  /// @return A beautiful formatted name that can be used in plots
  inline std::string prettifyParamName(const std::string &origName) const {
    return prettifyName(origName, "parameters");
  };

  /// @brief Convert hideous and vulgar internal sample name into a beautiful presentable name
  /// @details The pretty sample names should be specified in the style config file
  /// @param origName The "internal" name used to uniquely identify the sample inside the plotting
  /// code
  /// @return A beautiful formatted name that can be used in plots
  inline std::string prettifySampleName(const std::string &origName) const {
    return prettifyName(origName, "samples");
  };

  /// @brief Convert hideous and vulgar internal kinematics name into a beautiful presentable name
  /// @details The pretty kinematics names should be specified in the style config file
  /// @param origName The "internal" name used to uniquely identify the kinematics inside the plotting
  /// code
  /// @return A beautiful formatted name that can be used in plots
  inline std::string prettifyKinematicName(const std::string &origName) const {
    return prettifyName(origName, "kinematics");
  };

  // style setting options
  /// @brief Set the root colour palette to one of the default root pallettes as defined in (root
  /// docs)[https://root.cern.ch/doc/master/classTColor.html#C05]
  /// @param rootPlotStyle The index of the palette as defined by root
  void setPalette(const int rootPlotStyle) const;

  /// @brief Set the root colour palette to one of the ones defined in the style config
  /// @param rootPlotStyle The name of the palette you want to use, should be the same as it appears
  /// in the style config
  void setPalette(const std::string& configStyleName) const;

  /// @brief Set the style of a TH1 to one of the styles defined in the style config
  /// @param hist The TH1 that you wish to modify
  /// @param styleName The name of the style you want to use, as it appears in the config file
  void setTH1Style(TH1 *hist, const std::string& styleName) const;

private:
  YAML::Node _styleConfig;

  // helper to basically just read fancy names from config  
  std::string prettifyName(const std::string &origName, const std::string &nameType) const;

};

} // namespace MaCh3Plotting
