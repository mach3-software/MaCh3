#pragma once

// C++ includes
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>

// ROOT includes
#include "TColor.h"
#include "TH1.h"
#include "TStyle.h"

// MaCh3 Includes
#include "manager/MaCh3Logger.h"
#include "manager/YamlHelper.h"
#include "manager/MaCh3Exception.h"

// Other MaCh3Plotting stuff
#include "inputManager.h"
#include "styleManager.h"

namespace MaCh3Plotting {
/// @brief The main class to be used in plotting scripts.
/// @details When it comes to plotting, this guys in charge, the main man, the head honcho, the big
/// cheese. If it's a plot you need, this is the guy you call. You just call him in your scripts and
/// let him and his goons StyleManager, and InputManager worry about the details, no questions
/// asked. This class handles all the command line inputs and manages the other managers. Intention
/// with this class is that you should be able to create one from your custom plotting script, then
/// from that, be able to get all manner of inputs from the InputManager contained within, and set
/// any style options you like from the StyleManager. Also a hopefully not too distant dream is to
/// wrap this up in python so it is usable in .py scripts to take advantage of nice existing
/// plotting libraries for e.g. MCMC plotting.
/// @author Ewan Miller
class PlottingManager {
public:
  // EM: cant make these static as std::getenv("MACH3") not known at compile time
  const std::string DEFAULT_TRANSLATION_CONFIG =
      std::string(std::getenv("MACH3")) +
      "/plotting/universalTranslator.yaml"; //!< The default translation config to be used by
                                            //!< PlottingManager instances when instantiating the
                                            //!< InputManager.
  const std::string DEFAULT_STYLE_CONFIG =
      std::string(std::getenv("MACH3")) +
      "/plotting/StyleConfig.yaml"; //!< The default style config config to be used by
                                    //!< PlottingManager instances when instantiating the
                                    //!< StyleManager.
  const std::string DEFAULT_PLOTTING_CONFIG =
      std::string(std::getenv("MACH3")) +
      "/plotting/PlottingConfig.yaml"; //!< The default plotting config config to use when
                                       //!< instantiating new PlottingManager objects.

  /// @brief Construct a new PlottingManager using default plottingConfig config.
  /// @return Constructed PlottingManager instance.
  PlottingManager();

  /// @brief Construct a new PlottingManager using specified plottingConfig config.
  /// @param PlottingConfigName The config file file defining executables specific options, and
  /// other minor manager related options.
  /// @return Constructed PlottingManager instance.
  PlottingManager(const std::string &PlottingConfigName);

  /// @brief initalise this PlottingManager.
  void initialise();

  ~PlottingManager() {
    MACH3LOG_DEBUG("##### Deleting PlottingManager Instance #####");
  }

  /// @brief Parse command line arguments.
  /// @param argc The number of command line arguments.
  /// @param argv The arguments themselves.
  void parseInputs(int argc, char* const *argv);


  /// @brief Parse vector of command line arguments.
  /// @details This mainly just exists for the sake of the python binding.
  /// @param argv The arguments to parse.
  inline void parseInputsVec(std::vector<std::string> argv) {
    std::vector<char *> charVec;
    MACH3LOG_DEBUG("Parsing Inputs :: was given vector:");
    for( const std::string &arg : argv ) 
    {
      charVec.push_back( const_cast<char *>(arg.c_str()) );
      MACH3LOG_DEBUG("  - {}", arg );
    }
    parseInputs(int(argv.size()), charVec.data());
  }

  /// @brief Describe an option you want to add to the PlottingManager which can be read in from the
  /// command line and retrieved later with getOption(). This should be done before calling
  /// parseInputs().
  void addUserOption();

  /// @brief Retrieve a command line option you specified using addOption.
  std::string getUserOption(std::string option);

  /// @brief Print a usage message for the current executable.
  void usage();

  /// @brief Parse string of labels into a vector of strings.
  /// @param labelString string of labels of the form "label1;label2;...".
  /// @param labelVec vector that the individual label strings will be placed into.
  void parseFileLabels(std::string labelString, std::vector<std::string> &labelVec);

  /// @brief Parse and set the output file name, if extension specified, check its one root
  /// supports, if not, default to pdf.
  /// @param fileName the name of the output file.
  void setOutFileName(std::string fileName);

  /// @brief Internally set the name of the executable that manager is being used in.
  /// @param execName Name of the current executable, will also need to be defined in the plotting
  /// config file.
  void setExec(std::string execName);

  /// @brief Get a specific option from the config for this executable.
  /// @tparam T the type of parameter expected for this option, e.g. std::string.
  /// @param option The option that you want from the config.
  /// @return The specified value for the option.
  template <typename T> T getOption(std::string option) { return _execOptions[option].as<T>(); }
  YAML::Node getOption(std::string option) { return _execOptions[option]; }

  // ############# getters ##############
  /// @name General getters
  /// @{
  const std::string getFileName(int i) { return _fileNames[i]; }

  const std::string getFileLabel(int i) { return _fileLabels[i]; }

  const std::string getDrawOptions() { return _extraDrawOptions; }

  /// @brief Get the straight up output file name with no bells or whistles, just the file
  /// extension.
  /// @return The straight up output file name.
  const std::string getOutputName() { return _outputName; }

  /// @brief Get the output name but can specify a siffix to add to the name, before the file
  /// extension.
  /// @param suffix The suffix to add to the file name.
  /// @return Output file name with suffix added before the extension.
  const std::string getOutputName(const std::string &suffix);

  const std::vector<std::string> getFileNames() { return _fileNames; }

  const std::vector<std::string> getFileLabels() { return _fileLabels; }

  size_t getNFiles() { return _fileNames.size(); }

  bool getSplitBySample() { return _splitBySample; }

  bool getPlotRatios() { return _plotRatios; }

  bool getDrawGrid() { return _drawGrid; }

  /// @}

  // for managers contained in this manager
  /// @brief Get the StyleManager contained within this PlottingManager, for doing style related
  /// things.
  const StyleManager &style() { return *_styleMan; }

  /// @brief Get the InputManager contained within this PlottingManager, for doing input related
  /// things.
  const InputManager &input() { return *_inputMan; }

private:
  // name of the config file to read configs from
  std::string _configFileName = std::string(std::getenv("MACH3")) + "/plotting/PlottingConfig.yaml";
  // the parsed config
  YAML::Node _plottingConfig;
  // the config object holding the options for the current executable
  YAML::Node _execOptions;

  // input file names and labels
  std::vector<std::string> _fileNames;
  std::vector<std::string> _fileLabels;
  std::vector<std::string> _defaultFileLabels;

  // other string options
  std::string _outputName = "Plot.pdf";
  std::string _extraDrawOptions = "";

  // Generic plotting options
  bool _splitBySample = false;
  bool _plotRatios = false;
  bool _drawGrid = false;

  // other Manager objects
  std::unique_ptr<StyleManager> _styleMan;
  std::unique_ptr<InputManager> _inputMan;
};
} // namespace MaCh3Plotting
