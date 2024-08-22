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
#include "manager/YamlHelper.h"
#include "manager/MaCh3Logger.h"

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
      "/plotting/StyleConfig.yaml"; //!< The default style config config to be used by PlottingManager
                                    //!< instances when instantiating the StyleManager.
  const std::string DEFAULT_PLOTTING_CONFIG =
      std::string(std::getenv("MACH3")) +
      "/plotting/PlottingConfig.yaml"; //!< The default plotting config config to use when
                                       //!< instantiating new PlottingManager objects.

  /// @brief Construct a new PlottingManager using default plottingConfig config.
  /// @return Constructed PlottingManager instance.
  PlottingManager();

  /// @brief Construct a new PlottingManager using specified plottingConfig config.
  /// @param PlottingConfigName The config file file defining executables specific options, and other minor
  /// manager related options.
  /// @return Constructed PlottingManager instance.
  PlottingManager(std::string PlottingConfigName);

  /// @brief initalise this PlottingManager.
  void Initialise();

  ~PlottingManager() {
    delete styleMan;
    delete inputMan;
  }

  /// @brief Parse command line arguments.
  /// @param argc The number of command line arguments.
  /// @param argv The arguments themselves.
  void ParseInputs(int argc, char **argv);

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
  void SetExec(std::string execName);

  /// @brief Get a specific option from the config for this executable.
  /// @tparam T the type of parameter expected for this option, e.g. std::string.
  /// @param option The option that you want from the config.
  /// @return The specified value for the option.
  template <typename T> T GetOption(std::string option) {
    return _execOptions[option].as<T>();
  }
  YAML::Node GetOption(std::string option) 
  { 
    return _execOptions[option]; 
  }

  // ############# getters ##############
  /// @name General getters
  /// @{
  const std::string GetFileName(int i) 
	{
		return FileNames[i];
	}

  const std::string GetFileLabel(int i) 
	{
		return FileLabels[i];
	}

  const std::string GetDrawOptions() 
	{
		return extraDrawOptions;
	}

  /// @brief Get the straight up output file name with no bells or whistles, just the file
  /// extension.
  /// @return The straight up output file name.
  const std::string GetOutputName() 
	{
		return OutputName;
	}

  /// @brief Get the output name but can specify a siffix to add to the name, before the file
  /// extension.
  /// @param suffix The suffix to add to the file name.
  /// @return Output file name with suffix added before the extension.
  const std::string GetOutputName(std::string suffix);

  const std::vector<std::string> GetFileNames() 
	{
		return FileNames;
	}

  const std::vector<std::string> GetFileLabels() 
	{
		return FileLabels;
	}


  const int GetNFiles() 
	{
		return (int)FileNames.size();
	}


  const bool GetSplitBySample() 
	{
		return splitBySample;
	}

  const bool GetPlotRatios() 
	{
		return plotRatios;
	}

  const bool GetDrawGrid() 
	{
		return drawGrid;
	}

  /// @}

  // for managers contained in this manager
  /// @brief Get the StyleManager contained within this PlottingManager, for doing style related
  /// things.
  const StyleManager *Style() 
	{
		return styleMan;
	}

  /// @brief Get the InputManager contained within this PlottingManager, for doing input related
  /// things.
  const InputManager *Input() 
	{
		return inputMan;
	}


private:
  // name of the config file to read configs from
  std::string _configFileName = std::string(std::getenv("MACH3")) + "/plotting/PlottingConfig.yaml";
  // the parsed config
  YAML::Node _plottingConfig;
  // the config object holding the options for the current executable
  YAML::Node _execOptions;

  // input file names and labels
  std::vector<std::string> FileNames;
  std::vector<std::string> FileLabels;
  std::vector<std::string> FileLabels_default;

  // other string options
  std::string OutputName = "Plot.pdf";
  std::string extraDrawOptions = "";

  // Generic plotting options
  bool splitBySample = false;
  bool plotRatios = false;
  bool drawGrid = false;

  // other Manager objects
  StyleManager *styleMan = NULL;
  InputManager *inputMan = NULL;
};
} // namespace MaCh3Plotting
