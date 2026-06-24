#include "PlottingManager.h"

namespace MaCh3Plotting {
// this is the constructor using the default plotting config file
PlottingManager::PlottingManager() {
  // set config file name
  _configFileName = DEFAULT_PLOTTING_CONFIG;
}

// this is the constructor with user specified config
PlottingManager::PlottingManager(const std::string &PlottingConfigName) {
  // parse the config file
  _configFileName = PlottingConfigName;
}

/// - Read the plotting config file
/// - Instantiate a StyleManager using the style config file specified in the plotting config file,
/// or if none was provided, using DEFAULT_STYLE_CONFIG.
/// - Instantiate an InputManager using the translation config file specified in the plotting config
/// file, or if none was provided, using DEFAULT_TRANSLATION_CONFIG.
/// - Add all files specified in this PlottingManagers _fileNames vector to the new InputManager that
/// was just created
/// @warning This should always be called *After* parseInputs() unless you are
/// manually specifying all input file names and config file names in your drawing application.
void PlottingManager::initialise() {
  /// @todo should add some kind of validataConfigs() method to got throught all of the specified
  /// config files and make sure that all provided options are valid and all necessary options are provided
  /// as it can be pretty annoying and difficult to identify what's going wrong when yaml just fails
  /// to find an option at runtime
  _plottingConfig = M3OpenConfig(_configFileName);

  MACH3LOG_DEBUG("Initialising PlottingManager with plotting congif {}", _configFileName);
  // read options from the config
  YAML::Node managerOptions = _plottingConfig["ManagerOptions"];

  std::string translationConfig = managerOptions["translationConfig"].as<std::string>();
  if (translationConfig == "")
  {
    translationConfig = DEFAULT_TRANSLATION_CONFIG;
    MACH3LOG_DEBUG("PlottingManager: Using default translation config file name");
  }
  
  MACH3LOG_DEBUG("PlottingManager: Using translation config file: {}", translationConfig);

  std::string styleConfig = managerOptions["styleConfig"].as<std::string>();
  if (styleConfig == "")
  {
    styleConfig = DEFAULT_STYLE_CONFIG;
    MACH3LOG_DEBUG("PlottingManager: Using default style config file name");
  }

  MACH3LOG_DEBUG("PlottingManager: Using style config file: {}", styleConfig);

  // create the StyleManager
  _styleMan = std::make_unique<StyleManager>(styleConfig);

  // create the InputManager and add all the files to it
  _inputMan = std::make_unique<InputManager>(translationConfig);
  for (std::string fileName : _fileNames)
  {
    _inputMan->addFile(fileName);
  }
}

/// Takes in c style command line arguments and parses particular general ones that are useful for a
/// wide variety of plotting scripts.
///
/// *Please also keep these options in mind and use them in your own plotting scripts*.
///
/// General cmd line interface for a plotting script will look like:
/// @code
/// ./ScriptName [options] [FileName1 FileName2 FileName3 ...]
///
/// [options]
///   -o <_outputName> name of the file to output plots to
///   -l <FileLabelStr> string of labels for each file in the format Label1;Label2;... to be parsed
///   by parseFileLabels(), number of labels should match exactly the number of files -c
///   <PlottingConfig> specify a plotting config file to be used in place of the default one
///   when initialising this manager -d <DrawOptions> extra draw options to be passed when making
///   plots, see https://root.cern/doc/master/classTHistPainter.html#HP01a for examples
///
///   -s where possible, will try to split plots made by sample
///   -r where possible, make ratio plots comparing the different input files
///   -g where possible, draw a grid on the produced plots for improved readability
///
///   -h print the help message for the current executable and exit
///
/// unnamed options are treated as files to be processed
/// @endcode
/// @todo make this able to return any un-parsed arguments so that user can specify their own
/// arguments for use in their plotting scripts
void PlottingManager::parseInputs(int argc, char * const *argv) {
  // parse the inputs
  int c;
  while ((c = getopt(argc, argv, "o:l:d:c:srgh")) != -1)
  {
    if (c < 0)
      break;
    switch (c)
    {
    case 'o': {
      setOutFileName(optarg);
      break;
    }
    case 's': {
      _splitBySample = true;
      break;
    }
    case 'r': {
      _plotRatios = true;
      break;
    }
    case 'g': {
      _drawGrid = true;
      break;
    }
    case 'h': {
      usage();
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    case 'l': {
      parseFileLabels(optarg, _fileLabels);
      MACH3LOG_INFO("Specified file labels: ");
      for (std::string label : _fileLabels)
        MACH3LOG_INFO("  {}", label);
      break;
    }
    case 'c': {
      _configFileName = optarg;
      break;
    }
    case 'd': {
      _extraDrawOptions = optarg;
      break;
    }
    }
  }

  MACH3LOG_INFO("Input files provided");
  // optind is for the extra arguments that are not parsed by the program
  for (; optind < argc; optind++)
  {
    _fileNames.push_back(argv[optind]);
    MACH3LOG_INFO(argv[optind]);
    _defaultFileLabels.push_back(argv[optind]);
  }
  
  if (_splitBySample)
    MACH3LOG_INFO("Splitting by sample");

  if (_plotRatios && _fileNames.size() == 0)
  {
    MACH3LOG_ERROR("you specified -r <_plotRatios> = true but didnt specify any files to compare against, was this a mistake?");
  }

  if (_fileLabels.size() == 0)
  {
    MACH3LOG_INFO("No file labels specified, will just use the file names");
    _fileLabels = _defaultFileLabels;
  }

  if ((_fileLabels.size() != 0) && (_fileLabels.size() != _fileNames.size()))
  {
    MACH3LOG_ERROR("hmmm, you gave me {} labels but {} files", _fileLabels.size(), _fileNames.size());
  }

  initialise();
}

/// @todo Would be good to add functionality to this to allow user to add their own options.
/// @todo The way I can think to do this would be have fn addUserOption() to set the options,
/// defining the cmd line option (e.g. -x), the name of the option, and maybe some description of
/// the option to be used in the help message
/// @todo can then store these options in some map or
/// something to be retrieved later by getUserOption()
void PlottingManager::addUserOption() {
  /// @todo Implement this.
}

std::string PlottingManager::getUserOption(std::string option) {
  /// @todo Implement this.
  (void) option;
  return "";
}

void PlottingManager::usage() {
  /// @todo Implement this.
  /// @todo could add some function to allow user to specify the help message for their particular
  /// script, then auto generate what the cmd line syntax looks like based on user specified
  /// options?
}

/// Will check the provided saveName for file extensions, if it is one of .pdf or .eps, then just
/// use the provided string as the full output name. If no file extension is specified, append .pdf
/// so plots will be saved as pdf. if some other file extension is specified, replace with .pdf as
/// only .pdf and .eps support printing multiple plots to one file in root.
void PlottingManager::setOutFileName(const std::string& saveName) {
  if (saveName.find(".") == std::string::npos)
  {
    _outputName = saveName + ".pdf";
    return;
  }

  std::string ext = saveName;
  // if there are .'s then remove everything before them until we only have the file extension left
  while (ext.find(".") != std::string::npos)
    ext.erase(0, ext.find(".") + 1);
  if (ext == "pdf" || ext == "ps" || ext == "eps")
  {
    _outputName = saveName;
    return;
  }

  MACH3LOG_WARN("file extension '{}' that you provided doesnt support multiple plots in one file", ext);
  MACH3LOG_WARN("should be one of .pdf, .eps .ps, will use pdf");
  _outputName = saveName + ".pdf";
}

/// Output file name, including the file extension will be returned, but with specified suffix after
/// the name but before the extension. This is useful for e.g. saving multiple LLH scan types to
/// separate files: can specify suffix "_PriotLLH" will return OutputName_PriorLLH.ext
/// @todo Make this support .root files too
const std::string PlottingManager::getOutputName(const std::string &suffix) {
  std::string ext = _outputName;
  std::string name = _outputName;

  size_t dotPos = 0;
  while (ext.find(".") != std::string::npos)
  {
    dotPos += ext.find(".");
    ext.erase(0, ext.find(".") + 1);
  }

  name.erase(dotPos, ext.size() + 1);
  return name + suffix + "." + ext;
}

/// This is used by e.g. getOption() so the PlottingManager knows where to look for the option in
/// the plotting config file.
void PlottingManager::setExec(const std::string& execName) {
  MACH3LOG_DEBUG("Setting internal exec name to {}", execName);
  _execOptions = _plottingConfig[execName];
  if (_outputName == "Plot.pdf")
    setOutFileName(getOption<std::string>("defaultOutputName"));
}

void PlottingManager::parseFileLabels(std::string labelString, std::vector<std::string> &labelVec) {
  // take in a string defining labels of the form "label1;label2;...;labelN" and parse it into a
  // vector containing the labels

  size_t end = labelString.find(";");
  while (end != std::string::npos)
  { // Loop until no delimiter is left in the string.
    labelVec.push_back(labelString.substr(0, end));
    labelString.erase(labelString.begin(), labelString.begin() + end + 1);
    end = labelString.find(";");
  }
  labelVec.push_back(labelString.substr(0, end));
}
} // namespace MaCh3Plotting
