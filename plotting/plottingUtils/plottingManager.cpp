#include "plottingManager.h"

namespace MaCh3Plotting {
// this is the constructor using the default plotting toml
PlottingManager::PlottingManager() {
  // set toml card name
  configFileName = DEFAULT_PLOTTING_TOML;
}

// this is the constructor with user specified toml
PlottingManager::PlottingManager(std::string PlottingTomlName) {
  // parse the toml card
  configFileName = PlottingTomlName;
}

/// - Read the plotting config toml
/// - Instantiate a StyleManager using the style config toml specified in the plotting config toml,
/// or if none was provided, using DEFAULT_STYLE_TOML.
/// - Instantiate an InputManager using the translation config toml specified in the plotting config
/// toml, or if none was provided, using DEFAULT_TRANSLATION_TOML.
/// - Add all files specified in this PlottingManagers FileNames vector to the new InputManager that
/// was just created 
/// @warning This should always be called *After* ParseInputs() unless you are
/// manually specifying all input file names and config toml names in your drawing application.
void PlottingManager::Initialise() {
  card_toml = toml_h::parse_card(configFileName);

  // read options from the config
  toml::value managerOptions = toml_h::find(card_toml, "ManagerOptions");

  std::string translationToml = toml_h::find<std::string>(managerOptions, "translationToml");
  if (translationToml == "")
    translationToml = DEFAULT_TRANSLATION_TOML;

  std::string styleToml = toml_h::find<std::string>(managerOptions, "styleToml");
  if (styleToml == "")
    styleToml = DEFAULT_STYLE_TOML;

  // create the StyleManager
  styleMan = new StyleManager(styleToml);

  // create the InputManager and add all the files to it
  inputMan = new InputManager(translationToml);
  for (std::string fileName : FileNames)
  {
    inputMan->addFile(fileName);
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
///   -o <OutputName> name of the file to output plots to
///   -l <FileLabelStr> string of labels for each file in the format Label1;Label2;... to be parsed
///   by parseFileLabels(), number of labels should match exactly the number of files -c
///   <PlottingConfig> specify a plotting config toml file to be used in place of the default one
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
void PlottingManager::ParseInputs(int argc, char **argv) {
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
      splitBySample = true;
      break;
    }
    case 'r': {
      plotRatios = true;
      break;
    }
    case 'g': {
      drawGrid = true;
      break;
    }
    case 'h': {
      usage();
      throw;
    }
    case 'l': {
      parseFileLabels(optarg, FileLabels);
      std::cout << "INFO: Specified file labels {";
      for (std::string label : FileLabels)
        std::cout << label << ", ";
      std::cout << "}" << std::endl;
      break;
    }
    case 'c': {
      configFileName = optarg;
      break;
    }
    case 'd': {
      extraDrawOptions = optarg;
      break;
    }
    }
  }

  std::cout << std::endl << "Input files provided: \n{";
  // optind is for the extra arguments that are not parsed by the program
  for (; optind < argc; optind++)
  {
    FileNames.push_back(argv[optind]);
    std::cout << argv[optind] << ", ";
    FileLabels_default.push_back(argv[optind]);
  }
  std::cout << "}" << std::endl << std::endl;

  if (splitBySample)
    std::cout << "Splitting by sample" << std::endl;

  if (plotRatios && FileNames.size() == 0)
  {
    std::cerr << "ERROR: you specified -r <plotRatios> = true but didnt specify any files to "
                 "compare against, was this a mistake?"
              << std::endl;
    throw;
  }

  if (FileLabels.size() == 0)
  {
    std::cout << "INFO: No file labels specified, will just use the file names" << std::endl;
    FileLabels = FileLabels_default;
  }

  if ((FileLabels.size() != 0) && (FileLabels.size() != FileNames.size()))
  {
    std::cerr << "ERROR: hmmm, you gave me " << FileLabels.size() << " labels but "
              << FileNames.size() << " files" << std::endl;
    std::cerr << "       that doesn\'t seem right to me, did you forget a file? or a label maybe?"
              << std::endl;
    throw;
  }

  Initialise();
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
  return "";
}

void PlottingManager::usage() {
  /// @todo Inplement this.
  /// @todo could add some function to allow user to specify the help message for their particular
  /// script, then auto generate what the cmd line syntax looks like based on user specified
  /// options?
}

/// Will check the provided saveName for file extenstions, if it is one of .pdf or .eps, then just
/// use the provided string as the full output name. If no file extension is specified, append .pdf
/// so plots will be saved as pdf. if some other file extension is specified, replace with .pdf as
/// only .pdf and .eps support printing multiple plots to one file in root.
void PlottingManager::setOutFileName(std::string saveName) {
  if (saveName.find(".") == std::string::npos)
  {
    OutputName = saveName + ".pdf";
    return;
  }

  std::string ext = saveName;
  // if there are .'s then remove everything before them until we only have the file extension left
  while (ext.find(".") != std::string::npos)
    ext.erase(0, ext.find(".") + 1);
  if (ext == "pdf" || ext == "ps" || ext == "eps")
  {
    OutputName = saveName;
    return;
  }

  std::cout << "WARNING: file extension " << ext
            << " that you provided doesnt support multiple plots in one file" << std::endl;
  std::cout << "        should be one of .pdf, .eps .ps, will use pdf" << std::endl;
  OutputName = saveName + ".pdf";
}

/// Output file name, including the file extension will be returned, but with specified suffix after
/// the name but before the extension. This is useful for e.g. saving multiple LLH scan types to
/// separate files: can specify suffix "_PriotLLH" will return OutputName_PriorLLH.ext 
/// @todo Make this support .root files too
const std::string PlottingManager::GetOutputName(std::string suffix) {
  std::string ext = std::string(OutputName);
  std::string name = std::string(OutputName);

  int dotPos = 0;
  while (ext.find(".") != std::string::npos)
  {
    dotPos += ext.find(".");
    ext.erase(0, ext.find(".") + 1);
  }

  name.erase(dotPos, ext.size() + 1);
  return name + suffix + "." + ext;
}

/// This is used by e.g. GetOption() so the PlottingManager knows where to look for the option in
/// the plotting config file.
void PlottingManager::SetExec(std::string execName) {
  ExecOptions = toml_h::find(card_toml, execName);
  if (OutputName == "Plot.pdf")
    setOutFileName(GetOption<std::string>("defaultOutputName"));
}

void PlottingManager::parseFileLabels(std::string labelString, std::vector<std::string> &labelVec) {
  // take in a string defining labels of the form "label1;label2;...;labelN" and parse it into a
  // vector containing the labels

  int end = labelString.find(";");
  while (end != -1)
  { // Loop until no delimiter is left in the string.
    labelVec.push_back(labelString.substr(0, end));
    labelString.erase(labelString.begin(), labelString.begin() + end + 1);
    end = labelString.find(";");
  }
  labelVec.push_back(labelString.substr(0, end));
}
} // namespace MaCh3Plotting