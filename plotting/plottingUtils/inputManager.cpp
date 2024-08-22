#include "inputManager.h"

// EM: will move this somewhere more sensible
#define __BAD_FLOAT__ -999.999

namespace MaCh3Plotting {
// this is the constructor with user specified translation config file
InputManager::InputManager(std::string translationConfigName) {
  // read the config file
  _translatorConfig = YAML::LoadFile(translationConfigName);

  // split up the parts of the config for easier access later
  _fitterSpecConfig = _translatorConfig["FitterSpec"];
  _parametersConfig = _translatorConfig["Parameters"];
  _samplesConfig = _translatorConfig["Samples"];

  // go through each of the parameters and sample arrays and keep track of which ones we've been
  // told about
  for (std::string parameter : _parametersConfig["Parameters"].as<std::vector<std::string>>())
  {
    knownParameters.push_back(parameter);
  }

  for (std::string sample : _samplesConfig["Samples"].as<std::vector<std::string>>())
  {
    knownSamples.push_back(sample);
  }
}

/// Open an input file and add to the manager, consists of:
///  - Initialise a new InputFile using the specified file
///  - Get info about the file, like what fitter it came from, what it contains e.g. LLH scans,
///  processed post fit parameters etc.
///  - Load up the data from the file (LLH scans etc.) and put them in a common format to be used by
///  plotting scripts
///  - Push back a pointer to the InputFile objcet to the vector of files known to this
///  InputManager.
void InputManager::addFile(std::string fileName) {
  InputFile *fileInfo = new InputFile(fileName);

  // EM: need to be done in this order since fillFileData needs to know info about the file, e.g.
  // fitter and what things are in it
  fillFileInfo(fileInfo);
  fillFileData(fileInfo);

  fileVec.push_back(fileInfo);
}

/// If printLevel is "summary", will loop through all the files known to this InputManager and call
/// InputFile::Summarise(). If printLevel is "dump", will print all parameters known to this
/// InputManager and then loop through all input files and call InputFile::Dump().
void InputManager::Print(std::string printLevel) const {
  std::cout << std::endl << "Printing contents of InputManager instance:" << std::endl;

  if (printLevel == "dump")
  {
    std::cout << "parameters known to this manager: " << std::endl;
    for (std::string param : knownParameters)
    {
      std::cout << "  " << param << std::endl;
    }
  }

  int fileCount = 0;
  for (InputFile *file : fileVec)
  {
    std::cout << " For file " << fileCount << std::endl;
    if (printLevel == "summary")
    {
      file->Summarise();
    } else if (printLevel == "dump")
    { file->Dump(); }
    fileCount++;
  }

  std::cout << std::endl;
}

const float InputManager::GetPostFitError(int fileNum, std::string paramName,
                                          std::string errorType) const {
  const InputFile *inputFileDef = fileVec[fileNum];

  // set default type if not specified
  if (errorType == "")
    errorType = inputFileDef->defaultErrorType;

  if (inputFileDef->postFitErrors.find(errorType) == inputFileDef->postFitErrors.end())
  {
    std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << std::endl;
    std::cerr << "  requested error type, "
              << "\"" << errorType << "\" does not exist in the specified file: " << std::endl;
    std::cerr << "    " << inputFileDef->fileName << std::endl;
    std::cerr << "    at index " << fileNum << std::endl;
    throw;
  }

  if (inputFileDef->postFitErrors.at(errorType).find(paramName) !=
      inputFileDef->postFitErrors.at(errorType).end())
  {
    return inputFileDef->postFitErrors.at(errorType).at(paramName);
  }
  return __BAD_FLOAT__;
}

const float InputManager::GetPostFitValue(int fileNum, std::string paramName,
                                          std::string errorType) const {
  const InputFile *inputFileDef = fileVec[fileNum];

  // set default type if not specified
  if (errorType == "")
    errorType = inputFileDef->defaultErrorType;

  if (inputFileDef->postFitErrors.find(errorType) == inputFileDef->postFitErrors.end())
  {
    std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << std::endl;
    std::cerr << "  requested error type, "
              << "\"" << errorType << "\" does not exist in the specified file: " << std::endl;
    std::cerr << "    " << inputFileDef->fileName << std::endl;
    std::cerr << "    at index " << fileNum << std::endl;
    throw;
  }

  if (inputFileDef->postFitValues.at(errorType).find(paramName) !=
      inputFileDef->postFitValues.at(errorType).end())
  {
    return inputFileDef->postFitValues.at(errorType).at(paramName);
  }
  return __BAD_FLOAT__;
}

// ##################################################################
// ################## End of public interface #######################
// ##################################################################

std::vector<std::string> InputManager::parseLocation(std::string locationString, fitterEnum fitter,
                                                     fileTypeEnum fileType, std::string parameter,
                                                     std::string sample) const {
  std::vector<std::string> tokens;

  std::size_t pos = 0;
  std::string toReplace = "";

  // loop through the raw location string and replace any instance of {PARAMETER} with the name of
  // the parameter
  toReplace = "{PARAMETER}";
  while ((pos = locationString.find(toReplace)) != std::string::npos)
  {
    locationString.replace(pos, toReplace.length(),
                           getFitterSpecificParamName(fitter, fileType, parameter));
  }

  // same again with {SAMPLE}
  toReplace = "{SAMPLE}";
  while ((pos = locationString.find(toReplace)) != std::string::npos)
  {
    locationString.replace(pos, toReplace.length(),
                           getFitterSpecificSampleName(fitter, fileType, sample));
  }

  // Now we go through and look for ":" and split the location into two parts
  // first part should be the directory to look for the parameter
  // lats part should be the end of the name of the object to look for
  std::string token;
  std::string delimiter = ":";
  while ((pos = locationString.find(delimiter)) != std::string::npos)
  {
    token = locationString.substr(0, pos);
    tokens.push_back(token);
    locationString.erase(0, pos + delimiter.length());
  }
  tokens.push_back(locationString);

  // should only have found 2 parts, anything more means there was too many
  // ":" in the location string
  if (tokens.size() > 2)
  {
    std::cerr << "Too many : tokens in location string" << std::endl;
    throw;
  }

  return tokens;
}

TObject *InputManager::findRootObject(InputFile *fileDef,
                                      std::vector<std::string> locationVec) const {
  TObject *object = nullptr;

  // if vector only has one element, just interpret it as the absolute path to the object
  if (locationVec.size() == 1)
  {
    object = fileDef->file->Get(locationVec[0].c_str());
  }

  // if vector has two elements, interpret the first as the directory and the second as a string to
  // match with the end of the objects in the file
  else if (locationVec.size() == 2)
  {
    TDirectoryFile *directory = (TDirectoryFile *)fileDef->file->Get(locationVec[0].c_str());
    size_t nMatchingObjects = 0;

    // let's make sure that the directory itself exists
    if (directory == NULL)
    {
      object = nullptr;
    }

    else
    {
      // loop through the keys in the directory and find objects whose name matches the specified
      // pattern
      TIter next(directory->GetListOfKeys());
      while (TKey *key = (TKey *)next())
      {
        if (strEndsWith(std::string(key->GetName()), locationVec[1]))
        {
          object = directory->Get(key->GetName());
          nMatchingObjects++;
        }
      }
    }

    // check that only one object matched the pattern
    if (nMatchingObjects > 1)
    {
      std::cerr << "Too many objects match the pattern specified by " << locationVec[0];
      if (locationVec.size() == 2)
        std::cerr << locationVec[1];
      std::cerr << std::endl;
      std::cerr << "Found " << nMatchingObjects << " matching objects, should just be one"
                << std::endl;
      throw;
    }

  }

  // Vector too big!!
  else
  {
    std::cerr << "Invalid object location vector" << std::endl;
    std::cerr
        << "Should have two elements: [ (directory to look in), (end of the name of the object) ]"
        << std::endl;
    std::cerr << "Or one element that is just the absolute path to the object" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  return object;
}


bool InputManager::findBySampleLLH(InputFile *inputFileDef, std::string parameter,
                                   fitterEnum fitter, std::string sample, bool setInputFileScan) {

  YAML::Node thisFitterSpec_config = _fitterSpecConfig[convertFitterNames(fitterEnum(fitter))];

  // EM: Get where the by sample LLH scan for this parameter *should* live if it exists
  YAML::Node testLLHConfig = thisFitterSpec_config["bySample_LLH"];
  std::vector<std::string> testLLHRawLocations =
      testLLHConfig["location"].as<std::vector<std::string>>();

  std::string LLHObjType = thisFitterSpec_config["LLHObjectType"].as<std::string>();

  // EM: Now look for the parameter in this folder
  TObject *LLHObj = nullptr;
  for (std::string rawLocation : testLLHRawLocations)
  {
    LLHObj =
        findRootObject(inputFileDef, parseLocation(rawLocation, fitter, kLLH, parameter, sample));
    if (LLHObj != nullptr)
    {
      break;
    }
  }

  // EM: If it's not in there we can return here
  if (LLHObj == nullptr)
    return false;

  // EM: If specified, we set the object in the InputFile object
  if (setInputFileScan)
  {
    TGraph *LLHGraph = new TGraph();

    if (LLHObjType == "TH1D")
    {
      LLHGraph = new TGraph((TH1D *)LLHObj);
      delete LLHObj;
    } else if (LLHObjType == "TGraph")
    {
      LLHGraph = new TGraph(*(TGraph *)LLHObj);
      delete LLHObj;
    } else
    {
      std::cerr << "ERROR: uknown type of LLH object specified: " << LLHObjType
                << ", can tell me how to read it here" << std::endl;
      std::cerr << "       " << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    inputFileDef->LLHScansBySample_map[sample][parameter] = LLHGraph;
  }

  return true;
}

bool InputManager::findPostFitParamError(InputFile *inputFileDef, std::string parameter,
                                         fitterEnum fitter, std::string errorType,
                                         bool setInputFileError) {
  std::string specificName = getFitterSpecificParamName(fitter, kPostFit, parameter);
  YAML::Node thisFitterSpec_config = _fitterSpecConfig[convertFitterNames(fitterEnum(fitter))];

  // EM: Get which hist this parameter lives in from the config
  YAML::Node postFitErrorTypes = thisFitterSpec_config["postFitErrorTypes"];
  YAML::Node specificErrorType = postFitErrorTypes[errorType];
  std::vector<std::string> postFitLocations =
      specificErrorType["Loc"].as<std::vector<std::string>>();

  // EM: If the parameter has a specified list of locations then override the default one
  std::vector<std::string> postFitLocations_override;
  if (getFitterSpecificParamOption<std::vector<std::string>>(fitter, "postFitLoc",
                                                             &postFitLocations_override, parameter))
  {
    postFitLocations = postFitLocations_override;
  }

  for (std::string postFitLoc : postFitLocations)
  {
    TH1D *postFitErrors = (TH1D *)inputFileDef->file->Get(postFitLoc.c_str());

    // EM: the postfit hist for this parameter isn't in this file
    if (postFitErrors == NULL)
      continue;

    // EM: check the config to see if a suffix is supplied for this
    std::string suffix = thisFitterSpec_config["defaultPostFitSuffix"].as<std::string>();
    getFitterSpecificParamOption<std::string>(fitter, "postFitSuffix", &suffix, parameter);

    // EM: Loop through the hist to see if it contains the parameter we're looking for
    for (int binIdx = 0; binIdx <= postFitErrors->GetNbinsX(); binIdx++)
    {
      std::string binLabel = std::string(postFitErrors->GetXaxis()->GetBinLabel(binIdx));
      if (strEndsWith(binLabel, specificName + suffix))
      {
        if (setInputFileError)
        {
          // EM: if specified, we fill the postfit error TH1D in the provided InputFile object
          inputFileDef->postFitErrors[errorType][parameter] = postFitErrors->GetBinError(binIdx);
          inputFileDef->postFitValues[errorType][parameter] = postFitErrors->GetBinContent(binIdx);
        }

        delete postFitErrors;
        return true;
      }
    }
    delete postFitErrors;
  }
  // EM: Didn't find the parameter in any of the specified locations
  return false;
}

// EM: Lots of room for improvement in fillFileInfo and fillFileData, should be split up into more
// methods, currently a lot of copy pasting
void InputManager::fillFileInfo(InputFile *inputFileDef, bool printThoughts) {
  // use the contents of the file to decide which fitter it came from and what type of file it is
  if (printThoughts)
    std::cout << "Checking contents of file " << inputFileDef->fileName << std::endl;

  for (int i = 1; i < kNFitters - 1; i++)
  {

    fitterEnum fitter = fitterEnum(i);

    // flag for whether or not the current fitterEnum is the correct one
    bool foundFitter = false;
    if (printThoughts)
      std::cout << "Checking if this is a " << convertFitterNames(fitter) << " file" << std::endl;

    // EM: get the configuration specifying what the output of this fitter looks like
    if (!_fitterSpecConfig[convertFitterNames(fitter)])
    {
      std::cerr << "ERROR " << __FILE__ << ":" << __LINE__ << std::endl;
      std::cerr << "translation config doesnt contain a definition for fitter "
                << convertFitterNames(fitter) << std::endl;
      throw;
    }

    YAML::Node thisFitterSpec_config = _fitterSpecConfig[convertFitterNames(fitter)];

    size_t numLLHParams;

    // ##### Look for LLH scans in the input #####
    // check for all 3 LLH directory types
    for (std::string LLHType : {"sample", "penalty", "total"})
    {

      if (printThoughts)
        std::cout << ".... searching for " << LLHType << " LLH scans... ";

      // vector of all the possible locations that we might find LLH scans for this type of LLH
      YAML::Node testLLHConfig = thisFitterSpec_config[LLHType + "_LLH"];
      std::vector<std::string> testLLHRawLocations =
          testLLHConfig["location"].as<std::vector<std::string>>();

      // counter for the total number of parameters we find scans for
      numLLHParams = 0;

      std::vector<std::string> enabledLLHParams;

      for (std::string parameter : knownParameters)
      {

        inputFileDef->availableParams_map_LLH[LLHType][parameter] = false;

        // check if we find the parameter at any of the locations we think it should be at
        for (std::string rawLocation : testLLHRawLocations)
        {
          if (findRootObject(inputFileDef, parseLocation(rawLocation, fitter, kLLH, parameter)) !=
              nullptr)
          {
            numLLHParams++;
            enabledLLHParams.push_back(parameter);
            inputFileDef->availableParams_map_LLH[LLHType][parameter] = true;
            break; // <- we've found it, no point checking the rest of the locations
          }
        }
      }

      if (printThoughts)
        std::cout << "Found " << numLLHParams << std::endl;

      if (numLLHParams > 0)
      {
        foundFitter = true;
        inputFileDef->hasLLHScans = true;
        inputFileDef->availableParams_LLH = enabledLLHParams;
        inputFileDef->hasLLHScans_map[LLHType] = true;
      }

      /// \todo add a check here to make sure all the scans that are in the file are being picked up
      /// by the reader, and warn if any are not being used
    }

    if (printThoughts)
      std::cout << std::endl;

    // ##### now look for processed post fit errors #####
    if (printThoughts)
    {
      std::cout << "....searching for Post Fit Parameters (";

      if (thisFitterSpec_config["defaultPostFitErrorType"])
      {
        std::cout << "Default type specified with possible locations: ";
        YAML::Node postFitErrorSpec = thisFitterSpec_config["postFitErrorTypes"];
        YAML::Node defaultErrorType =
            postFitErrorSpec[thisFitterSpec_config["defaultPostFitErrorType"].as<std::string>()];
        std::vector<std::string> locations = defaultErrorType["Loc"].as<std::vector<std::string>>();
        for (std::string loc : locations)
          std::cout << std::endl << loc;
        std::cout << ")" << std::endl;
      } else
      {
        std::cout << "No default location specified)";
      }
      std::cout << std::endl;
    }

    int numPostFitParams = 0;
    std::vector<std::string> enabledPostFitParams;

    std::string defaultErrorType =
        thisFitterSpec_config["defaultPostFitErrorType"].as<std::string>();
    for (std::string parameter : knownParameters)
    {
      if (findPostFitParamError(inputFileDef, parameter, fitterEnum(i), defaultErrorType))
      {
        numPostFitParams++;
        enabledPostFitParams.push_back(parameter);
      }
    }

    if (printThoughts)
      std::cout << "........Found " << numPostFitParams << std::endl;

    if (numPostFitParams > 0)
    {
      foundFitter = true;
      inputFileDef->hasPostFitErrors = true;
      inputFileDef->availableParams_postFitErrors = enabledPostFitParams;
    }

    // ######### Now look for LLH scans broken down by sample #########
    // This is really just a check to see if the number of by sample LLH scans is the same as normal
    // LLH scans
    if (printThoughts)
      std::cout << "....Searching for LLH scans broken down by sample" << std::endl;

    for (std::string sample : knownSamples)
    {
      size_t numLLHBySampleParams = 0;
      for (std::string parameter : knownParameters)
      {
        inputFileDef->availableParams_map_LLHBySample[sample][parameter] = false;
        if (findBySampleLLH(inputFileDef, parameter, fitterEnum(i), sample))
        {
          inputFileDef->availableParams_map_LLHBySample[sample][parameter] = true;
          numLLHBySampleParams++;
        }
      }

      if (numLLHBySampleParams != 0)
      {
        inputFileDef->availableSamples_LLH.push_back(sample);

        if (printThoughts)
          std::cout << "........Found " << numLLHBySampleParams << " LLH scans for sample "
                    << sample << std::endl;

        if ((numLLHParams != numLLHBySampleParams))
        {
          std::cerr << "ERROR: hmmmmm something weird is happening here" << std::endl;
          std::cerr << "       I have " << numLLHBySampleParams << " LLH scans for sample "
                    << sample << " But " << numLLHParams << " total LLH parameters " << std::endl;
          std::cerr
              << "       If you are at peace with this and understand why this is you can come "
                 "here and remove the throw, otherwise you should investigate whats going on here"
              << std::endl;
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          throw;
        }
      }
    }

    /*
    switch (i) {
        // any other weird fitter specific conditions/ edge cases should go in here
    }
    */

    if (foundFitter)
    {
      if (printThoughts)
        std::cout << "This is a " << convertFitterNames(fitterEnum(i)) << " File!" << std::endl;
      inputFileDef->fitter = fitterEnum(i);
      return;
    }
  }

  // if we didn't return above then the fitter type wasn't found
  if (printThoughts)
    std::cout << "I don't know what kinda fitter this came from, will proceed with caution"
              << std::endl;
  inputFileDef->fitter = fitterEnum::kNFitters;
}

void InputManager::fillFileData(InputFile *inputFileDef, bool printThoughts) {
  // load in the data in the file using the info gotten above
  // EM: a lot of this is copy paste from above, could be better organised
  if (printThoughts)
    std::cout << "....getting data from file " << inputFileDef->fileName << std::endl;

  YAML::Node thisFitterSpec_config = _fitterSpecConfig[convertFitterNames(inputFileDef->fitter)];

  // set the default post fit error type so we can read it as default later
  inputFileDef->defaultErrorType =
      thisFitterSpec_config["defaultPostFitErrorType"].as<std::string>();

  // ########### First fill up the LLH vectors ############
  for (std::string LLHType : {"sample", "penalty", "total"})
  {
    if (!inputFileDef->hasLLHScans_map.at(LLHType))
      continue;

    // vector of all the possible locations that we might find LLH scans for this type of LLH
    YAML::Node testLLHConfig = thisFitterSpec_config[LLHType + "_LLH"];
    std::vector<std::string> testLLHRawLocations =
        testLLHConfig["location"].as<std::vector<std::string>>();

    // get the expected root object type of the llh scans
    std::string LLHObjType = thisFitterSpec_config["LLHObjectType"].as<std::string>();

    // EM: now get the objects from the file
    for (std::string parameter : inputFileDef->availableParams_LLH)
    {

      TObject *LLHObj = nullptr;

      // check the locations for the object
      for (std::string rawLocation : testLLHRawLocations)
      {
        LLHObj = findRootObject(inputFileDef,
                                parseLocation(rawLocation, inputFileDef->fitter, kLLH, parameter));
        // if we found it then break out of the loop
        if (LLHObj != nullptr)
          break;
      }

      // double check that we actually found it.. just in case
      if (LLHObj == nullptr)
      {
        std::cerr << "Hmmm, something seems to have gone wrong and I couldnt find the " << LLHType
                  << " LLH scan for " << parameter << " when attempting to read the data for it"
                  << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      }

      // now convert it to a TGraph
      TGraph *LLHGraph = new TGraph();

      // EM: maybe have the type of the LLH object specified in the config to know what to downcast
      // it to??
      if (LLHObjType == "TH1D")
      {
        LLHGraph = new TGraph((TH1D *)LLHObj);
        delete LLHObj;
      }

      else if (LLHObjType == "TGraph")
      {
        LLHGraph = new TGraph(*(TGraph *)LLHObj);
        delete LLHObj;
      }

      else
      {
        std::cerr << "ERROR: uknown type of LLH object specified: " << LLHObjType
                  << ", can tell me how to read it here" << std::endl;
        std::cerr << "       " << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }

      inputFileDef->LLHScans_map[LLHType][parameter] = LLHGraph;
    }
  }

  // ####### Get the by sample LLH scans #######
  for (std::string parameter : inputFileDef->availableParams_LLH)
  {
    for (std::string sample : knownSamples)
    {
      findBySampleLLH(inputFileDef, parameter, inputFileDef->fitter, sample, true);
    }
  }

  // ####### Get the processed post fit errors #######
  for (std::string parameter : inputFileDef->availableParams_postFitErrors)
  {
    std::vector<std::string> availableErrorTypes =
        thisFitterSpec_config["AvailablePostFitErrorTypes"].as<std::vector<std::string>>();
    for (std::string errorType : availableErrorTypes)
    {
      findPostFitParamError(inputFileDef, parameter, inputFileDef->fitter, errorType, true);
    }
  }
}

// ########################################################################################
// ############# Everything below here is things I would like to get rid of ###############
// ########################################################################################

std::string InputManager::MaCh3ToBANFFName(std::string mach3ParamName) {

  // KS: Need to convert MaCh3 name to BANFF name
  std::string BANFFName = mach3ParamName;

  if (BANFFName.find("_sam") == BANFFName.size() - 4)
  {
    BANFFName.erase(BANFFName.find("_sam"), 4);
    BANFFName += "_scan_sample";
  }

  if (BANFFName.find("_full") == BANFFName.size() - 5)
  {
    BANFFName.erase(BANFFName.find("_full"), 5);
    BANFFName += "_scan_total";
  }

  if (BANFFName.find("_xs") == BANFFName.size() - 3)
  {
    BANFFName.erase(BANFFName.find("_xs"), 3);
    BANFFName += "_scan_penalty";
  }

  // KS: BANFF has special case for Flux bins so we have to hardcode it a bit :(
  if (mach3ParamName.rfind("b_", 0) == 0)
    BANFFName = BANFFfluxName(BANFFName);

  if (!mach3ParamName.find("ndd_"))
    return NOT_FOUND_STR;

  return BANFFName;
}

std::string InputManager::BANFFfluxName(std::string mach3ParamName) {
  // EM: Keep hard coding of BANFF flux parameter names because I cant think of a way of simplifying
  // the name conversion :(

  // KS: Need to convert MaCh3 name to BANFF name
  std::string BANFFName = mach3ParamName;
  std::string Delete = "_sam";

  std::string::size_type i = BANFFName.find(Delete);

  if (i != std::string::npos)
    BANFFName.erase(i, Delete.length());

  BANFFName += "_scan_sample";

  // KS: MaCh3 names are b_0_ up to b_99_
  int FluxNumber = 0;
  char TEMPstring = BANFFName.at(3);
  const char blarb = '_';
  if (TEMPstring == blarb) // if BANFFName.at(3) == '_'
  {
    FluxNumber = int(BANFFName[2] - '0');
    BANFFName.erase(0, 1); // b
    BANFFName.erase(0, 1); // b_
    BANFFName.erase(0, 1); // b_X
  } else
  {

    FluxNumber += 10 * (int)(BANFFName[2] - '0');
    FluxNumber += (int)(BANFFName[3] - '0');
    BANFFName.erase(0, 1); // b
    BANFFName.erase(0, 1); // b_
    BANFFName.erase(0, 1); // b_X
    BANFFName.erase(0, 1); // b_XX
  }
  std::string BANFFprefix;

  std::cout << " DEBUG FluxNumber " << FluxNumber << " BANFFName " << BANFFName << std::endl;

  bool isSKflux = false;
  // if above 50 it is SK flux
  if (FluxNumber >= 50)
  {
    isSKflux = true;
    FluxNumber = FluxNumber - 50;
  }

  if (FluxNumber < 11)
    BANFFprefix = ("NDNuModeNumu" + std::to_string(FluxNumber));

  else if (FluxNumber - 11 < 5)
    BANFFprefix = "NDNuModeNumub" + std::to_string(FluxNumber - 11);

  else if (FluxNumber - 11 - 5 < 7)
    BANFFprefix = "NDNuModeNue" + std::to_string(FluxNumber - 11 - 5);

  else if (FluxNumber - 11 - 5 - 7 < 2)
    BANFFprefix = "NDNuModeNueb" + std::to_string(FluxNumber - 11 - 5 - 7);

  else if (FluxNumber - 11 - 5 - 7 - 2 < 5)
    BANFFprefix = "NDANuModeNumu" + std::to_string(FluxNumber - 11 - 5 - 7 - 2);

  else if (FluxNumber - 11 - 5 - 7 - 2 - 5 < 11)
    BANFFprefix = "NDANuModeNumub" + std::to_string(FluxNumber - 11 - 5 - 7 - 2 - 5);

  else if (FluxNumber - 11 - 5 - 7 - 2 - 5 - 11 < 2)
    BANFFprefix = "NDANuModeNue" + std::to_string(FluxNumber - 11 - 5 - 7 - 2 - 5 - 11);

  else if (FluxNumber - 11 - 5 - 7 - 2 - 5 - 11 - 2 < 7)
    BANFFprefix = "NDANuModeNueb" + std::to_string(FluxNumber - 11 - 5 - 7 - 2 - 5 - 11 - 2);

  // Now replace the "ND" part with "SK" to make the SK names
  if (isSKflux)
    BANFFprefix.replace(BANFFprefix.find("ND"), 2, std::string("SK"));

  std::cout << " DEBUG BANFFName " << BANFFName << " BANFFprefix " << BANFFprefix << std::endl;

  BANFFprefix += BANFFName;

  return BANFFprefix;
}

// EM: annoyingly, need a different function for gundam names for postfit errors than LLH scans as
// file structure is totally different, and so we also need to find the gundam file directory that
// the parameter is in
std::string InputManager::getGundamParamName_PostFit(std::string mach3ParamName,
                                                     TDirectoryFile *fitDir, TH1D *PostFitError_ret,
                                                     Int_t *parIdx) {
  // get the name of the parameter, as it appears in the Hesse and Migrad histogram labels
  std::string gundamName = NOT_FOUND_STR;
  TDirectoryFile *errorDir = (TDirectoryFile *)fitDir->Get("errors");
  TH1D *PostFitErrors;
  int gundamBinIdx;

  if ((mach3ParamName.find("b_") != std::string::npos) && mach3ParamName.size() < 5)
  { // we got a flux parameter
    TDirectoryFile *systDir = (TDirectoryFile *)errorDir->Get("Flux Systematics");
    TDirectoryFile *valuesDir = (TDirectoryFile *)systDir->Get("values");
    PostFitErrors = (TH1D *)valuesDir->Get("postFitErrors_TH1D");

    std::string mach3ParamName_reduced = std::string(mach3ParamName);
    mach3ParamName_reduced.erase(mach3ParamName_reduced.find("b_"), std::string("b_").length());
    gundamName = "#" + mach3ParamName_reduced;

    // EM: get the bin index
    for (gundamBinIdx = 0; gundamBinIdx <= PostFitErrors->GetNbinsX(); gundamBinIdx++)
    {
      std::string binLabel = std::string(PostFitErrors->GetXaxis()->GetBinLabel(gundamBinIdx));
      // EM: we check that mach3ParamName_reduced is at the very end of binLabel, so that we
      // correctly find double digit numbers, e.g. if we a re looking for OOFV_1, it will be found
      // in OOFV_11, OOFV_12 etc.
      if (binLabel == gundamName)
        break;
    }
  }

  else if ((mach3ParamName.find("ND Det ") != std::string::npos))
  { // we got a ND Cov parameter
    TDirectoryFile *systDir =
        (TDirectoryFile *)errorDir->Get("ND280 Detector Systematics (variations parameters)");
    TDirectoryFile *valuesDir = (TDirectoryFile *)systDir->Get("values");
    PostFitErrors = (TH1D *)valuesDir->Get("postFitErrors_TH1D");

    std::string mach3ParamName_reduced = std::string(mach3ParamName);
    mach3ParamName_reduced.erase(mach3ParamName_reduced.find("ND Det "),
                                 std::string("ND Det ").length());
    gundamName = "#" + mach3ParamName_reduced;

    // EM: get the bin index
    for (gundamBinIdx = 0; gundamBinIdx <= PostFitErrors->GetNbinsX(); gundamBinIdx++)
    {
      std::string binLabel = std::string(PostFitErrors->GetXaxis()->GetBinLabel(gundamBinIdx));
      // EM: we check that mach3ParamName_reduced is at the very end of binLabel, so that we
      // correctly find double digit numbers, e.g. if we a re looking for OOFV_1, it will be found
      // in OOFV_11, OOFV_12 etc.
      if (binLabel == gundamName)
        break;
    }
  }

  else
  {
    // EM: Need to check through all the parameter type sub-directories
    //     not sure what the deal is with "Cross-Section FSI Systematics" directory, its
    bool found = false;
    for (std::string paramType : {"Cross-Section (binned) Systematics", "Cross-Section Systematics",
                                  "ND280 Detector Systematics (splined detector parameters)",
                                  "ND280 Detector Systematics (variations parameters)"})
    {
      TDirectoryFile *systDir = (TDirectoryFile *)errorDir->Get(paramType.c_str());
      TDirectoryFile *valuesDir = (TDirectoryFile *)systDir->Get("values");
      PostFitErrors = (TH1D *)valuesDir->Get("postFitErrors_TH1D");

      // EM: strip out the NDS_ at the beginning of the near detector systematic parameteters, this
      // isnt there in gundam
      std::string mach3ParamName_reduced = std::string(mach3ParamName);
      if (mach3ParamName_reduced.find("NDS_") != std::string::npos)
        mach3ParamName_reduced.erase(mach3ParamName_reduced.find("NDS_"),
                                     std::string("NDS_").length());

      else if (mach3ParamName.find("EB_dial_") != std::string::npos)
      {
        mach3ParamName_reduced.erase(mach3ParamName_reduced.find("EB_dial_"),
                                     std::string("EB_dial_").length());
        mach3ParamName_reduced = "EB_bin_" + mach3ParamName_reduced;
      }

      else if (mach3ParamName.find("Alpha") != std::string::npos)
      { mach3ParamName_reduced = "EB_alpha"; }

      for (gundamBinIdx = 0; gundamBinIdx <= PostFitErrors->GetNbinsX(); gundamBinIdx++)
      {
        std::string binLabel = std::string(PostFitErrors->GetXaxis()->GetBinLabel(gundamBinIdx));
        // EM: we check that mach3ParamName_reduced is at the very end of binLabel, so that we
        // correctly find double digit numbers, e.g. if we a re looking for OOFV_1, it will be found
        // in OOFV_11, OOFV_12 etc.
        if (binLabel.find(mach3ParamName_reduced) != std::string::npos &&
            binLabel.find(mach3ParamName_reduced) ==
                binLabel.size() - mach3ParamName_reduced.size())
        {
          gundamName = binLabel;

          found = true;
          break;
        }
      }

      if (found)
        break;
    }
  }
  // std::string gundamName_reduced = getGundamParamName_LLHScan(mach3ParamName, gundamFile);

  // EM: fill the optional return values, the error TH1D that this parameter was found in and the
  // index where it was found
  if (PostFitError_ret != NULL)
    *PostFitError_ret = *PostFitErrors;
  if (parIdx != NULL)
    *parIdx = gundamBinIdx;

  delete PostFitErrors;

  return gundamName;
}
} // namespace MaCh3Plotting