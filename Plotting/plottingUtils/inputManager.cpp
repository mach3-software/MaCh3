#include "inputManager.h"


namespace MaCh3Plotting {
// this is the constructor with user specified translation config file
InputManager::InputManager(const std::string &translationConfigName) {
  // read the config file
  _translatorConfig = M3OpenConfig(translationConfigName);

  MACH3LOG_DEBUG("InputManager: have loaded translation config file");
  
  // split up the parts of the config for easier access later
  _fitterSpecConfig = _translatorConfig["FitterSpec"];
  _parametersConfig = _translatorConfig["Parameters"];
  _samplesConfig = _translatorConfig["Samples"];

  // check the config file and get which parameters, samples, and fitters we've been told about
  _knownFitters = Get<std::vector<std::string>>(_fitterSpecConfig["fitters"], __FILE__, __LINE__);
  _knownParameters = Get<std::vector<std::string>>(_parametersConfig["Parameters"], __FILE__, __LINE__);
  _knownSamples = Get<std::vector<std::string>>(_samplesConfig["Samples"], __FILE__, __LINE__);

  MACH3LOG_DEBUG("Will now check the specified parameters for tags");

  // loop through all parameters and get their tags
  for ( const std::string &param: _knownParameters )
  {
    std::vector<std::string> tags;

    MACH3LOG_DEBUG("Looking for tags for parameter {}", param);
    if ( _parametersConfig[param] )
    {
      if ( _parametersConfig[param]["tags"] )
      {
        tags = Get<std::vector<std::string>>(_parametersConfig[param]["tags"], __FILE__, __LINE__);
        MACH3LOG_DEBUG("  - Found {}!", tags.size());
      }
    }
    _paramToTagsMap[param] = tags;
  }

  MACH3LOG_DEBUG("Will now check the specified samples for tags");
  
  // same again for samples
  for ( const std::string &samp: _knownSamples )
  {
    std::vector<std::string> tags;

    MACH3LOG_DEBUG("Looking for tags for sample {}", samp);
    if ( _samplesConfig[samp])
    {
      if ( _samplesConfig[samp]["tags"] )
      {
        tags = Get<std::vector<std::string>>(_samplesConfig[samp]["tags"], __FILE__, __LINE__);
        MACH3LOG_DEBUG("  - Found {}!", tags.size());
      }
    }
    _sampleToTagsMap[samp] = tags;
  }
}

/// Open an input file and add to the manager, consists of:
///  -  a new InputFile using the specified file
///  - Get info about the file, like what fitter it came from, what it contains e.g. LLH scans,
///  processed post fit parameters etc.
///  - Load up the data from the file (LLH scans etc.) and put them in a common format to be used by
///  plotting scripts
///  - Push back a pointer to the InputFile objcet to the vector of files known to this
///  InputManager.
void InputManager::addFile(const std::string &fileName) {
  _fileVec.emplace_back(fileName);

  // EM: need to be done in this order since fillFileData needs to know info about the file, e.g.
  // fitter and what things are in it
  InputFile &fileInfo = _fileVec.back();
  fillFileInfo(fileInfo);
  fillFileData(fileInfo);
}

/// If printLevel is "summary", will loop through all the files known to this InputManager and call
/// InputFile::Summarise(). If printLevel is "dump", will print all parameters known to this
/// InputManager and then loop through all input files and call InputFile::Dump().
void InputManager::print(const std::string &printLevel) const {
  MACH3LOG_INFO("Printing contents of InputManager instance:");

  if (printLevel == "dump")
  {
    MACH3LOG_INFO("parameters known to this manager: ");
    for (std::string param : _knownParameters)
    {
      MACH3LOG_INFO("  ");
    }
  }

  int fileCount = 0;
  for (const InputFile &file : _fileVec)
  {
    MACH3LOG_INFO(" For file {}", fileCount);
    if (printLevel == "summary")
    {
      file.Summarise();
    } else if (printLevel == "dump")
    { file.Dump(); }
    fileCount++;
  }

  MACH3LOG_INFO("");
}

double InputManager::getPostFitError(const int fileNum, const std::string &paramName,
                                     std::string errorType) const {
  const InputFile &inputFileDef = getFile(fileNum);

  // set default type if not specified
  if (errorType == "")
    errorType = inputFileDef.defaultErrorType;

  if (inputFileDef.postFitErrors.find(errorType) == inputFileDef.postFitErrors.end())
  {
    MACH3LOG_CRITICAL("Requested error type, {} does not exist in the specified file: ", errorType);
    MACH3LOG_CRITICAL("  {}", inputFileDef.fileName);
    MACH3LOG_CRITICAL("  at index {}", fileNum);

    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (inputFileDef.postFitErrors.at(errorType).find(paramName) !=
      inputFileDef.postFitErrors.at(errorType).end())
  {
    return inputFileDef.postFitErrors.at(errorType).at(paramName);
  }

  MACH3LOG_WARN("Didn't find {} post fit error for {}. Returning {}", errorType, paramName, M3::_BAD_DOUBLE_);

  return M3::_BAD_DOUBLE_;
}

double InputManager::getPostFitValue(const int fileNum, const std::string &paramName,
                                     std::string errorType) const {
  const InputFile &inputFileDef = getFile(fileNum);

  // set default type if not specified
  if (errorType == "")
    errorType = inputFileDef.defaultErrorType;

  if (inputFileDef.postFitErrors.find(errorType) == inputFileDef.postFitErrors.end())
  {
    MACH3LOG_CRITICAL("Requested error type, {} does not exist in the specified file: ", errorType);
    MACH3LOG_CRITICAL("  {}", inputFileDef.fileName);
    MACH3LOG_CRITICAL("  at index {}", fileNum);

    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (inputFileDef.postFitValues.at(errorType).find(paramName) !=
      inputFileDef.postFitValues.at(errorType).end())
  {
    return inputFileDef.postFitValues.at(errorType).at(paramName);
  }
  
  MACH3LOG_WARN("Didn't find {} post fit value for {}. Returning {}", errorType, paramName, M3::_BAD_DOUBLE_);

  return M3::_BAD_DOUBLE_;
}


// ##################################################################
// ################## End of public interface #######################
// ##################################################################

std::vector<std::string> InputManager::getTaggedValues(const std::vector<std::string> &values, 
                                                       const std::unordered_map<std::string, std::vector<std::string>> &tagMap, 
                                                       const std::vector<std::string> &tags, std::string checkType) const {
  // check that checkType is valid
  if(
    checkType != "all" &&
    checkType != "any" &&
    checkType != "exact"
  )
  {  
    MACH3LOG_ERROR("Invalid tag check type specified: {}. Will instead use the default type: 'all'", checkType);
    checkType = "all";
  }

  // If no tags were specified, take this to mean that anything should be a match
  if (tags.size() == 0) return values;

  std::vector<std::string> retVec;

  MACH3LOG_DEBUG("Getting tagged values using checkType {}", checkType);
  for ( std::string val: values )
  {
    MACH3LOG_DEBUG("Checking tags of {}", val);
    // get the tags that this value has
    const std::vector<std::string> &valTags = tagMap.at(val);

    // we can skip it if it has no tags
    if (valTags.size() == 0) continue;

    // count how many of the specified tags match the tags of the current value 
    unsigned int tagCount = 0;
    for ( const std::string &tag: tags )
    {
      if ( std::find( valTags.begin(), valTags.end(), tag ) != valTags.end() ) 
      {    
        MACH3LOG_DEBUG("  - Matched tag {} !", tag);
        tagCount ++;
      }
    }

    // now decide if we include the current value based on the check type
    if ( checkType == "all" )
    {
      if ( tagCount == tags.size() ) retVec.push_back(val);
    }
    else if ( checkType == "any" )
    {
      if ( tagCount > 0 ) retVec.push_back(val);
    }
    else if ( checkType == "exact" )
    {
      // EM: note that this will break if duplicate tags are specified in either vector... so please don't do that
      if ( tagCount == valTags.size() ) retVec.push_back(val);
    }
  }
  MACH3LOG_DEBUG("Found {} values matching the specified tags", retVec.size());
  return retVec;
}

std::vector<std::string> InputManager::parseLocation(const std::string &rawLocationString, std::string &fitter,
                                                     fileTypeEnum fileType, const std::string &parameter,
                                                     const std::string &sample, const std::string &parameter2) const {
  std::string locationString(rawLocationString);
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

  // loop through the raw location string and replace any instance of {PARAMETER2} with the name of
  // the second specified parameter
  toReplace = "{PARAMETER2}";
  while ((pos = locationString.find(toReplace)) != std::string::npos)
  {
    locationString.replace(pos, toReplace.length(),
                           getFitterSpecificParamName(fitter, fileType, parameter2));
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
    throw MaCh3Exception(__FILE__ , __LINE__, "Too many : tokens in location string: " + rawLocationString);
  }

  return tokens;
}

std::shared_ptr<TObject> InputManager::findRootObject(const InputFile &fileDef,
                                      const std::vector<std::string> &locationVec) const {
  std::shared_ptr<TObject> object = nullptr;

  // if vector only has one element, just interpret it as the absolute path to the object
  if (locationVec.size() == 1)
  {
    object = std::shared_ptr<TObject>(fileDef.file->Get(locationVec[0].c_str()));
  }

  // if vector has two elements, interpret the first as the directory and the second as a string to
  // match with the end of the objects in the file
  else if (locationVec.size() == 2)
  {
    TDirectoryFile *directory = fileDef.file->Get<TDirectoryFile>(locationVec[0].c_str());
    size_t nMatchingObjects = 0;

    // let's make sure that the directory itself exists
    if (directory == nullptr)
    {
      object = nullptr;
    }
    else
    {
      // loop through the keys in the directory and find objects whose name matches the specified
      // pattern
      TIter next(directory->GetListOfKeys());
      while (TKey *key = static_cast<TKey*>(next()))
      {
        if (strEndsWith(std::string(key->GetName()), locationVec[1]))
        {
          object = std::shared_ptr<TObject>(directory->Get(key->GetName()));
          nMatchingObjects++;
        }
      }
    }
    // check that only one object matched the pattern
    if (nMatchingObjects > 1)
    {
      MACH3LOG_CRITICAL("Too many objects match the pattern specified by {} {}", locationVec[0], locationVec.size()==2 ? locationVec[1] : "");
      MACH3LOG_CRITICAL("Found {} matching objects, should just be one", nMatchingObjects);

      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }
  else // Vector too big!!
  {
    MACH3LOG_CRITICAL("Invalid object location vector");
    MACH3LOG_CRITICAL("Should have two elements: [ (directory to look in), (end of the name of the object) ]");
    MACH3LOG_CRITICAL("Or one element that is just the absolute path to the object");
    
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return object;
}

bool InputManager::findBySampleLLH(InputFile &inputFileDef, const std::string &parameter,
                                   std::string &fitter, const std::string &sample, bool setInputFileScan) {
  YAML::Node thisFitterSpec_config = _fitterSpecConfig[fitter];

  // EM: Get where the by sample LLH scan for this parameter *should* live if it exists
  YAML::Node testLLHConfig = thisFitterSpec_config["bySample_LLH"];
  auto testLLHRawLocations = Get<std::vector<std::string>>(testLLHConfig["location"], __FILE__, __LINE__);

  auto LLHObjType = Get<std::string>(thisFitterSpec_config["LLHObjectType"], __FILE__, __LINE__);
  // EM: Now look for the parameter in this folder
  std::shared_ptr<TObject> LLHObj = nullptr;
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
    std::shared_ptr<TGraph> LLHGraph = std::make_shared<TGraph>();

    if (LLHObjType == "TH1D")
    {
      LLHGraph = std::make_shared<TGraph>(static_cast<TH1D*>(LLHObj.get()));
    } else if (LLHObjType == "TGraph")
    {
      LLHGraph = std::shared_ptr<TGraph>(static_cast<TGraph*>(LLHObj->Clone()));
    } else
    {
      throw MaCh3Exception(__FILE__ , __LINE__, "uknown type of LLH object specified: " + LLHObjType);
    }

    inputFileDef.LLHScansBySample_map[sample][parameter] = LLHGraph;
  }
  return true;
}

// check the input file for raw MCMC step values for a particular parameter
bool InputManager::findRawChainSteps(InputFile &inputFileDef, const std::string &parameter, std::string &fitter, bool setInputBranch) const {
  // we'll assume for now that the chain is in the form of a TTree and the branch names are parameter names
  bool wasFound = false;

  // make sure that the filedef object has all the necessary stuff to read from the posterior tree
  if ( (inputFileDef.mcmcProc != nullptr) && (inputFileDef.posteriorTree != nullptr) )
  {
    const std::vector<TString> branchNames = inputFileDef.mcmcProc->GetBranchNames();

    std::string specificName = getFitterSpecificParamName(fitter, kMCMC, parameter);

    // loop over possible parameters and compare names
    for ( int paramIdx = 0; paramIdx < inputFileDef.mcmcProc->GetNParams() -1 ; paramIdx ++ )
    {
      TString title;
      double prior, priorError; // <- will be discarded
      inputFileDef.mcmcProc->GetNthParameter(paramIdx, prior, priorError, title);

      if ( strEndsWith(title.Data(), specificName) )
      {
        wasFound = true;
        if ( setInputBranch )
        {
          // EM: should probably use MCMCProcessor for this so we can use caching, gpu etc.
          inputFileDef.MCMCstepParamsMap[parameter] = new double( M3::_BAD_DOUBLE_ ); // <- initialise the parameter step values
          inputFileDef.posteriorTree->SetBranchAddress( branchNames[paramIdx], inputFileDef.MCMCstepParamsMap.at(parameter) );
        }
        break;
      }
    }
  }
  return wasFound;
}

// check the input file for processed 1d posteriors for a particular parameter
bool InputManager::find1dPosterior(InputFile &inputFileDef, const std::string &parameter, std::string &fitter, bool setFileData) const {
  bool wasFound = false;
  
  YAML::Node thisFitterSpec_config = _fitterSpecConfig[fitter];

  if ( thisFitterSpec_config["1dPosteriors"] ) {
    auto rawLocations = Get<std::vector<std::string>>(thisFitterSpec_config["1dPosteriors"]["location"],
                                                      __FILE__, __LINE__);
    for ( const std::string &rawLoc : rawLocations)
    {
      std::shared_ptr<TH1D> posterior1d = std::static_pointer_cast<TH1D>(findRootObject(inputFileDef, parseLocation(rawLoc, fitter, kMCMC, parameter)));

      if ( posterior1d != nullptr )
      {
        wasFound = true;

        if ( setFileData )
        {
          inputFileDef.posteriors1d_map[parameter] = std::make_shared<TGraph>(posterior1d.get());
        }
        break;
      }
    }
  }
  return wasFound;
}

bool InputManager::findPostFitParamError(InputFile &inputFileDef, const std::string &parameter,
                                         std::string &fitter, const std::string &errorType,
                                         const bool setInputFileError) {
  std::string specificName = getFitterSpecificParamName(fitter, kPostFit, parameter);
  YAML::Node thisFitterSpec_config = _fitterSpecConfig[fitter];

  // EM: Get which hist this parameter lives in from the config
  YAML::Node postFitErrorTypes = thisFitterSpec_config["postFitErrorTypes"];
  YAML::Node specificErrorType = postFitErrorTypes[errorType];
  auto postFitLocations = Get<std::vector<std::string>>(specificErrorType["location"], __FILE__, __LINE__);

  // EM: If the parameter has a specified list of locations then override the default one
  std::vector<std::string> postFitLocations_override;
  if (getFitterSpecificParamOption<std::vector<std::string>>(fitter, "postFitLoc",
                                                             postFitLocations_override, parameter))
  {
    postFitLocations = postFitLocations_override;
  }

  for (std::string postFitLoc : postFitLocations)
  {
    std::shared_ptr<TH1D> postFitErrors = std::static_pointer_cast<TH1D>(findRootObject(inputFileDef, parseLocation(postFitLoc, fitter, kPostFit, parameter)));

    // EM: the postfit hist for this parameter isn't in this file
    if (postFitErrors == nullptr)
      continue;

    // EM: Loop through the hist to see if it contains the parameter we're looking for
    for (int binIdx = 0; binIdx <= postFitErrors->GetNbinsX(); binIdx++)
    {
      std::string binLabel = std::string(postFitErrors->GetXaxis()->GetBinLabel(binIdx));
      if (strEndsWith(binLabel, specificName))
      {
        if (setInputFileError)
        {
          // EM: if specified, we fill the postfit error TH1D in the provided InputFile object
          inputFileDef.postFitErrors[errorType][parameter] = postFitErrors->GetBinError(binIdx);
          inputFileDef.postFitValues[errorType][parameter] = postFitErrors->GetBinContent(binIdx);
        }

        return true;
      }
    }
  }
  // EM: Didn't find the parameter in any of the specified locations
  return false;
}

// EM: Lots of room for improvement in fillFileInfo and fillFileData, should be split up into more
// methods, currently a lot of copy pasting
void InputManager::fillFileInfo(InputFile &inputFileDef, const bool printThoughts) {
  /// @todo Would like to be able to specify what kind of file and what fitter an input is from on the command like:
  /// e.g. like `plotApp [options] <fileName1>;<fileType>;<fitterName>... and only try to auto-detect it if its not specified, 
  /// this would save some time and would also be very helpful in situations where we can't auto-detect e.g. if there is some kind of 
  ///overlap in the file structure between two fitters
  
  // use the contents of the file to decide which fitter it came from and what type of file it is
  if (printThoughts)
    MACH3LOG_INFO("Checking contents of file {}", inputFileDef.fileName);

  for (std::string fitter: _knownFitters)
  {
    // flag for whether or not the current fitter is the correct one
    bool foundFitter = false;
    if (printThoughts)
      MACH3LOG_INFO("Checking if this is a {} file", fitter);

    // EM: get the configuration specifying what the output of this fitter looks like
    if (!_fitterSpecConfig[fitter])
    {
      throw MaCh3Exception(__FILE__ , __LINE__, "translation config doesnt contain a definition for fitter " + fitter);
    }

    YAML::Node thisFitterSpec_config = _fitterSpecConfig[fitter];

    size_t numLLHParams;

    // ##### Look for LLH scans in the input #####
    // check for all 3 LLH directory types
    for (std::string LLHType : {"sample", "penalty", "total"})
    {
      if (printThoughts)
        MACH3LOG_INFO(".... searching for {} LLH scans... ", LLHType);

      // vector of all the possible locations that we might find LLH scans for this type of LLH
      YAML::Node testLLHConfig = thisFitterSpec_config[LLHType + "_LLH"];
      auto testLLHRawLocations = Get<std::vector<std::string>>(testLLHConfig["location"], __FILE__, __LINE__);

      // counter for the total number of parameters we find scans for
      numLLHParams = 0;

      std::vector<std::string> enabledLLHParams;

      for (const std::string &parameter : _knownParameters)
      {
        MACH3LOG_DEBUG("     - for {}", parameter);
        inputFileDef.availableParams_map_LLH[LLHType][parameter] = false;

        // check if we find the parameter at any of the locations we think it should be at
        for (const std::string &rawLocation : testLLHRawLocations)
        {
          if (findRootObject(inputFileDef, parseLocation(rawLocation, fitter, kLLH, parameter)) !=
              nullptr)
          {
            numLLHParams++;
            enabledLLHParams.push_back(parameter);
            inputFileDef.availableParams_map_LLH[LLHType][parameter] = true;
            MACH3LOG_DEBUG("      FOUND!");
            break; // <- we've found it, no point checking the rest of the locations
          }
        }
      }

      if (printThoughts)
        MACH3LOG_INFO(".... Found {}", numLLHParams);

      if (numLLHParams > 0)
      {
        foundFitter = true;
        inputFileDef.hasLLHScans = true;
        inputFileDef.availableParams_LLH = enabledLLHParams;
        inputFileDef.hasLLHScans_map[LLHType] = true;
      }

      /// \todo add a check here to make sure all the scans that are in the file are being picked up
      /// by the reader, and warn if any are not being used
    }

    if (printThoughts)
      MACH3LOG_INFO("");

    // ##### now look for processed post fit errors #####
    if (printThoughts)
    {
      MACH3LOG_INFO("....searching for Post Fit Parameters");

      if (thisFitterSpec_config["defaultPostFitErrorType"])
      {
        MACH3LOG_DEBUG("  Default type specified with possible locations: ");
        YAML::Node postFitErrorSpec = thisFitterSpec_config["postFitErrorTypes"];
        YAML::Node defaultErrorType =
            postFitErrorSpec[thisFitterSpec_config["defaultPostFitErrorType"].as<std::string>()];
        auto locations = Get<std::vector<std::string>>(defaultErrorType["location"], __FILE__, __LINE__);
        for (std::string loc : locations)
          MACH3LOG_DEBUG(loc);
      } 
      else
      {
        MACH3LOG_DEBUG("  No default location specified");
      }
    }

    int numPostFitParams = 0;
    std::vector<std::string> enabledPostFitParams;

    auto defaultErrorType =
        Get<std::string>(thisFitterSpec_config["defaultPostFitErrorType"], __FILE__, __LINE__);
    for (std::string parameter : _knownParameters)
    {
      if (findPostFitParamError(inputFileDef, parameter, fitter, defaultErrorType))
      {
        numPostFitParams++;
        enabledPostFitParams.push_back(parameter);
      }
    }

    if (printThoughts)
      MACH3LOG_INFO(".... Found {}", numPostFitParams);

    if (numPostFitParams > 0)
    {
      foundFitter = true;
      inputFileDef.hasPostFitErrors = true;
      inputFileDef.availableParams_postFitErrors = enabledPostFitParams;
    }

    // ######### Now look for LLH scans broken down by sample #########
    // This is really just a check to see if the number of by sample LLH scans is the same as normal
    // LLH scans
    if (printThoughts)
      MACH3LOG_INFO("....Searching for LLH scans broken down by sample");

    for (const std::string &sample : _knownSamples)
    {
      size_t numLLHBySampleParams = 0;
      for (const std::string &parameter : _knownParameters)
      {
        inputFileDef.availableParams_map_LLHBySample[sample][parameter] = false;
        if (findBySampleLLH(inputFileDef, parameter, fitter, sample))
        {
          inputFileDef.availableParams_map_LLHBySample[sample][parameter] = true;
          numLLHBySampleParams++;
        }
      }

      if (numLLHBySampleParams != 0)
      {
        inputFileDef.availableSamples_LLH.push_back(sample);

        if (printThoughts)
          MACH3LOG_INFO("........ Found {} LLH scans for sample {}", numLLHBySampleParams, sample);

        if ((numLLHParams != numLLHBySampleParams))
        {
          MACH3LOG_ERROR("hmmmmm something weird is happening here");
          MACH3LOG_ERROR("  I have {} LLH scans for sample {}", numLLHBySampleParams, sample);
          MACH3LOG_ERROR("  But {} parameters with Total_LLH scals", numLLHParams);
        }
      }
    }

    // ######### Now for the main event: Look for MCMC chains and processed posteriors ###########
    
    // EM: if "MCMCsteps" was defined for this fitter, we assume that it is a MaCh3 raw MCMC file
    // thus it needs an MCMCProcessor to read from it. This isn't super general and it would probably
    // be good to have some additional "isMaCh3" option that can decide whether or not to use MCMCProcessor
    // but hey ho it's good enough for now
    if ( thisFitterSpec_config["MCMCsteps"] )
    {
      MACH3LOG_DEBUG("Initialising MCMCProcessor for the input file");
      auto posteriorTreeRawLocations =
            Get<std::vector<std::string>>(thisFitterSpec_config["MCMCsteps"]["location"], __FILE__, __LINE__);

      TTree *postTree = nullptr;
      for ( const std::string &rawLoc: posteriorTreeRawLocations )
      {
        MACH3LOG_DEBUG("  - Looking for MCMC chain parameter values at: {}", rawLoc);

        postTree = inputFileDef.file->Get<TTree>(rawLoc.c_str());

        if ( postTree != nullptr )
        {
          inputFileDef.mcmcProc = new MCMCProcessor(inputFileDef.fileName);
          inputFileDef.mcmcProc->Initialise();

          MACH3LOG_DEBUG("  - FOUND!");
          break;
        }
      }
    
      if ( postTree != nullptr )
      {
        inputFileDef.posteriorTree = postTree;
        inputFileDef.nMCMCentries = int(postTree->GetEntries());
      }
    }

    if (printThoughts)
    {
      MACH3LOG_INFO("....Searching for MCMC related things");
    }
    
    size_t num1dPosteriors = 0;
    std::vector<std::string> enabled1dPosteriorParams;

    size_t numMCMCchainParams = 0;
    std::vector<std::string> enabledMCMCchainParams;

    for ( const std::string &parameter : _knownParameters )
    {
      MACH3LOG_DEBUG("     - for {}", parameter);
      // check for 1d post processing posterior.q
      inputFileDef.availableParams_map_1dPosteriors[parameter] = false;
      if ( thisFitterSpec_config["1dPosteriors"] && find1dPosterior(inputFileDef, parameter, fitter) )
      {
        MACH3LOG_DEBUG("       Found 1d processed posterior!");
        enabled1dPosteriorParams.push_back(parameter);
        num1dPosteriors++;
        inputFileDef.availableParams_map_1dPosteriors[parameter] = true;
      }
      // now check for parameters in chain 
      inputFileDef.availableParams_map_MCMCchain[parameter] = false;
      if ( thisFitterSpec_config["MCMCsteps"] && findRawChainSteps(inputFileDef, parameter, fitter) )
      {
        MACH3LOG_DEBUG("       Found raw MCMC steps!");
        enabledMCMCchainParams.push_back(parameter);
        numMCMCchainParams++;
        inputFileDef.availableParams_map_MCMCchain[parameter] = true;
      }
    }

    if  (num1dPosteriors > 0 )
    {
      foundFitter = true;
      inputFileDef.has1dPosteriors = true;
      inputFileDef.availableParams_1dPosteriors = enabled1dPosteriorParams;

      if (printThoughts)
        MACH3LOG_INFO("........ Found {} 1d processed posteriors", num1dPosteriors);
    }
    if ( numMCMCchainParams > 0 )
    {
      foundFitter = true;
      inputFileDef.hasMCMCchain = true;
      inputFileDef.availableParams_MCMCchain = enabledMCMCchainParams;

      if (printThoughts)
        MACH3LOG_INFO("........ Found {} parameters in MCMC chain", numMCMCchainParams);
    }

    /*
    switch (i) {
        // any other weird fitter specific conditions/ edge cases should go in here
    }
    */

    if (foundFitter)
    {
      if (printThoughts)
      {
        MACH3LOG_INFO("This is a {} file!", fitter);
      }
      inputFileDef.fitter = fitter;
      return;
    }
  }

  // if we didn't return above then the fitter type wasn't found
  if (printThoughts)
    MACH3LOG_WARN("I don't know what kinda fitter this came from, will proceed with caution");
  inputFileDef.fitter = "UNKNOWN_FITTER";
}

void InputManager::fillFileData(InputFile &inputFileDef, const bool printThoughts) {
  // load in the data in the file using the info gotten above
  // EM: a lot of this is copy paste from above, could be better organised
  if (printThoughts)
    MACH3LOG_INFO("....getting data from file {}", inputFileDef.fileName);

  YAML::Node thisFitterSpec_config = _fitterSpecConfig[inputFileDef.fitter];

  // set the default post fit error type so we can read it as default later
  inputFileDef.defaultErrorType =
      Get<std::string>(thisFitterSpec_config["defaultPostFitErrorType"], __FILE__, __LINE__);

  // ########### First fill up the LLH vectors ############
  for (const std::string LLHType : {"sample", "penalty", "total"})
  {
    if (!inputFileDef.hasLLHScans_map.at(LLHType))
      continue;

    // vector of all the possible locations that we might find LLH scans for this type of LLH
    YAML::Node testLLHConfig = thisFitterSpec_config[LLHType + "_LLH"];
    const std::vector<std::string> testLLHRawLocations =
        Get<std::vector<std::string>>(testLLHConfig["location"], __FILE__, __LINE__);

    // get the expected root object type of the llh scans
    auto LLHObjType = Get<std::string>(thisFitterSpec_config["LLHObjectType"], __FILE__, __LINE__);
    // EM: now get the objects from the file
    for (const std::string &parameter : inputFileDef.availableParams_LLH)
    {
      std::shared_ptr<TObject> LLHObj = nullptr;

      // check the locations for the object
      for (const std::string &rawLocation : testLLHRawLocations)
      {
        LLHObj = findRootObject(inputFileDef,
                                parseLocation(rawLocation, inputFileDef.fitter, kLLH, parameter));
        // if we found it then break out of the loop
        if (LLHObj != nullptr)
          break;
      }

      // double check that we actually found it.. just in case
      if (LLHObj == nullptr)
      {
        MACH3LOG_ERROR("Hmmm, something seems to have gone wrong and I couldnt find the {} LLH scan for {} when attempting to read the data for it", LLHType, parameter);
        MACH3LOG_ERROR("This will very likely cause segfaults and sadness");
      }

      // now convert it to a TGraph
      std::shared_ptr<TGraph> LLHGraph = std::make_shared<TGraph>();

      // EM: maybe have the type of the LLH object specified in the config to know what to downcast
      // it to??
      if (LLHObjType == "TH1D")
      {
        LLHGraph = std::make_shared<TGraph>(static_cast<TH1D*>(LLHObj.get()));
      }

      else if (LLHObjType == "TGraph")
      {
        LLHGraph = std::shared_ptr<TGraph>(static_cast<TGraph*>(LLHObj->Clone()));
      }

      else
      {
        MACH3LOG_CRITICAL("ERROR: uknown type of LLH object specified: {}", LLHObjType);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }

      inputFileDef.LLHScans_map[LLHType][parameter] = LLHGraph;
    }
  }

  // ####### Get the by sample LLH scans #######
  for (const std::string &parameter : inputFileDef.availableParams_LLH)
  {
    for (const std::string &sample : _knownSamples)
    {
      findBySampleLLH(inputFileDef, parameter, inputFileDef.fitter, sample, true);
    }
  }

  // ####### Get the processed post fit errors #######
  for (const std::string &parameter : inputFileDef.availableParams_postFitErrors)
  {
    auto availableErrorTypes = Get<std::vector<std::string>>(thisFitterSpec_config["AvailablePostFitErrorTypes"],
                                                             __FILE__, __LINE__);

    for (const std::string &errorType : availableErrorTypes)
    {
      findPostFitParamError(inputFileDef, parameter, inputFileDef.fitter, errorType, true);
    }
  }

  // ########## Get the MCMC related posteriors ###########
  for (const std::string &parameter : inputFileDef.availableParams_MCMCchain)
  {
    findRawChainSteps(inputFileDef, parameter, inputFileDef.fitter, true);
  }
  for (const std::string &parameter : inputFileDef.availableParams_1dPosteriors)
  {
    find1dPosterior(inputFileDef, parameter, inputFileDef.fitter, true);
  }
}
} // namespace MaCh3Plotting
