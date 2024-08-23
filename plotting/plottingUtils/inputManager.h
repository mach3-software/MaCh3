#pragma once

// MaCh3 Includes
#include "manager/MaCh3Logger.h"
#include "manager/YamlHelper.h"
#include "manager/MaCh3Exception.h"

// Other plotting includes
#include "plottingUtils.h"

namespace MaCh3Plotting {
/// @brief Possible fitters that a file can come from.
/// @todo I don't think this is actually necessary, could make it so list of possible fitters is
/// taken straight from which ones are defined in the translation config file.
enum fitterEnum {
  kMaCh3,   //!< Came from MaCh3, not used but might be useful in cases where MaCh3_ND and MaCh3_FD
            //!< are indistinguishable?
  kMaCh3ND, //!< Came from MaCh3 near detector executable
  kMaCh3FD, //!< Came from MaCh3 far detector executable
  kGundam,  //!< Came from GUNDAM
  kBANFF,   //!< Came from BANFF (Not fully implemented)
  // kValor, //!< Came from Valor (Not at all implemented)
  // kPTheta, //!< Came from PTheta (Not at all implemented)

  kNFitters //!< Number of fitters known
};

/// @brief Types of possible file that can be read.
enum fileTypeEnum {
  kLLH,         //!< Log Likelihood scan
  kPostFit,     //!< Processed post fit errors
  kMarkovChain, //!< MCMC chain
  kSigmaVar,    //!< Sigma variations

  kNFileTypes //!< Number of types of file
};

struct InputFile {
  /// @brief Struct which wraps around the actual input file and also holds general information, and
  /// data from the file to be used by other classes.
  /// @details This is just a little guy thats
  /// intended to be just a plain old data type which simply holds and organises information and
  /// data read from the underlying input file, just some methods to inspect the contents of the
  /// object. As such it does not implement methods to manipulate such data, which is all done via
  /// the InputManager() class. Reason for this object is to try to smooth over the cracks that come
  /// from very different file structures used by different fitters and provide a common format to
  /// be used for plot making. Intended use of this objects is only via InputManager which contains
  /// the methods to fill and manipulate it.

  /// @brief Create InputFile instance based on a specified file.
  /// @param fName The name of the file to open.
  /// @return The constructed InputFile.
  /// @details Constructor which opens the specified file and sets default values for information
  /// parameters. Currently this assumes the specified file is a root file, and opens it as a TFile.
  /// In future it would be nice to extend this to be able to read different types of files, perhaps
  /// also moving the opening of the file to a different place to keep this struct nice and simple,
  /// perhaps into the InputManager.
  InputFile(const std::string &fName) {
    fileName = fName;

    file = std::make_shared<TFile>(fileName.c_str());
    hasLLHScans = false;
    hasPostFitErrors = false;
    hasSigmaVars = false;

    for (std::string LLHType : {"sample", "penalty", "total"})
    {
      hasLLHScans_map[LLHType] = false;
    }
  }

  /// @brief Destructor.
  ~InputFile() {
    LLHScans_map.clear();
    file->Close();
  }

  /// @brief Print out a small summary of what is contained in the file.
  /// @todo Could add some flag to this to check if the relevant information has actually been
  /// filled already and if not give some warning or print only some of the values.
  void Summarise() const {
    MACH3LOG_INFO("### Input File Summary ###");
    MACH3LOG_INFO("  Root file loc: {}", fileName);
    MACH3LOG_INFO("  Came from fitter: {}", fitter);
    MACH3LOG_INFO("  N LLH scans: {}", availableParams_LLH.size());
    MACH3LOG_INFO("  N Processed post fit errors: {}", availableParams_postFitErrors.size());
  }

  /// @brief Print out a more detailed summary of what is contained in the file.
  /// @todo Same as for Summarise(), could add a check that the file info has actually been filled.
  void Dump() const {
    Summarise();
    MACH3LOG_INFO("Available LLH parameters: ");
    for (std::string paramName : availableParams_LLH)
    {
      MACH3LOG_INFO("  {}", paramName);
    }
    MACH3LOG_INFO("Available Post Fit errors: ");
    for (std::string paramName : availableParams_postFitErrors)
    {
      MACH3LOG_INFO("  {}", paramName);
    }
  }

  std::shared_ptr<TFile> file;          //!< Pointer to the underlying file for this InputFile instance.
  std::string fileName; //!< The location of the underlying file.

  fitterEnum
      fitter; //!< Which fitter this file came from, detected by InputManager::fillFileInfo().

  /// @todo I think it would be nice to store all the InputFile data in some non root, maybe custom
  /// c++ types, and have separate reader classes for potential different input file formats that
  /// can convert to these types.
  /// @todo Currently all fitters use root but maybe some wonderful day
  /// in the far far future this wont be the case. Doing things this way we could maintain the
  /// possible future option of moving away from root???

  // EM: info on the LLH scans
  bool hasLLHScans; //!< Whether or not this file contains any log likelihood scans.
  std::vector<std::string>
      availableParams_LLH; //!< The parameters that this file contains likelihood scans for.
  std::unordered_map<std::string, bool>
      hasLLHScans_map; //!< Whether this file contains specific types of likelihood scans, e.g.
                       //!< hasLLHScans_map.at("prior").
  std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<TGraph>>>
      LLHScans_map; //!< The actual graphs of the likelihood scans, organised as
                    //!< LLHScans_map.at(LLHType).at(parameter).
  std::unordered_map<std::string, std::unordered_map<std::string, bool>>
      availableParams_map_LLH; //!< Whether this file has LLH scans of a particular type for a
                               //!< particular parameter, organised as
                               //!< availableParams_map_LLH.at(LLHType).at(parameter)
  // EM: maybe don't really need the "availableParams.." "hasLLHScans...", could just implement
  // these as functions that look for the given parameter in the existing maps, might make things
  // nicer

  // EM: info on LLH scans broken down by sample
  bool hasLLHScansBySample; //!< Whether or not this file contains log likelihood scans broken down
                            //!< per sample.
  std::vector<std::string> availableSamples_LLH; //!< The samples that this file contains individual
                                                 //!< likelihood scans for.
  std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<TGraph>>>
      LLHScansBySample_map; //!< The actual graphs of the likelihood scans, organised as
                            //!< LLHScans_map.at(sample).at(parameter).
  std::unordered_map<std::string, std::unordered_map<std::string, bool>>
      availableParams_map_LLHBySample; //!< Whether this file has LLH scans of a particular type for
                                       //!< a particular parameter, organised as
                                       //!< availableParams_map_LLH.at(sample).at(parameter)

  // EM: info on post fit errors
  bool hasPostFitErrors; //!< Whether or not this file contains any processed post fit errors.
  std::vector<std::string>
      availableParams_postFitErrors; //!< The parameters that this file has post fit errors for.
  std::unordered_map<std::string, std::unordered_map<std::string, double>>
      postFitErrors; //!< The actual errors on the parameters, organised as
                     //!< postFitErrors.at(errorType).at(parameter).
  std::unordered_map<std::string, std::unordered_map<std::string, double>>
      postFitValues; //!< The post fit values of the parameters, organised as
                     //!< postFitValues.at(errorType).at(parameter).
  std::string defaultErrorType;
  // TH1D *postFitErrors; //!< Pointer to the TH1D with named bins that contain the post fit value
  // and errors.

  // EM: almost certainly won't want to load all of these into memory at the start
  bool hasSigmaVars; //!< Whether or not this file contains Sigma variations.
};

/// @brief This guy talks to the input files and is intended to smooth over the fact that our OA
/// fitters use such a beautiful and diverse range of output file structures.
/// @details Intended use
/// for this class is to be called from inside of a PlottingManager object which is what the user is
/// intended to interact with, though in principle it should be strong and independent enough to be
/// used all by itself. The interface to this object is intentionally pretty small as it is intended
/// to deal with the complexities of the different input files pretty automatically. Usage should
/// follow something like this:
/// @code InputManager man(translationConfig.yaml) man.addFile(filename1)
/// man.addFile(filename2)
/// ...
///
/// then use various getter functions like GetLLHScan() to access the data in the files to make
/// plots.
/// @endcode
/// @todo A lot of string comparisons going on in the code for the post fit errors,
/// would maybe be good to implement some kind of indexing for parameters and just have a map
/// between the parameter names and this index to be used internally by the class.
/// @todo Add code to read MCMC from input file
/// @todo Add code to read Sigma variations. There are usually so many of
/// these that we might not want to read them all into memory at once when adding the file but maybe
/// just check to see what ones are in the file, then be able to read them on the fly later.
/// @todo Add code to read in Violin plots from the postfit error output files.
class InputManager {
public:
  const std::string NOT_FOUND_STR =
      "__PARAM_NAME_NOT_FOUND__"; //!< the default string to return if something can't be found

  /// @brief Construct a new InputManager using specified fitter translation config file.
  /// @param translationConfigName The config file defining the fitter file structures, fit
  /// parameter, and what the parameters are called in each fitter.
  /// @return Constructed InputManager instance.
  InputManager(const std::string &translationConfigName);

  /// @brief Add a new InputFile object to this input manager.
  /// @param fileName The name of the file to read.
  void addFile(const std::string &fileName);

  /// @brief Convert from fitterEnum to the name of the fitter.
  /// @param fitterId The fitter ID to convert.
  /// @return The name of the fitter, or "UNKNOWN_FITTER" if the fitterId does not match any known
  /// fitter.
  static const std::string convertFitterNames(fitterEnum fitterId) {
    switch (fitterId)
    {
    case kMaCh3ND:
      return "MaCh3_ND";
    case kMaCh3FD:
      return "MaCh3_FD";
    case kGundam:
      return "GUNDAM";
    case kBANFF:
      return "BANFF";
      // case kValor: return "Valor";
      // case kPTheta: return "PTheta";

    default:
      return "UNKNOWN_FITTER";
    }
  }

  /// @brief Convert from fileTypeEnum to the name of the file type.
  /// @param fileType The file type ID to convert.
  /// @return The name of the file type, or "UNKNOWN_FILE_TYPE" if the fileType does not match any
  /// known file type.
  static const std::string convertFileTypeNames(fileTypeEnum fileType) {
    switch (fileType)
    {
    case kLLH:
      return "LLH";
    case kPostFit:
      return "PostFit";
    case kMarkovChain:
      return "MarkovChain";
    case kSigmaVar:
      return "SigmaVar";

    default:
      return "UNKNOWN_FILE_TYPE";
    }
  }

  /// @brief Print out what this Inputmanager instance knows about.
  /// @param printLevel The level of detail to go into in the summary printout: options are
  /// "summary", and "dump".
  void Print(const std::string &printLevel = "summary") const;

  // FNs to get to the tasty tasty data stored in the files

  /// @brief Get the log likelihood scan for a particular parameter from a particular input file.
  /// @param fileNum The index of the file that you would like to get the scan from.
  /// @param paramName The name of the parameter whose LLH scan you would like.
  /// @param LLHType The type of likelihood scan you would like, e.h. total, prior etc.
  /// @tparam T The type you would like your scan returned as, currently only TGraph and TH1D are
  /// supported
  /// @return The graph of the likelihood scan.
  TGraph GetLLHScan_TGraph(int fileNum, std::string paramName, std::string LLHType) const {
    if (!GetEnabled_LLH(fileNum, paramName, LLHType))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for parameter {}", fileNum, paramName);
      MACH3LOG_WARN("am returning an empty TGraph");
      return TGraph(); 
    }
    return *fileVec[fileNum].LLHScans_map.at(LLHType).at(paramName);
  }

  TH1D GetLLHScan_TH1D(int fileNum, std::string paramName, std::string LLHType) const {
    if (!GetEnabled_LLH(fileNum, paramName, LLHType))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for parameter {}", fileNum, paramName);
      MACH3LOG_WARN("am returning an empty TH1D");
      return TH1D(); 
    }
    return TGraphToTH1D(*fileVec[fileNum].LLHScans_map.at(LLHType).at(paramName));
  }

  TGraph GetSampleSpecificLLHScan_TGraph(int fileNum, std::string paramName,
                                         std::string sample) const {
    
    if (!GetEnabled_LLHBySample(fileNum, paramName, sample))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for sample {} for parameter {}", fileNum, sample, paramName);
      MACH3LOG_WARN("am returning an empty TGraph");
      return TGraph();
    }
    return *fileVec[fileNum].LLHScansBySample_map.at(sample).at(paramName);
  }

  TH1D GetSampleSpecificLLHScan_TH1D(int fileNum, std::string paramName, std::string sample) const {
    if (!GetEnabled_LLHBySample(fileNum, paramName, sample))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for sample {} for parameter {}", fileNum, sample, paramName);
      MACH3LOG_WARN("am returning an empty TH1D");
      return TH1D();
    }
    return TGraphToTH1D(*fileVec[fileNum].LLHScansBySample_map.at(sample).at(paramName));
  }

  /// @brief Get the log likelihood scan for a particular parameter, for a specific sample, from a
  /// particular input file.
  /// @param fileNum The index of the file that you would like to get the scan from.
  /// @param paramName The name of the parameter whose LLH scan you would like.
  /// @param sample The sample that you would like the LLH scan for.
  /// @tparam T The type you would like your scan returned as, currently only TGraph and TH1D are
  /// supported
  /// @return The graph of the likelihood scan.
  template <typename T = TGraph>
  // Default template (i.e. the template specified was not one of the ones implemented below)
  T GetSampleSpecificLLHScan(int fileNum, std::string paramName, std::string sample) const {
    throw MaCh3Exception(__FILE__, __LINE__, "uuuuuh sorry but im not sure how to convert likelihood scans to the specified type");
  }

  /// @brief Get whether or not a particular parameter has an LLH scan in a particular input file.
  /// @param fileNum The index of the file that you would like to know about.
  /// @param paramName The name of the parameter whose LLH scan you would like to check for.
  /// @param LLHType The type of likelihood scan you would like to check for, e.h. total, prior etc.
  /// @return true if scan exists, false if not.
  inline bool GetEnabled_LLH(int fileNum, std::string paramName,
                             std::string LLHType = "total") const {
    return fileVec[fileNum].availableParams_map_LLH.at(LLHType).at(paramName);
  }

  /// @brief Get whether or not a particular parameter has an LLH scan in a particular input file
  /// for a particular sample.
  /// @param fileNum The index of the file that you would like to know about.
  /// @param paramName The name of the parameter whose LLH scan you would like to check for.
  /// @param sample The sample to check.
  /// @return true if scan exists, false if not.
  inline bool GetEnabled_LLHBySample(int fileNum, std::string paramName, std::string sample) const {
    return fileVec[fileNum].availableParams_map_LLHBySample.at(sample).at(paramName);
  }

  /// @brief Get the post fit error for a particular parameter from a particular input file.
  /// @param fileNum The index of the file that you would like to get the value from.
  /// @param paramName The name of the parameter whose error you would like.
  /// @param errorType The type of error to return, usually "prior" will be defined, but other
  /// possible types will be fitter dependent, e.g. "gauss" or "hpd" for MaCh3. If not specified,
  /// will use the default one, as specified in the fitter definition config.
  /// @return The error on the specified parameter.
  const float GetPostFitError(int fileNum, const std::string &paramName, std::string errorType = "") const;

  /// @brief Get the post fit value for a particular parameter from a particular input file.
  /// @param fileNum The index of the file that you would like to get the value from.
  /// @param paramName The name of the parameter whose value you would like.
  /// @param errorType The type of error to return, usually "prior" will be defined, but other
  /// possible types will be fitter dependent, e.g. "gauss" or "hpd" for MaCh3. If not specified,
  /// will use the default one, as specified in the fitter definition config.
  /// @return The value of the specified parameter.
  const float GetPostFitValue(int fileNum, const std::string &paramName, std::string errorType = "") const;

  /// @name General Getters
  /// @{
  inline const std::vector<std::string> &GetKnownParameters() const { return knownParameters; }
  inline const std::vector<std::string> &GetKnownSamples() const { return knownSamples; }
  inline const int GetNInputFiles() const { return fileVec.size(); }
  /// @}

  /// @name File Specific Getters
  /// @{
  inline InputFile GetFile(int fileId) const { return fileVec[fileId]; }
  inline std::string TranslateName(int fileId, fileTypeEnum fileType, std::string paramName) const {
    return getFitterSpecificParamName(fileVec[fileId].fitter, fileType, paramName);
  }
  inline const std::vector<std::string> &GetKnownLLHParameters(int fileId) const {
    return fileVec[fileId].availableParams_LLH;
  }
  inline const std::vector<std::string> &GetKnownLLHSamples(int fileId) const {
    return fileVec[fileId].availableSamples_LLH;
  }
  inline const std::vector<std::string> &GetKnownPostFitParameters(int fileId) const {
    return fileVec[fileId].availableParams_postFitErrors;
  }
  /// @}

private:
  // Helper function to parse a root file location string
  // any instance of {PARAMETER} gets replaced with fitter specific parameter name
  // similar for {SAMPLE}
  // any other token that gets added in future can also just be added here and added to end of the
  // argument list The string will also be split by the ":" delimeter into a directory and a name to
  // look for in that directory
  std::vector<std::string> parseLocation(const std::string &locationString, fitterEnum fitter,
                                         fileTypeEnum fileType, const std::string &parameter = "",
                                         const std::string &sample = "") const;

  // helper function to look for an object defined by a vector of strings [ (directory to look in),
  // (end of the name of the object) ] this just calls roots TDirectoryFile->Get() function so
  // behaves in the same way:
  //   if the object exists then it is loaded and returned
  //   if it doesn't exist then just returns nullptr
  std::shared_ptr<TObject> findRootObject(const InputFile &fileDef, const std::vector<std::string> &locationVec) const;

  // Check the input file for a particular post fit error for a particular parameter. if
  // setInputFileError, will add the information to the InputFiles postFitErrors histogram.
  bool findPostFitParamError(InputFile &inputFileDef, const std::string &parameter, fitterEnum fitter,
                             const std::string &errorType, bool setInputFileError = false);

  // Check the input file for a particular post fit error for a particular parameter. if
  // setInputFileError, will add the information to the InputFiles postFitErrors histogram.
  bool findBySampleLLH(InputFile &inputFileDef, const std::string &parameter, fitterEnum fitter,
                       const std::string &sample, bool setInputFileScan = false);

  // fns tp read an input file
  void fillFileInfo(InputFile &inputFileDef, bool printThoughts = true);
  void fillFileData(InputFile &inputFileDef, bool printThoughts = true);

  // helper function to read from the translation config file to get an option for a particular sub
  // node (e.g. Parameters or Samples) returns true and sets "ret" if the option is specified for this
  // fitter and parameter otherwise returns false
  template <typename T>
  inline bool getFitterSpecificOption(fitterEnum fitter, std::string option, T &ret, const std::string &parameter,
                               YAML::Node subConfig) const{
    if (subConfig[parameter])
    {

      // EM: this is config definition of fitter specific names for this parameter
      YAML::Node paramTranslation = subConfig[parameter];

      if (paramTranslation[convertFitterNames(fitter)])
      {
        // EM: then this is definition of how to find parameter in the specified fitter
        YAML::Node fitterParamTranslation = paramTranslation[convertFitterNames(fitter)];

        if (fitterParamTranslation[option])
        {
          ret = fitterParamTranslation[option].as<T>();
          return true;
        }
      }
    }
    return false;
  }

  // Specialised option getter for parameters
  template <typename T>
  inline bool getFitterSpecificParamOption(fitterEnum fitter, const std::string &option, T &ret,
                                           const std::string &parameter) const {
    return getFitterSpecificOption<T>(fitter, option, ret, parameter, _parametersConfig);
  }

  // specialised option getter for samples
  template <typename T>
  inline bool getFitterSpecificSampleOption(fitterEnum fitter, std::string option, T &ret,
                                            std::string parameter) const {
    return getFitterSpecificOption<T>(fitter, option, ret, parameter, _samplesConfig);
  }

  // helper function to read from the translation config to get the parameter name for a specific
  // fitter and file type
  inline std::string getFitterSpecificParamName(fitterEnum fitter, fileTypeEnum fileType,
                                                const std::string &parameter) const {
    std::string specificName;
    if (getFitterSpecificParamOption<std::string>(fitter, convertFileTypeNames(fileType),
                                                  specificName, parameter))
    {
      return specificName;
    }
    // EM: Default to just return the specified name
    return parameter;
  }

  // helper function to read from the translation config file to get the sample name for a specific
  // fitter and file type
  inline std::string getFitterSpecificSampleName(fitterEnum fitter, fileTypeEnum fileType,
                                                 const std::string &sample) const {
    std::string specificName;
    if (getFitterSpecificSampleOption<std::string>(fitter, convertFileTypeNames(fileType),
                                                   specificName, sample))
    {
      return specificName;
    }
    // EM: Default to just return the specified name
    return sample;
  }

  // helper fn to test if string "str" ends with other string "ending"
  inline bool strEndsWith(std::string str, std::string ending) const {
    uint pos = str.find(ending);
    return (pos == str.length() - ending.length());
  }

  // vector of InputFile objects that this manager is responsible for
  std::vector<InputFile> fileVec;

  // all parameters which are known to this InputManager: all the ones defined in the translation
  // config used to create it
  std::vector<std::string> knownParameters;

  // all samples which are known to this InputManager: all the ones defined in the translation
  // config used to create it
  std::vector<std::string> knownSamples;

  // map the enum to actual fitter names
  std::map<fitterEnum, std::string> fitterNames;

private:
  // the configs defining the translation of parameters between fitters and also the directory
  // structure for each fitter
  YAML::Node _translatorConfig;
  YAML::Node _fitterSpecConfig;
  YAML::Node _parametersConfig;
  YAML::Node _samplesConfig;
};
} // namespace MaCh3Plotting