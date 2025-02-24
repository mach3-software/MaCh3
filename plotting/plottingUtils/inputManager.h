#pragma once

// MaCh3 Includes
#include "manager/MaCh3Logger.h"
#include "manager/YamlHelper.h"
#include "manager/MaCh3Exception.h"
#include "mcmc/MCMCProcessor.h"
#include "SampleHandler/Structs.h"

// Other plotting includes
#include "plottingUtils.h"

namespace MaCh3Plotting {

/// @brief Types of possible file that can be read.
enum fileTypeEnum {
  kLLH,         //!< Log Likelihood scan
  kPostFit,     //!< Processed post fit errors
  kMCMC,        //!< MCMC chain
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

  // NO COPYING!!
  InputFile( const InputFile& ) = delete;
  // Moving ok
  InputFile( InputFile&& ) = default;
  // EM: ^^ there should only really be one instance of an InputFile for each
  // underlying root file. If copying is allowed then it can lead to bad bad things
  // like pointers getting deleted and then trying to be accessed by a copied 
  // InputFile instance

  /// @brief Destructor.
  ~InputFile() {
    MACH3LOG_DEBUG("###### Deleting InputFile Object holding file ######");
    LLHScans_map.clear();
  }

  /// @brief Close out the underlying root file
  /// @details Should only be done once this InputFile is *done* with, should only really ever be done by the InputManager that holds this object 
  void Close(){
    MACH3LOG_DEBUG("[InputFile] closing file {}", fileName);
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

  /// ptr to an MCMCProcessor instance to be used if this is a MaCh3 input file
  MCMCProcessor *mcmcProc = nullptr;
  TTree *posteriorTree = nullptr;

  std::shared_ptr<TFile> file;          //!< Pointer to the underlying file for this InputFile instance.
  std::string fileName; //!< The location of the underlying file.

  std::string fitter; //!< Which fitter this file came from, detected by InputManager::fillFileInfo().

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

  // Stuff relating to MCMC
  /// Whether or not the file has processed 1d posteriors
  bool has1dPosteriors;
  /// Whether or not the file has unprocessed MCMC chain steps
  bool hasMCMCchain;

  /// The number of steps in the MCMC chain
  int nMCMCentries;

  /// whether or not specific parameters exist in the MCMC posterior chain
  std::unordered_map<std::string, bool> availableParams_map_MCMCchain;
  std::unordered_map<std::string, bool> availableParams_map_1dPosteriors;

  std::vector<std::string> availableParams_1dPosteriors;
  std::vector<std::string> availableParams_MCMCchain;

  std::unordered_map<std::string, double*> MCMCstepParamsMap;
  std::unordered_map<std::string, int> MCMCstepTreeIndicesMap;

  std::unordered_map<std::string, std::shared_ptr<TGraph>> posteriors1d_map;
  
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
/// then use various getter functions like getLLHScan() to access the data in the files to make
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

  /// @brief Destructor
  /// @details Close out all the files that the manager is responsible for
  ~InputManager()
  {
    MACH3LOG_DEBUG("##### Deleting InputManager Instance #####");
    for (InputFile &file: _fileVec)
    {
      file.Close();
    }

    _fileVec.clear();
  }

  // NO COPYING!
  InputManager(const InputManager&) = delete;
  // moving is ok
  InputManager(InputManager&&) = default;
  // EM: ^^ Do this as there should only ever really be one instance of the 
  // InputManager (should maybe make it a singleton actually)
  // if copied then it can lead to things being deleted that we really don't
  // want to be deleted

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
    case kMCMC:
      return "MCMC";
    case kSigmaVar:
      return "SigmaVar";
    case kNFileTypes:
      return "NFileTypes";
    default:
      return "UNKNOWN_FILE_TYPE";
    }
  }

  /// @brief Print out what this Inputmanager instance knows about.
  /// @param printLevel The level of detail to go into in the summary printout: options are
  /// "summary", and "dump".
  void print(const std::string &printLevel = "summary") const;

  // FNs to get to the tasty tasty data stored in the files

  /// @brief Get the log likelihood scan data for a particular parameter from a particular input file.
  /// @param fileNum The index of the file you want the data from.
  /// @param paramName The parameter you want the information about.
  /// @param LLHType The type of log likelihood scan you want (e.g. total, penalty, etc.)
  /// @return A vector of vectors containing the LLH scan data. First entry is x axis, 2nd is y axis
  std::vector<std::vector<double>> getLLHScan(int fileNum, std::string paramName, std::string LLHType) const {
    if (!getEnabledLLH(fileNum, paramName, LLHType))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for parameter {}", fileNum, paramName);
      MACH3LOG_WARN("am returning an empty vector");
      return std::vector<std::vector<double>>(2); 
    }
    return TGraphToVector(*_fileVec[fileNum].LLHScans_map.at(LLHType).at(paramName));
  }

  /// @brief Get the MCMC chain entry in an InputFile.
  /// @param fileNum The index of the file you want the data from.
  /// @param entry The entry to get.
  void getMCMCentry(int fileNum, int entry) const {
    // EM: the const here is a little bit of a lie since GetEntry will in fact modify 
    //     the pointers to the stored data for the MCMC chain values but hopefully this should be ok
    const InputFile &file = _fileVec[fileNum];

    if( entry > file.nMCMCentries )
    {
      MACH3LOG_ERROR("Trying to access entries beyond what exist in file {}. No-can-do!", file.fileName);
    }
    else
    {
      MACH3LOG_TRACE("Getting entry {} in MCMC tree for file at index {}", entry, fileNum);
      file.posteriorTree->GetEntry(entry);
      MACH3LOG_TRACE("  Got successfuly");
    }
  }

  /// @brief Get the parameter value for the current step for a particular parameter from a particular input file.
  /// @param fileNum The index of the file you want the data from.
  /// @param paramName The parameter you want the value of.
  double getMCMCvalue(int fileNum, std::string paramName) const {
    if (!getEnabledMCMCchain(fileNum, paramName))
    {
      MACH3LOG_WARN("file at index {} does not have an MCMC entry for parameter {}", fileNum, paramName);
      MACH3LOG_WARN("am returning a bad float");
      return M3::_BAD_DOUBLE_;
    }

    return *_fileVec[fileNum].MCMCstepParamsMap.at(paramName);
  }

  /// @brief Get the 1d posterior particular parameter from a particular input file.
  /// @param fileNum The index of the file you want the data from.
  /// @param paramName The parameter you want the information about.
  /// @return A vector of vectors containing the posterior data. First entry is x axis (i.e. the parameter values), 2nd is y axis
  std::vector<std::vector<double>> get1dPosterior(int fileNum, std::string paramName) const {
    if (!getEnabled1dPosteriors(fileNum, paramName))
    {
      MACH3LOG_WARN("file at index {} does not have a 1d posterior for parameter {}", fileNum, paramName);
      MACH3LOG_WARN("am returning an empty vector");
      return std::vector<std::vector<double>>(2);
    }
    return TGraphToVector(*_fileVec[fileNum].posteriors1d_map.at(paramName));
  }

  /// @brief Get the number of entries in the MCMC chain in a particular file.
  /// @param fileNum The index of the file you want the number of steps from.
  int getnMCMCentries(int fileNum) const {
    return _fileVec[fileNum].nMCMCentries;
  }

  /// @brief Get the log likelihood scan for a particular parameter from a particular input file.
  /// @param fileNum The index of the file that you would like to get the scan from.
  /// @param paramName The name of the parameter whose LLH scan you would like.
  /// @param LLHType The type of likelihood scan you would like, e.h. total, prior etc.
  /// @tparam T The type you would like your scan returned as, currently only TGraph and TH1D are
  /// supported
  /// @return The graph of the likelihood scan.
  inline TGraph getLLHScan_TGraph(int fileNum, std::string paramName, std::string LLHType) const {
    if (!getEnabledLLH(fileNum, paramName, LLHType))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for parameter {}", fileNum, paramName);
      MACH3LOG_WARN("am returning an empty TGraph");
      return TGraph(); 
    }
    return *_fileVec[fileNum].LLHScans_map.at(LLHType).at(paramName);
  }

  inline TH1D getLLHScan_TH1D(int fileNum, std::string paramName, std::string LLHType) const {
    if (!getEnabledLLH(fileNum, paramName, LLHType))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for parameter {}", fileNum, paramName);
      MACH3LOG_WARN("am returning an empty TH1D");
      return TH1D(); 
    }
    return TGraphToTH1D(*_fileVec[fileNum].LLHScans_map.at(LLHType).at(paramName));
  }

  /// @brief Get the log likelihood scan for a particular parameter, for a specific sample, from a
  /// particular input file.
  /// @param fileNum The index of the file that you would like to get the scan from.
  /// @param paramName The name of the parameter whose LLH scan you would like.
  /// @param sample The sample that you would like the LLH scan for.
  /// @return A vector of vectors containing the LLH scan data. First entry is x axis, 2nd is y axis.
  std::vector<std::vector<double>> getSampleSpecificLLHScan(int fileNum, std::string paramName, std::string sample) const {
    if (!getEnabledLLHBySample(fileNum, paramName, sample))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for sample {} for parameter {}", fileNum, sample, paramName);
      MACH3LOG_WARN("am returning an empty vector");
      return std::vector<std::vector<double>>(2); 
    }
    return TGraphToVector(*_fileVec[fileNum].LLHScansBySample_map.at(sample).at(paramName));
  }

  inline TGraph getSampleSpecificLLHScan_TGraph(int fileNum, std::string paramName,
                                         std::string sample) const {
    if (!getEnabledLLHBySample(fileNum, paramName, sample))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for sample {} for parameter {}", fileNum, sample, paramName);
      MACH3LOG_WARN("am returning an empty TGraph");
      return TGraph();
    }
    return *_fileVec[fileNum].LLHScansBySample_map.at(sample).at(paramName);
  }

  inline TH1D getSampleSpecificLLHScan_TH1D(int fileNum, std::string paramName, std::string sample) const {
    if (!getEnabledLLHBySample(fileNum, paramName, sample))
    {
      MACH3LOG_WARN("file at index {} does not have LLH scan for sample {} for parameter {}", fileNum, sample, paramName);
      MACH3LOG_WARN("am returning an empty TH1D");
      return TH1D();
    }
    return TGraphToTH1D(*_fileVec[fileNum].LLHScansBySample_map.at(sample).at(paramName));
  }

  /// @brief Get whether or not a particular parameter has an LLH scan in a particular input file.
  /// @param fileNum The index of the file that you would like to know about.
  /// @param paramName The name of the parameter whose LLH scan you would like to check for.
  /// @param LLHType The type of likelihood scan you would like to check for, e.h. total, prior etc.
  /// @return true if scan exists, false if not.
  inline bool getEnabledLLH(int fileNum, std::string paramName,
                             std::string LLHType = "total") const {
    return _fileVec[fileNum].availableParams_map_LLH.at(LLHType).at(paramName);
  }

  /// @brief Get whether or not a particular parameter has an LLH scan in a particular input file
  /// for a particular sample.
  /// @param fileNum The index of the file that you would like to know about.
  /// @param paramName The name of the parameter whose LLH scan you would like to check for.
  /// @param sample The sample to check.
  /// @return true if scan exists, false if not.
  inline bool getEnabledLLHBySample(int fileNum, std::string paramName, std::string sample) const {
    return _fileVec[fileNum].availableParams_map_LLHBySample.at(sample).at(paramName);
  }

  /// @brief Get whether or not a particular parameter has MCMC chain entries in a particular input file.
  /// @param fileNum The index of the file that you would like to know about.
  /// @param paramName The name of the parameter whose LLH scan you would like to check for.
  /// @return true if scan exists, false if not.
  inline bool getEnabledMCMCchain(int fileNum, std::string paramName) const {
    return _fileVec[fileNum].availableParams_map_MCMCchain.at(paramName);
  }

  /// @brief Get whether or not a particular parameter has 1d posteriors in a particular input file.
  /// @param fileNum The index of the file that you would like to know about.
  /// @param paramName The name of the parameter whose 1d posterior you would like to check for.
  /// @return true if scan exists, false if not.
  inline bool getEnabled1dPosteriors(int fileNum, std::string paramName) const {
    return _fileVec[fileNum].availableParams_map_1dPosteriors.at(paramName);
  }

  /// @brief Get the post fit error for a particular parameter from a particular input file.
  /// @param fileNum The index of the file that you would like to get the value from.
  /// @param paramName The name of the parameter whose error you would like.
  /// @param errorType The type of error to return, usually "prior" will be defined, but other
  /// possible types will be fitter dependent, e.g. "gauss" or "hpd" for MaCh3. If not specified,
  /// will use the default one, as specified in the fitter definition config.
  /// @return The error on the specified parameter.
  double getPostFitError(int fileNum, const std::string &paramName, std::string errorType = "") const;

  /// @brief Get the post fit value for a particular parameter from a particular input file.
  /// @param fileNum The index of the file that you would like to get the value from.
  /// @param paramName The name of the parameter whose value you would like.
  /// @param errorType The type of error to return, usually "prior" will be defined, but other
  /// possible types will be fitter dependent, e.g. "gauss" or "hpd" for MaCh3. If not specified,
  /// will use the default one, as specified in the fitter definition config.
  /// @return The value of the specified parameter.
  double getPostFitValue(int fileNum, const std::string &paramName, std::string errorType = "") const;

  /// @name General Getters
  /// @{
  inline const std::vector<std::string> &getKnownParameters() const { return _knownParameters; }
  inline const std::vector<std::string> &getKnownSamples() const { return _knownSamples; }
  inline size_t getNInputFiles() const { return _fileVec.size(); }

  /// @brief Get all parameters which have some set of tags
  /// @param tags The tags to check for 
  /// @param checkType The type of check to perform:
  /// checkType can be one of:
  ///   - all: All parameters which have *all* of the specified tags will be returned 
  ///   - any: All parameters which have *any* of the specified tags will be returned 
  ///   - exact: All of the parameters which have *exactly* the specified tags will be returned
  inline std::vector<std::string> getTaggedParameters(const std::vector<std::string> &tags, std::string checkType = "all") const {
    return getTaggedValues(_knownParameters, _paramToTagsMap, tags, checkType);
  }

  /// @brief Get all samples which have some set of tags
  /// @param tags The tags to check for 
  /// @param checkType The type of check to perform:
  /// checkType can be one of:
  ///   - all: All samples which have *all* of the specified tags will be returned 
  ///   - any: All samples which have *any* of the specified tags will be returned 
  ///   - exact: All of the samples which have *exactly* the specified tags will be returned
  inline std::vector<std::string> getTaggedSamples(const std::vector<std::string> &tags, std::string checkType = "all") const {
    return getTaggedValues(_knownSamples, _sampleToTagsMap, tags, checkType);
  }
  /// @}

  /// @name File Specific Getters
  /// @{
  inline InputFile const &getFile(int fileId) const { return _fileVec[fileId]; }
  inline std::string translateName(int fileId, fileTypeEnum fileType, std::string paramName) const {
    return getFitterSpecificParamName(_fileVec[fileId].fitter, fileType, paramName);
  }
  inline const std::vector<std::string> &getKnownLLHParameters(int fileId) const {
    return _fileVec[fileId].availableParams_LLH;
  }
  inline const std::vector<std::string> &getKnownLLHSamples(int fileId) const {
    return _fileVec[fileId].availableSamples_LLH;
  }
  inline const std::vector<std::string> &getKnownPostFitParameters(int fileId) const {
    return _fileVec[fileId].availableParams_postFitErrors;
  }
  inline const std::vector<std::string> &getKnownMCMCParameters(int fileId) const {
    return _fileVec[fileId].availableParams_MCMCchain;
  }
  inline const std::vector<std::string> &getKnown1dPosteriorParameters(int fileId) const {
    return _fileVec[fileId].availableParams_1dPosteriors;
  }
  /// @}

private:
  // Helper function to get tagged values from a vector of values
  // specify the initial list of *values*, the map of values to their tags, the tags to check,
  // and the type of check to perform (see getTaggedParameter() for details)
  std::vector<std::string> getTaggedValues(const std::vector<std::string> &values, 
                                           const std::unordered_map<std::string, 
                                           std::vector<std::string>> &tagMap, 
                                           const std::vector<std::string> &tags, std::string checkType) const;

  // Helper function to parse a root file location string
  // any instance of {PARAMETER} gets replaced with fitter specific parameter name
  // similar for {SAMPLE}
  // any other token that gets added in future can also just be added here and added to end of the
  // argument list The string will also be split by the ":" delimeter into a directory and a name to
  // look for in that directory
  std::vector<std::string> parseLocation(const std::string &locationString, std::string &fitter,
                                         fileTypeEnum fileType, const std::string &parameter = "",
                                         const std::string &sample = "", const std::string &parameter2 = "") const;

  // helper function to look for an object defined by a vector of strings [ (directory to look in),
  // (end of the name of the object) ] this just calls roots TDirectoryFile->Get() function so
  // behaves in the same way:
  //   if the object exists then it is loaded and returned
  //   if it doesn't exist then just returns nullptr
  std::shared_ptr<TObject> findRootObject(const InputFile &fileDef, const std::vector<std::string> &locationVec) const;

  // Check the input file for a particular post fit error for a particular parameter. if
  // setInputFileError, will add the information to the InputFiles postFitErrors histogram.
  bool findPostFitParamError(InputFile &inputFileDef, const std::string &parameter, std::string &fitter,
                             const std::string &errorType, bool setInputFileError = false);

  // Check the input file for a particular post fit error for a particular parameter. if
  // setInputFileError, will add the information to the InputFiles postFitErrors histogram.
  bool findBySampleLLH(InputFile &inputFileDef, const std::string &parameter, std::string &fitter,
                       const std::string &sample, bool setInputFileScan = false);

  // check the input file for raw MCMC step values for a particular parameter
  bool findRawChainSteps(InputFile &inputFileDef, const std::string &parameter, std::string &fitter, bool setInputBranch = false ) const ;

  // check the input file for processed 1d posteriors for a particular parameter
  bool find1dPosterior(InputFile &inputFileDef, const std::string &parameter, std::string &fitter, bool setFileData = false) const ;

  // fns tp read an input file
  void fillFileInfo(InputFile &inputFileDef, bool printThoughts = true);
  void fillFileData(InputFile &inputFileDef, bool printThoughts = true);

  // helper function to read from the translation config file to get an option for a particular sub
  // node (e.g. Parameters or Samples) returns true and sets "ret" if the option is specified for this
  // fitter and parameter otherwise returns false
  template <typename T>
  inline bool getFitterSpecificOption(const std::string &fitter, const std::string &option, T &ret, const std::string &parameter,
                               YAML::Node subConfig) const{
    if (subConfig[parameter])
    {
      // EM: this is config definition of fitter specific names for this parameter
      YAML::Node paramTranslation = subConfig[parameter];

      if (paramTranslation[fitter])
      {
        // EM: then this is definition of how to find parameter in the specified fitter
        YAML::Node fitterParamTranslation = paramTranslation[fitter];

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
  inline bool getFitterSpecificParamOption(const std::string &fitter, const std::string &option, T &ret,
                                           const std::string &parameter) const {
    return getFitterSpecificOption<T>(fitter, option, ret, parameter, _parametersConfig);
  }

  // specialised option getter for samples
  template <typename T>
  inline bool getFitterSpecificSampleOption(const std::string &fitter, std::string option, T &ret,
                                            std::string parameter) const {
    return getFitterSpecificOption<T>(fitter, option, ret, parameter, _samplesConfig);
  }

  // helper function to read from the translation config to get the parameter name for a specific
  // fitter and file type
  inline std::string getFitterSpecificParamName(const std::string &fitter, fileTypeEnum fileType,
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
  inline std::string getFitterSpecificSampleName(const std::string &fitter, fileTypeEnum fileType,
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
  inline bool strEndsWith(const std::string& str, const std::string& ending) const {
    if (str.size() >= ending.size()) {
      return str.compare(str.size() - ending.size(), ending.size(), ending) == 0;
    }
    return false;
  }

private:
  // all parameters which are known to this InputManager: all the ones defined in the translation
  // config used to create it
  std::vector<std::string> _knownParameters;

  // all samples which are known to this InputManager: all the ones defined in the translation
  // config used to create it
  std::vector<std::string> _knownSamples;

  // hold the names of the fitters known to this manager
  std::vector<std::string> _knownFitters;

  // map parameter names to their specified tags
  std::unordered_map<std::string, std::vector<std::string>> _paramToTagsMap;
  // map sample names to their specified tags
  std::unordered_map<std::string, std::vector<std::string>> _sampleToTagsMap;
  
  // the configs defining the translation of parameters between fitters and also the directory
  // structure for each fitter
  YAML::Node _translatorConfig;
  YAML::Node _fitterSpecConfig;
  YAML::Node _parametersConfig;
  YAML::Node _samplesConfig;

  // vector of InputFile objects that this manager is responsible for
  std::vector<InputFile> _fileVec;
};
} // namespace MaCh3Plotting
