#include "Parameters/ParameterHandlerGeneric.h"

// ********************************************
// ETA - YAML constructor
// this will replace the root file constructor but let's keep it in
// to do some validations
ParameterHandlerGeneric::ParameterHandlerGeneric(const std::vector<std::string>& YAMLFile, std::string name,
                                                 double threshold, int FirstPCA, int LastPCA) {
// ********************************************
  MACH3LOG_INFO("Constructing instance of ParameterHandler using");
  inputFile = YAMLFile[0];
  matrixName = name;
  pca = true;

  doSpecialStepProposal = false;
  // Not using adaptive by default
  use_adaptive = false;
  for(unsigned int i = 0; i < YAMLFile.size(); i++)
  {
    MACH3LOG_INFO("{}", YAMLFile[i]);
  }
  MACH3LOG_INFO("as an input");

  if (threshold < 0 || threshold >= 1) {
    MACH3LOG_INFO("Principal component analysis but given the threshold for the principal components to be less than 0, or greater than (or equal to) 1. This will not work");
    MACH3LOG_INFO("Please specify a number between 0 and 1");
    MACH3LOG_INFO("You specified: ");
    MACH3LOG_INFO("Am instead calling the usual non-PCA constructor...");
    pca = false;
  }

  InitialiseFromConfig(YAMLFile);
  // Call the innocent helper function
  if (pca) ConstructPCA(threshold, FirstPCA, LastPCA);

  InitParametersTypeFromConfig();

  //ETA - again this really doesn't need to be hear...
  for (int i = 0; i < _fNumPar; i++)
  {
    // Sort out the print length
    if(int(_fNames[i].length()) > PrintLength) PrintLength = int(_fNames[i].length());
  } // end the for loop

  MACH3LOG_DEBUG("Constructing instance of ParameterHandler");
  InitParameters();
  // Print
  Print();
}

// ********************************************
void ParameterHandlerGeneric::LoadAndMergeYAML(const std::vector<std::string>& YAMLFile,
                                               std::map<std::pair<int,int>, std::unique_ptr<TMatrixDSym>>& ThrowSubMatrixOverrides) {
// ********************************************
  int running_num_file_pars = 0;

  _fYAMLDoc["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for(unsigned int i = 0; i < YAMLFile.size(); i++)
  {
    YAML::Node YAMLDocTemp = M3OpenConfig(YAMLFile[i]);

    if (YAMLDocTemp["ThrowMatrixOverride"]) { // LP: this allows us to put in
      // proposal matrix overrides per
      // parameter-containing file, add
      // the block diagonal proposal
      // matrix to a list and overwrite
      // the throw matrix after set up.
      auto filename =
      YAMLDocTemp["ThrowMatrixOverride"]["file"].as<std::string>();
      TFile *submatrix_file = TFile::Open(filename.c_str());

      auto matrixname =
      YAMLDocTemp["ThrowMatrixOverride"]["matrix"].as<std::string>();
      std::unique_ptr<TMatrixDSym> submatrix{
        submatrix_file->Get<TMatrixDSym>(matrixname.c_str())};
        if (!submatrix) {
          MACH3LOG_CRITICAL("Covariance matrix {} doesn't exist in file: {}",
                            matrixname, filename);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        auto numrows = submatrix->GetNrows();
        // LP: the -1 here is because we specify the last index for consistency
        // with PCAHandler, not the first index after the end as is more common
        // throughout computer science...
        ThrowSubMatrixOverrides[{running_num_file_pars,
          running_num_file_pars + (numrows - 1)}] = std::move(submatrix);

          // LP: check names by default, but have option to disable check if you
          // know what you're doing
          if (!bool(YAMLDocTemp["ThrowMatrixOverride"]["check_names"]) ||
            YAMLDocTemp["ThrowMatrixOverride"]["check_names"].as<bool>()) {
            auto nametree = submatrix_file->Get<TTree>("param_names");
          if (!nametree) {
            MACH3LOG_CRITICAL("TTree param_names doesn't exist in file: {}. Set "
            "ThrowMatrixOverride: {{ check_names: False }} to "
            "disable this check.",
            filename);
            throw MaCh3Exception(__FILE__, __LINE__);
          }
          std::string *param_name = nullptr;
          nametree->SetBranchAddress("name", &param_name);

          if (nametree->GetEntries() != int(YAMLDocTemp["Systematics"].size())) {
            MACH3LOG_CRITICAL("TTree param_names in file: {} has {} entries, but "
            "the corresponding yaml file only declares {} "
            "parameters. Set ThrowMatrixOverride: {{ "
              "check_names: False }} to disable this check.",
              filename, nametree->GetEntries(),
                              YAMLDocTemp["Systematics"].size());
            throw MaCh3Exception(__FILE__, __LINE__);
          }

          int pit = 0;
          for (const auto &param : YAMLDocTemp["Systematics"]) {
            nametree->GetEntry(pit++);
            auto yaml_pname = Get<std::string>(
              param["Systematic"]["Names"]["FancyName"], __FILE__, __LINE__);
            if ((*param_name) != yaml_pname) {
              MACH3LOG_CRITICAL(
                "TTree param_names in file: {} at entry {} has parameter {}, "
                "but "
                "the corresponding yaml parameter is named {}. Set "
                "ThrowMatrixOverride: {{ "
                  "check_names: False }} to disable this check.",
                  filename, pit, (*param_name), yaml_pname);
              throw MaCh3Exception(__FILE__, __LINE__);
            }
          }
            }
            submatrix_file->Close();
    }

    for (const auto& item : YAMLDocTemp["Systematics"]) {
      _fYAMLDoc["Systematics"].push_back(item);
      running_num_file_pars++;
    }
  }
}

// ********************************************
void ParameterHandlerGeneric::LoadCorrelationFromConfig(std::vector<std::map<std::string,double>>& Correlations,
                                                        std::map<std::string, int>& CorrNamesMap) {
// ********************************************
  // ETA Now that we've been through all systematic let's fill the covmatrix
  //This makes the root TCov from YAML
  TMatrixDSym* _fCovMatrix = new TMatrixDSym(_fNumPar);
  for(int j = 0; j < _fNumPar; j++) {
    (*_fCovMatrix)(j, j) = _fError[j]*_fError[j];
    //Get the map of parameter name to correlation from the Correlations object
    for (auto const& pair : Correlations[j]) {
      auto const& key = pair.first;
      auto const& val = pair.second;
      int index = -1;
      //If you found the parameter name then get the index
      if (CorrNamesMap.find(key) != CorrNamesMap.end()) {
        index = CorrNamesMap[key];
      } else {
        MACH3LOG_ERROR("Parameter {} not in list! Check your spelling?", key);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
      double Corr1 = val;
      double Corr2 = 0;
      if(Correlations[index].find(_fFancyNames[j]) != Correlations[index].end()) {
        Corr2 = Correlations[index][_fFancyNames[j]];
        //Do they agree to better than float precision?
        if(std::abs(Corr2 - Corr1) > FLT_EPSILON) {
          MACH3LOG_ERROR("Correlations are not equal between {} and {}", _fFancyNames[j], key);
          MACH3LOG_ERROR("Got : {} and {}", Corr2, Corr1);
          throw MaCh3Exception(__FILE__ , __LINE__ );
        }
      } else {
        MACH3LOG_ERROR("Correlation does not appear reciprocally between {} and {}", _fFancyNames[j], key);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
      (*_fCovMatrix)(j, index)= (*_fCovMatrix)(index, j) = Corr1*_fError[j]*_fError[index];
    }
  }

  //Now make positive definite
  MakePosDef(_fCovMatrix);
  SetCovMatrix(_fCovMatrix);
}

// ********************************************
// ETA An init function for the YAML constructor
// All you really need from the YAML file is the number of Systematics
void ParameterHandlerGeneric::InitialiseFromConfig(const std::vector<std::string>& YAMLFile) {
// ********************************************
  std::map<std::pair<int, int>, std::unique_ptr<TMatrixDSym>> ThrowSubMatrixOverrides;
  LoadAndMergeYAML(YAMLFile, ThrowSubMatrixOverrides);

  const int nThreads = M3::GetNThreads();
  //KS: set Random numbers for each thread so each thread has different seed
  //or for one thread if without MULTITHREAD
  random_number.reserve(nThreads);
  for (int iThread = 0; iThread < nThreads; iThread++) {
    random_number.emplace_back(std::make_unique<TRandom3>(0));
  }
  PrintLength = 35;

  // Set the covariance matrix
  _fNumPar = int(_fYAMLDoc["Systematics"].size());

  ReserveMemory(_fNumPar);

  int i = 0;
  std::vector<std::map<std::string,double>> Correlations(_fNumPar);
  std::map<std::string, int> CorrNamesMap;

  //ETA - read in the systematics. Would be good to add in some checks to make sure
  //that there are the correct number of entries i.e. are the _fNumPar for Names,
  //PreFitValues etc etc.
  for (auto const &param : _fYAMLDoc["Systematics"])
  {
    _fFancyNames[i] = Get<std::string>(param["Systematic"]["Names"]["FancyName"], __FILE__ , __LINE__);
    _fPreFitValue[i] = Get<double>(param["Systematic"]["ParameterValues"]["PreFitValue"], __FILE__ , __LINE__);
    _fIndivStepScale[i] = Get<double>(param["Systematic"]["StepScale"]["MCMC"], __FILE__ , __LINE__);
    _fError[i] = Get<double>(param["Systematic"]["Error"], __FILE__ , __LINE__);
    if(_fError[i] <= 0) {
      MACH3LOG_ERROR("Error for param {}({}) is negative and equal to {}", _fFancyNames[i], i, _fError[i]);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    //ETA - a bit of a fudge but works
    auto TempBoundsVec = GetBounds(param["Systematic"]["ParameterBounds"]);
    _fLowBound[i] = TempBoundsVec[0];
    _fUpBound[i] = TempBoundsVec[1];

    //ETA - now for parameters which are optional and have default values
    _fFlatPrior[i] = GetFromManager<bool>(param["Systematic"]["FlatPrior"], false, __FILE__ , __LINE__);

    // Allow to fix param, this setting should be used only for params which are permanently fixed like baseline, please use global config for fixing param more flexibly
    if(GetFromManager<bool>(param["Systematic"]["FixParam"], false, __FILE__ , __LINE__)) {
      SetFixParameter(_fFancyNames[i]);
    }

    if(param["Systematic"]["SpecialProposal"]) {
      EnableSpecialProposal(param["Systematic"]["SpecialProposal"], i);
    }

    //Fill the map to get the correlations later as well
    CorrNamesMap[param["Systematic"]["Names"]["FancyName"].as<std::string>()]=i;

    //Also loop through the correlations
    if(param["Systematic"]["Correlations"]) {
      for(unsigned int Corr_i = 0; Corr_i < param["Systematic"]["Correlations"].size(); ++Corr_i){
        for (YAML::const_iterator it = param["Systematic"]["Correlations"][Corr_i].begin(); it!=param["Systematic"]["Correlations"][Corr_i].end();++it) {
          Correlations[i][it->first.as<std::string>()] = it->second.as<double>();
        }
      }
    }
    i++;
  } // end loop over para
  if(i != _fNumPar) {
    MACH3LOG_CRITICAL("Inconsistent number of params in Yaml  {} vs {}, this indicate wrong syntax", i, i, _fNumPar);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // load correlation
  LoadCorrelationFromConfig(Correlations, CorrNamesMap);
  if (_fNumPar <= 0) {
    MACH3LOG_ERROR("ParameterHandler object has {} systematics!", _fNumPar);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  for(auto const & matovr : ThrowSubMatrixOverrides){
    SetSubThrowMatrix(matovr.first.first, matovr.first.second, *matovr.second);
  }

  Tunes = std::make_unique<ParameterTunes>(_fYAMLDoc["Systematics"]);

  MACH3LOG_INFO("Created covariance matrix from files: ");
  for(const auto &file : YAMLFile){
    MACH3LOG_INFO("{} ", file);
  }
  MACH3LOG_INFO("----------------");
  MACH3LOG_INFO("Found {} systematics parameters in total", _fNumPar);
  MACH3LOG_INFO("----------------");
}

// ********************************************
void ParameterHandlerGeneric::InitParametersTypeFromConfig() {
// ********************************************
  _fSystToGlobalSystIndexMap.resize(SystType::kSystTypes);

  _fParamType = std::vector<SystType>(_fNumPar);
  _ParameterGroup = std::vector<std::string>(_fNumPar);
  _fSampleNames = std::vector<std::vector<std::string>>(_fNumPar);

  //KS: We know at most how params we expect so reserve memory for max possible params. Later we will shrink to size to not waste memory. Reserving means slightly faster loading and possible less memory fragmentation.
  NormParams.reserve(_fNumPar);
  SplineParams.reserve(_fNumPar);
  FuncParams.reserve(_fNumPar);
  OscParams.reserve(_fNumPar);

  int i = 0;
  unsigned int ParamCounter[SystType::kSystTypes] = {0};
  //ETA - read in the systematics. Would be good to add in some checks to make sure
  //that there are the correct number of entries i.e. are the _fNumPars for Names,
  //PreFitValues etc etc.
  for (auto const &param : _fYAMLDoc["Systematics"])
  {
    _ParameterGroup[i] = Get<std::string>(param["Systematic"]["ParameterGroup"], __FILE__ , __LINE__);
    _fSampleNames[i] = GetFromManager<std::vector<std::string>>(param["Systematic"]["SampleNames"], {}, __FILE__, __LINE__);

    //Fill the map to get the correlations later as well
    auto ParamType = Get<std::string>(param["Systematic"]["Type"], __FILE__ , __LINE__);
    //Now load in variables for spline systematics only
    if (ParamType == SystType_ToString(SystType::kSpline))
    {
      //Set param type
      _fParamType[i] = SystType::kSpline;
      // Fill Spline info
      SplineParams.push_back(GetSplineParameter(param["Systematic"], i));

      //Insert the mapping from the spline index i.e. the length of _fSplineNames etc
      //to the Systematic index i.e. the counter for things like _fSampleID
      _fSystToGlobalSystIndexMap[SystType::kSpline].insert(std::make_pair(ParamCounter[SystType::kSpline], i));
      ParamCounter[SystType::kSpline]++;
    } else if(ParamType == SystType_ToString(SystType::kNorm)) {
      _fParamType[i] = SystType::kNorm;
      NormParams.push_back(GetNormParameter(param["Systematic"], i));
      _fSystToGlobalSystIndexMap[SystType::kNorm].insert(std::make_pair(ParamCounter[SystType::kNorm], i));
      ParamCounter[SystType::kNorm]++;
    } else if(ParamType == SystType_ToString(SystType::kFunc)){
      _fParamType[i] = SystType::kFunc;
      FuncParams.push_back(GetFunctionalParameters(param["Systematic"], i));
      _fSystToGlobalSystIndexMap[SystType::kFunc].insert(std::make_pair(ParamCounter[SystType::kFunc], i));
      ParamCounter[SystType::kFunc]++;
    } else if(ParamType == SystType_ToString(SystType::kOsc)){
      _fParamType[i] = SystType::kOsc;
      OscParams.push_back(GetOscillationParameters(param["Systematic"], i));
      _fSystToGlobalSystIndexMap[SystType::kOsc].insert(std::make_pair(ParamCounter[SystType::kOsc], i));
      ParamCounter[SystType::kOsc]++;
    } else{
      MACH3LOG_ERROR("Given unrecognised systematic type: {}", ParamType);
      std::string expectedTypes = "Expecting ";
      for (int s = 0; s < SystType::kSystTypes; ++s) {
        if (s > 0) expectedTypes += ", ";
        expectedTypes += SystType_ToString(static_cast<SystType>(s)) + "\"";
      }
      expectedTypes += ".";
      MACH3LOG_ERROR(expectedTypes);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    i++;
  } //end loop over params

  //KS We resized them above to all params to fight memory fragmentation, now let's resize to fit only allocated memory to save RAM
  NormParams.shrink_to_fit();
  SplineParams.shrink_to_fit();
  FuncParams.shrink_to_fit();
  OscParams.shrink_to_fit();
}

// ********************************************
ParameterHandlerGeneric::~ParameterHandlerGeneric() {
// ********************************************
  MACH3LOG_DEBUG("Deleting ParameterHandler");
}

// ********************************************
// DB Grab the Spline Names for the relevant SampleName
const std::vector<std::string> ParameterHandlerGeneric::GetSplineParsNamesFromSampleName(const std::string& SampleName) const {
// ********************************************
  std::vector<std::string> returnVec;
  for (auto &pair : _fSystToGlobalSystIndexMap[SystType::kSpline]) {
    auto &SplineIndex = pair.first;
    auto &SystIndex = pair.second;
    if (AppliesToSample(SystIndex, SampleName)) { //If parameter applies to required Sample
      returnVec.push_back(SplineParams[SplineIndex]._fSplineNames);
    }
  }
  return returnVec;
}

// ********************************************
const std::vector<SplineInterpolation> ParameterHandlerGeneric::GetSplineInterpolationFromSampleName(const std::string& SampleName) const {
// ********************************************
  std::vector<SplineInterpolation> returnVec;
  for (auto &pair : _fSystToGlobalSystIndexMap[SystType::kSpline]) {
    auto &SplineIndex = pair.first;
    auto &SystIndex = pair.second;

    if (AppliesToSample(SystIndex, SampleName)) { //If parameter applies to required SampleID
      returnVec.push_back(SplineParams.at(SplineIndex)._SplineInterpolationType);
    }
  }
  return returnVec;
}

// ********************************************
// DB Grab the Spline Modes for the relevant SampleName
const std::vector< std::vector<int> > ParameterHandlerGeneric::GetSplineModeVecFromSampleName(const std::string& SampleName) const {
// ********************************************
  std::vector< std::vector<int> > returnVec;
  //Need a counter or something to correctly get the index in _fSplineModes since it's not of length nPars
  //Should probably just make a std::map<std::string, int> for param name to FD spline index
  for (auto &pair : _fSystToGlobalSystIndexMap[SystType::kSpline]) {
    auto &SplineIndex = pair.first;
    auto &SystIndex = pair.second;
    if (AppliesToSample(SystIndex, SampleName)) { //If parameter applies to required SampleID
      returnVec.push_back(SplineParams.at(SplineIndex)._fSplineModes);
    }
  }
  return returnVec;
}

// ********************************************
// Get Norm params
NormParameter ParameterHandlerGeneric::GetNormParameter(const YAML::Node& param, const int Index) const {
// ********************************************
  NormParameter norm;

  GetBaseParameter(param, Index, norm);

  // ETA size 0 to mean apply to all
  // Ultimately all this information ends up in the NormParams vector
  norm.modes = GetFromManager<std::vector<int>>(param["Mode"], {}, __FILE__ , __LINE__);
  norm.pdgs = GetFromManager<std::vector<int>>(param["NeutrinoFlavour"], {}, __FILE__ , __LINE__);
  norm.preoscpdgs = GetFromManager<std::vector<int>>(param["NeutrinoFlavourUnosc"], {}, __FILE__ , __LINE__);
  norm.targets = GetFromManager<std::vector<int>>(param["TargetNuclei"], {}, __FILE__ , __LINE__);

  if(_fLowBound[Index] < 0.) {
    MACH3LOG_ERROR("Normalisation Parameter {} ({}), has lower parameters bound which can go below 0 and is equal {}",
                   GetParFancyName(Index), Index, _fLowBound[Index]);
    MACH3LOG_ERROR("Normalisation parameters can't go bellow 0 as this is unphysical");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return norm;
}

// ********************************************
// Get Base Param
void ParameterHandlerGeneric::GetBaseParameter(const YAML::Node& param, const int Index, TypeParameterBase& Parameter) const {
// ********************************************
  Parameter.name = GetParFancyName(Index);

  // Set the global parameter index of the normalisation parameter
  Parameter.index = Index;

  int NumKinematicCuts = 0;
  if(param["KinematicCuts"]) {
    NumKinematicCuts = int(param["KinematicCuts"].size());

    std::vector<std::string> TempKinematicStrings;
    std::vector<std::vector<std::vector<double>>> TempKinematicBounds;
    //First element of TempKinematicBounds is always -999, and size is then 3
    for(int KinVar_i = 0 ; KinVar_i < NumKinematicCuts ; ++KinVar_i) {
      //ETA: This is a bit messy, Kinematic cuts is a list of maps
      for (YAML::const_iterator it = param["KinematicCuts"][KinVar_i].begin();it!=param["KinematicCuts"][KinVar_i].end();++it) {
        TempKinematicStrings.push_back(it->first.as<std::string>());
        TempKinematicBounds.push_back(Get2DBounds(it->second));
      }
      if(TempKinematicStrings.size() == 0) {
        MACH3LOG_ERROR("Received a KinematicCuts node but couldn't read the contents (it's a list of single-element dictionaries (python) = map of pairs (C++))");
        MACH3LOG_ERROR("For Param {}", Parameter.name);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }//KinVar_i
    Parameter.KinematicVarStr = TempKinematicStrings;
    Parameter.Selection = TempKinematicBounds;
  }

  //Next ones are kinematic bounds on where normalisation parameter should apply
  //We set a bool to see if any bounds exist so we can short-circuit checking all of them every step
  bool HasKinBounds = false;

  if(Parameter.KinematicVarStr.size() > 0) HasKinBounds = true;

  Parameter.hasKinBounds = HasKinBounds;
  //End of kinematic bound checking
}

// ********************************************
// Grab the global syst index for the relevant SampleName
// i.e. get a vector of size nSplines where each entry is filled with the global syst number
const std::vector<int> ParameterHandlerGeneric::GetGlobalSystIndexFromSampleName(const std::string& SampleName, const SystType Type) const {
// ********************************************
  std::vector<int> returnVec;
  for (auto &pair : _fSystToGlobalSystIndexMap[Type]) {
    auto &SystIndex = pair.second;
    if (AppliesToSample(SystIndex, SampleName)) { //If parameter applies to required SampleID
      returnVec.push_back(SystIndex);
    }
  }
  return returnVec;
}

// ********************************************
// Grab the global syst index for the relevant SampleName
// i.e. get a vector of size nSplines where each entry is filled with the global syst number
const std::vector<int> ParameterHandlerGeneric::GetSystIndexFromSampleName(const std::string& SampleName,  const SystType Type) const {
// ********************************************
  std::vector<int> returnVec;
  for (auto &pair : _fSystToGlobalSystIndexMap[Type]) {
    auto &SplineIndex = pair.first;
    auto &systIndex = pair.second;
    if (AppliesToSample(systIndex, SampleName)) { //If parameter applies to required SampleID
      returnVec.push_back(SplineIndex);
    }
  }
  return returnVec;
}

// ********************************************
// Get Spline params
SplineParameter ParameterHandlerGeneric::GetSplineParameter(const YAML::Node& param, const int Index) const {
// ********************************************
  SplineParameter Spline;

  GetBaseParameter(param, Index, Spline);
  auto& SplinePar = param["SplineInformation"];
  //Now get the Spline interpolation type
  if (SplinePar["InterpolationType"]) {
    for(int InterpType = 0; InterpType < kSplineInterpolations ; ++InterpType){
      if(SplinePar["InterpolationType"].as<std::string>() == SplineInterpolation_ToString(SplineInterpolation(InterpType)))
        Spline._SplineInterpolationType = SplineInterpolation(InterpType);
    }
  } else { //KS: By default use TSpline3
    Spline._SplineInterpolationType = kTSpline3;
  }

  Spline._fSplineNames = Get<std::string>(SplinePar["SplineName"], __FILE__ , __LINE__);
  Spline._SplineKnotUpBound = GetFromManager<double>(SplinePar["SplineKnotUpBound"], M3::DefSplineKnotUpBound, __FILE__ , __LINE__);
  Spline._SplineKnotLowBound = GetFromManager<double>(SplinePar["SplineKnotLowBound"], M3::DefSplineKnotLowBound, __FILE__ , __LINE__);

  if(Spline._SplineKnotUpBound != M3::DefSplineKnotUpBound ||  Spline._SplineKnotLowBound != M3::DefSplineKnotLowBound) {
    MACH3LOG_WARN("Spline knot capping enabled with bounds [{}, {}]. For reliable fits, consider modifying the input generation instead.",
                  Spline._SplineKnotLowBound, Spline._SplineKnotUpBound);
  }
  //If there is no mode information given then this will be an empty vector
  Spline._fSplineModes = GetFromManager(SplinePar["Mode"], std::vector<int>(), __FILE__ , __LINE__);

  return Spline;
}

// ********************************************
bool ParameterHandlerGeneric::AppliesToSample(const int SystIndex, const std::string& SampleName) const {
// ********************************************
  // Empty means apply to all
  if (_fSampleNames[SystIndex].size() == 0) return true;

  // Check for unsupported wildcards in SampleName
  if (SampleName.find('*') != std::string::npos) {
    MACH3LOG_ERROR("Wildcards ('*') are not supported in sample name: '{}'", SampleName);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return M3::CaseInsensitiveMatchAny(SampleName, _fSampleNames[SystIndex]);
}

// ********************************************
// Get Func params
FunctionalParameter ParameterHandlerGeneric::GetFunctionalParameters(const YAML::Node& param, const int Index) const {
// ********************************************
  FunctionalParameter func;
  GetBaseParameter(param, Index, func);

  func.modes = GetFromManager<std::vector<int>>(param["Mode"], std::vector<int>(), __FILE__ , __LINE__);

  return func;
}

// ********************************************
// Get Osc params
OscillationParameter ParameterHandlerGeneric::GetOscillationParameters(const YAML::Node& param, const int Index) const {
// ********************************************
  OscillationParameter OscParamInfo;
  GetBaseParameter(param, Index, OscParamInfo);

  return OscParamInfo;
}

// ********************************************
// HH: Grab the Functional parameters for the relevant SampleName
const std::vector<FunctionalParameter> ParameterHandlerGeneric::GetFunctionalParametersFromSampleName(const std::string& SampleName) const {
// ********************************************
  return GetTypeParamsFromSampleName(_fSystToGlobalSystIndexMap[SystType::kFunc], FuncParams, SampleName);
}

// ********************************************
// DB Grab the Normalisation parameters for the relevant SampleName
const std::vector<NormParameter> ParameterHandlerGeneric::GetNormParsFromSampleName(const std::string& SampleName) const {
// ********************************************
  return GetTypeParamsFromSampleName(_fSystToGlobalSystIndexMap[SystType::kNorm], NormParams, SampleName);
}

// ********************************************
// KS Grab the Spline parameters for the relevant SampleName
const std::vector<SplineParameter> ParameterHandlerGeneric::GetSplineParsFromSampleName(const std::string& SampleName) const {
// ********************************************
  return GetTypeParamsFromSampleName(_fSystToGlobalSystIndexMap[SystType::kSpline], SplineParams, SampleName);
}

// ********************************************
template<typename ParamT>
std::vector<ParamT> ParameterHandlerGeneric::GetTypeParamsFromSampleName(const std::map<int, int>& indexMap,
                                                                         const std::vector<ParamT>& params, const std::string& SampleName) const {
// ********************************************
  std::vector<ParamT> returnVec;
  for (const auto& pair : indexMap) {
    const auto& localIndex = pair.first;
    const auto& globalIndex = pair.second;
    if (AppliesToSample(globalIndex, SampleName)) {
      returnVec.push_back(params[localIndex]);
    }
  }
  return returnVec;
}

// ********************************************
// DB Grab the number of parameters for the relevant SampleName
int ParameterHandlerGeneric::GetNumParamsFromSampleName(const std::string& SampleName, const SystType Type) const {
// ********************************************
  int returnVal = 0;
  IterateOverParams(SampleName,
    [&](int i) { return GetParamType(i) == Type; }, // Filter condition
    [&](int) { returnVal += 1; } // Action to perform if filter passes
  );
  return returnVal;
}

// ********************************************
// DB Grab the parameter names for the relevant SampleName
const std::vector<std::string> ParameterHandlerGeneric::GetParsNamesFromSampleName(const std::string& SampleName, const SystType Type) const {
// ********************************************
  std::vector<std::string> returnVec;
  IterateOverParams(SampleName,
    [&](int i) { return GetParamType(i) == Type; }, // Filter condition
    [&](int i) { returnVec.push_back(GetParFancyName(i)); } // Action to perform if filter passes
  );
  return returnVec;
}

// ********************************************
// DB DB Grab the parameter indices for the relevant SampleName
const std::vector<int> ParameterHandlerGeneric::GetParsIndexFromSampleName(const std::string& SampleName, const SystType Type) const {
// ********************************************
  std::vector<int> returnVec;
  IterateOverParams(SampleName,
    [&](int i) { return GetParamType(i) == Type; }, // Filter condition
    [&](int i) { returnVec.push_back(i); } // Action to perform if filter passes
  );
  return returnVec;
}

// ********************************************
template <typename FilterFunc, typename ActionFunc>
void ParameterHandlerGeneric::IterateOverParams(const std::string& SampleName, FilterFunc filter, ActionFunc action) const {
// ********************************************
  for (int i = 0; i < _fNumPar; ++i) {
    if ((AppliesToSample(i, SampleName)) && filter(i)) { // Common filter logic
      action(i); // Specific action for each function
    }
  }
}

// ********************************************
void ParameterHandlerGeneric::InitParameters() {
// ********************************************
  for (int i = 0; i < _fNumPar; ++i) {
    //ETA - set the name to be param_% as this is what ProcessorMCMC expects
    _fNames[i] = "param_"+std::to_string(i);

    // KS: Plenty of MaCh3 processing script rely on osc params having "fancy name" this is to maintain backward compatibility with them
    if(_fParamType[i] == kOsc) {
      _fNames[i] = _fFancyNames[i];

      if(_ParameterGroup[i] != "Osc"){
        MACH3LOG_ERROR("Parameter {}, is of type Oscillation but doesn't belong to Osc group", _fFancyNames[i]);
        MACH3LOG_ERROR("It belongs to {} group", _ParameterGroup[i]);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
    }
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wuseless-cast"
    // Set ParameterHandler parameters (Curr = current, Prop = proposed, Sigma = step)
    _fCurrVal[i] = _fPreFitValue[i];
    _fPropVal[i] = static_cast<M3::float_t>(_fCurrVal[i]);
    #pragma GCC diagnostic pop
  }
  Randomize();
  //KS: Transfer the starting parameters to the PCA basis, you don't want to start with zero..
  if (pca) {
    PCAObj->SetInitialParameters();
  }
}

// ********************************************
// Print everything we know about the inputs we're Getting
void ParameterHandlerGeneric::Print() const {
// ********************************************
  MACH3LOG_INFO("#################################################");
  MACH3LOG_INFO("Printing ParameterHandlerGeneric:");

  PrintGlobablInfo();

  PrintNormParams();

  PrintSplineParams();

  PrintFunctionalParams();

  PrintOscillationParams();

  PrintParameterGroups();

  MACH3LOG_INFO("Finished");
  MACH3LOG_INFO("#################################################");

  CheckCorrectInitialisation();
} // End

// ********************************************
void ParameterHandlerGeneric::PrintGlobablInfo() const {
// ********************************************
  MACH3LOG_INFO("============================================================================================================================================================");
  MACH3LOG_INFO("{:<5} {:2} {:<40} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<20} {:2} {:<10}", "#", "|", "Name", "|", "Prior", "|", "Error", "|", "Lower", "|", "Upper", "|", "StepScale", "|", "SampleNames", "|", "Type");
  MACH3LOG_INFO("------------------------------------------------------------------------------------------------------------------------------------------------------------");
  for (int i = 0; i < GetNumParams(); i++) {
    std::string ErrString = fmt::format("{:.2f}", _fError[i]);
    std::string SampleNameString = "";
    for (const auto& SampleName : _fSampleNames[i]) {
      if (!SampleNameString.empty()) {
        SampleNameString += ", ";
      }
      SampleNameString += SampleName;
    }
    MACH3LOG_INFO("{:<5} {:2} {:<40} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<20} {:2} {:<10}", i, "|", GetParFancyName(i), "|", _fPreFitValue[i], "|", "+/- " + ErrString, "|", _fLowBound[i], "|", _fUpBound[i], "|", _fIndivStepScale[i], "|", SampleNameString, "|", SystType_ToString(_fParamType[i]));
  }
  MACH3LOG_INFO("============================================================================================================================================================");
}

// ********************************************
void ParameterHandlerGeneric::PrintNormParams() const {
// ********************************************
  // Output the normalisation parameters as a sanity check!
  MACH3LOG_INFO("Normalisation parameters:  {}", NormParams.size());
  if(_fSystToGlobalSystIndexMap[SystType::kNorm].size() == 0) return;

  bool have_parameter_with_kin_bounds = false;

  //KS: Consider making some class producing table..
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┬────────────────────┬────────────────────┬────────────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│{3:20}│{4:20}│{5:20}│", "#", "Global #", "Name", "Int. mode", "Target", "pdg");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┼────────────────────┼────────────────────┼────────────────────┤");

  for (unsigned int i = 0; i < NormParams.size(); ++i)
  {
    std::string intModeString;
    for (unsigned int j = 0; j < NormParams[i].modes.size(); j++) {
      intModeString += std::to_string(NormParams[i].modes[j]);
      intModeString += " ";
    }
    if (NormParams[i].modes.empty()) intModeString += "all";

    std::string targetString;
    for (unsigned int j = 0; j < NormParams[i].targets.size(); j++) {
      targetString += std::to_string(NormParams[i].targets[j]);
      targetString += " ";
    }
    if (NormParams[i].targets.empty()) targetString += "all";

    std::string pdgString;
    for (unsigned int j = 0; j < NormParams[i].pdgs.size(); j++) {
      pdgString += std::to_string(NormParams[i].pdgs[j]);
      pdgString += " ";
    }
    if (NormParams[i].pdgs.empty()) pdgString += "all";

    MACH3LOG_INFO("│{: <4}│{: <10}│{: <40}│{: <20}│{: <20}│{: <20}│", i, NormParams[i].index, NormParams[i].name, intModeString, targetString, pdgString);

    if(NormParams[i].hasKinBounds) have_parameter_with_kin_bounds = true;
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┴────────────────────┴────────────────────┴────────────────────┘");

  if(have_parameter_with_kin_bounds) {
    MACH3LOG_INFO("Normalisation parameters KinematicCuts information");
    MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┬────────────────────┬────────────────────────────────────────┐");
    MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│{3:20}│{4:40}│", "#", "Global #", "Name", "KinematicCut", "Value");
    MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┼────────────────────┼────────────────────────────────────────┤");
    for (unsigned int i = 0; i < NormParams.size(); ++i)
    {
      //skip parameters with no KinematicCuts
      if(!NormParams[i].hasKinBounds) continue;

      const long unsigned int ncuts = NormParams[i].KinematicVarStr.size();
      for(long unsigned int icut = 0; icut < ncuts; icut++) {
        std::string kinematicCutValueString;
        for(const auto & value : NormParams[i].Selection[icut]) {
          for (const auto& v : value) {
            kinematicCutValueString += fmt::format("{:.2f} ", v);
          }
        }
        if(icut == 0)
          MACH3LOG_INFO("│{: <4}│{: <10}│{: <40}│{: <20}│{: <40}│", i, NormParams[i].index, NormParams[i].name, NormParams[i].KinematicVarStr[icut], kinematicCutValueString);
        else
          MACH3LOG_INFO("│{: <4}│{: <10}│{: <40}│{: <20}│{: <40}│", "", "", "", NormParams[i].KinematicVarStr[icut], kinematicCutValueString);
      }//icut
    }//i
    MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┴────────────────────┴────────────────────────────────────────┘");
  }
  else
    MACH3LOG_INFO("No normalisation parameters have KinematicCuts defined");
}

// ********************************************
void ParameterHandlerGeneric::PrintSplineParams() const {
// ********************************************
  MACH3LOG_INFO("Spline parameters: {}", _fSystToGlobalSystIndexMap[SystType::kSpline].size());
  if(_fSystToGlobalSystIndexMap[SystType::kSpline].size() == 0) return;
  MACH3LOG_INFO("=====================================================================================================================================================================");
  MACH3LOG_INFO("{:<4} {:<2} {:<40} {:<2} {:<40} {:<2} {:<20} {:<2} {:<20} {:<2} {:<20} {:<2}", "#", "|", "Name", "|", "Spline Name", "|", "Spline Interpolation", "|", "Low Knot Bound", "|", "Up Knot Bound", "|");
  MACH3LOG_INFO("---------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  for (auto &pair : _fSystToGlobalSystIndexMap[SystType::kSpline]) {
    auto &SplineIndex = pair.first;
    auto &GlobalIndex = pair.second;

    MACH3LOG_INFO("{:<4} {:<2} {:<40} {:<2} {:<40} {:<2} {:<20} {:<2} {:<20} {:<2} {:<20} {:<2}",
                  SplineIndex, "|", GetParFancyName(GlobalIndex), "|",
                  SplineParams[SplineIndex]._fSplineNames, "|",
                  SplineInterpolation_ToString(GetParSplineInterpolation(SplineIndex)), "|",
                  GetParSplineKnotLowerBound(SplineIndex), "|",
                  GetParSplineKnotUpperBound(SplineIndex), "|");
  }
  MACH3LOG_INFO("=====================================================================================================================================================================");
}

// ********************************************
void ParameterHandlerGeneric::PrintFunctionalParams() const {
// ********************************************
  MACH3LOG_INFO("Functional parameters: {}", _fSystToGlobalSystIndexMap[SystType::kFunc].size());
  if(_fSystToGlobalSystIndexMap[SystType::kFunc].size() == 0) return;
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│", "#", "Global #", "Name");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┤");
  for (auto &pair : _fSystToGlobalSystIndexMap[SystType::kFunc]) {
    auto &FuncIndex = pair.first;
    auto &GlobalIndex = pair.second;
    MACH3LOG_INFO("│{0:4}│{1:<10}│{2:40}│", std::to_string(FuncIndex), GlobalIndex, GetParFancyName(GlobalIndex));
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┘");
}

// ********************************************
void ParameterHandlerGeneric::PrintOscillationParams() const {
// ********************************************
  MACH3LOG_INFO("Oscillation parameters: {}", _fSystToGlobalSystIndexMap[SystType::kOsc].size());
  if(_fSystToGlobalSystIndexMap[SystType::kOsc].size() == 0) return;
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│", "#", "Global #", "Name");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┤");
  for (auto &pair : _fSystToGlobalSystIndexMap[SystType::kOsc]) {
    auto &OscIndex = pair.first;
    auto &GlobalIndex = pair.second;
    MACH3LOG_INFO("│{0:4}│{1:<10}│{2:40}│", std::to_string(OscIndex), GlobalIndex, GetParFancyName(GlobalIndex));
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┘");
}

// ********************************************
void ParameterHandlerGeneric::PrintParameterGroups() const {
// ********************************************
  // KS: Create a map to store the counts of unique strings, in principle this could be in header file
  std::unordered_map<std::string, int> paramCounts;

  std::for_each(_ParameterGroup.begin(), _ParameterGroup.end(),
                [&paramCounts](const std::string& param) {
                  paramCounts[param]++;
                });

  MACH3LOG_INFO("Printing parameter groups");
  // Output the counts
  for (const auto& pair : paramCounts) {
    MACH3LOG_INFO("Found {}: {} params", pair.second, pair.first);
  }
}

// ********************************************
std::vector<std::string> ParameterHandlerGeneric::GetUniqueParameterGroups() const {
// ********************************************
  std::unordered_set<std::string> uniqueGroups;

  // Fill the set with unique values
  for (const auto& param : _ParameterGroup) {
    uniqueGroups.insert(param);
  }

  // Convert to vector and return
  std::vector<std::string> result(uniqueGroups.begin(), uniqueGroups.end());
  return result;
}

// ********************************************
// KS: Check if matrix is correctly initialised
void ParameterHandlerGeneric::CheckCorrectInitialisation() const {
// ********************************************
  // KS: Lambda Function which simply checks if there are no duplicates in std::vector
  auto CheckForDuplicates = [](const std::vector<std::string>& names, const std::string& nameType, bool Warning) {
    std::unordered_map<std::string, size_t> seenStrings;
    for (size_t i = 0; i < names.size(); ++i) {
      const auto& name = names[i];
      if (seenStrings.find(name) != seenStrings.end()) {
        size_t firstIndex = seenStrings[name];
        if(Warning){
          MACH3LOG_WARN("There are two systematics with the same {} '{}', first at index {}, and again at index {}", nameType, name, firstIndex, i);
          MACH3LOG_WARN("Is this desired?");
          return;
        } else {
          MACH3LOG_CRITICAL("There are two systematics with the same {} '{}', first at index {}, and again at index {}", nameType, name, firstIndex, i);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
      seenStrings[name] = i;
    }
  };
  std::vector<std::string> SplineNameTemp(SplineParams.size());
  for(size_t it = 0; it < SplineParams.size(); it++){
    SplineNameTemp[it] = SplineParams[it]._fSplineNames;
  }
  // KS: Checks if there are no duplicates in fancy names etc, this can happen if we merge configs etc
  CheckForDuplicates(_fFancyNames, "_fFancyNames", false);
  CheckForDuplicates(SplineNameTemp, "_fSplineNames", true);
}

// ********************************************
// Function to set to prior parameters of a given group
void ParameterHandlerGeneric::SetGroupOnlyParameters(const std::vector< std::string>& Groups) {
// ********************************************
  for(size_t i = 0; i < Groups.size(); i++){
    SetGroupOnlyParameters(Groups[i]);
  }
}

// ********************************************
// Function to set to prior parameters of a given group
void ParameterHandlerGeneric::SetGroupOnlyParameters(const std::string& Group, const std::vector<double>& Pars) {
// ********************************************
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuseless-cast"
  // If empty, set the proposed to prior
  if (Pars.empty()) {
    for (int i = 0; i < _fNumPar; i++) {
      if(IsParFromGroup(i, Group)) _fPropVal[i] = static_cast<M3::float_t>(_fPreFitValue[i]);
    }
  } else{
    const size_t ExpectedSize = static_cast<size_t>(GetNumParFromGroup(Group));
    if (Pars.size() != ExpectedSize) {
      MACH3LOG_ERROR("Number of param in group {} is {}, while you passed {}", Group, ExpectedSize, Pars.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    int Counter = 0;
    for (int i = 0; i < _fNumPar; i++) {
      // If belongs to group set value from parsed vector, otherwise use propose value
      if(IsParFromGroup(i, Group)){
        _fPropVal[i] = static_cast<M3::float_t>(Pars[Counter]);
        Counter++;
      }
    }
  }
  // And if pca make the transfer
  if (pca) {
    PCAObj->TransferToPCA();
    PCAObj->TransferToParam();
  }
  #pragma GCC diagnostic pop
}

// ********************************************
// Set parameters to be fixed in a given group
void ParameterHandlerGeneric::SetFixGroupOnlyParameters(const std::string& Group) {
// ********************************************
  for (int i = 0; i < _fNumPar; ++i)
    if(IsParFromGroup(i, Group)) SetFixParameter(i);
}

// ********************************************
// Set parameters of several groups to be fixed
void ParameterHandlerGeneric::SetFixGroupOnlyParameters(const std::vector< std::string>& Groups) {
// ********************************************
  for(size_t i = 0; i < Groups.size(); i++)
    SetFixGroupOnlyParameters(Groups[i]);
}

// ********************************************
// Set parameters to be free in a given group
void ParameterHandlerGeneric::SetFreeGroupOnlyParameters(const std::string& Group) {
// ********************************************
  for (int i = 0; i < _fNumPar; ++i)
    if(IsParFromGroup(i, Group)) SetFreeParameter(i);
}

// ********************************************
// Set parameters of several groups to be fixed
void ParameterHandlerGeneric::SetFreeGroupOnlyParameters(const std::vector< std::string>& Groups) {
// ********************************************
  for(size_t i = 0; i < Groups.size(); i++)
    SetFreeGroupOnlyParameters(Groups[i]);
}

// ********************************************
// Checks if parameter belongs to a given group
bool ParameterHandlerGeneric::IsParFromGroup(const int i, const std::string& Group) const {
// ********************************************
  std::string groupLower = Group;
  std::string paramGroupLower = _ParameterGroup[i];

  // KS: Convert both strings to lowercase, this way comparison will be case insensitive
  std::transform(groupLower.begin(), groupLower.end(), groupLower.begin(), ::tolower);
  std::transform(paramGroupLower.begin(), paramGroupLower.end(), paramGroupLower.begin(), ::tolower);

  return groupLower == paramGroupLower;
}

// ********************************************
int ParameterHandlerGeneric::GetNumParFromGroup(const std::string& Group) const {
// ********************************************
  int Counter = 0;
  for (int i = 0; i < _fNumPar; i++) {
    if(IsParFromGroup(i, Group)) Counter++;
  }
  return Counter;
}

// ********************************************
// DB Grab the Normalisation parameters for the relevant sample name
std::vector<const M3::float_t*> ParameterHandlerGeneric::GetOscParsFromSampleName(const std::string& SampleName) const {
// ********************************************
  std::vector<const M3::float_t*> returnVec;
  for (const auto& pair : _fSystToGlobalSystIndexMap[SystType::kOsc]) {
    const auto& globalIndex = pair.second;
    if (AppliesToSample(globalIndex, SampleName)) {
      returnVec.push_back(RetPointer(globalIndex));
    }
  }
  return returnVec;
}

// ********************************************
// Dump Matrix to ROOT file, useful when we need to pass matrix info to another fitting group
void ParameterHandlerGeneric::DumpMatrixToFile(const std::string& Name) {
// ********************************************
  // KS: Ideally we remove it eventually...
  TH2D* CorrMatrix = GetCorrelationMatrix();
  M3::DumpParamHandlerToFile(_fNumPar,
                             _fPreFitValue,
                             _fError,
                             _fLowBound,
                             _fUpBound,
                             _fIndivStepScale,
                             _fFancyNames,
                             _fFlatPrior,
                             SplineParams,
                             covMatrix,
                             CorrMatrix,
                             Name);
  delete CorrMatrix;
  MACH3LOG_INFO("Finished dumping ParameterHandler object");
}
