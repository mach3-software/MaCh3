#include "covarianceXsec.h"

// ********************************************
// ETA - YAML constructor
// this will replace the root file constructor but let's keep it in
// to do some validations
covarianceXsec::covarianceXsec(const std::vector<std::string>& YAMLFile, const char *name, double threshold, int FirstPCA, int LastPCA)
               : covarianceBase(YAMLFile, name, threshold, FirstPCA, LastPCA){
// ********************************************

  InitXsecFromConfig();

  //ETA - again this really doesn't need to be hear...
  for (int i = 0; i < _fNumPar; i++)
  {
    // Sort out the print length
    if(_fNames[i].length() > PrintLength) PrintLength = _fNames[i].length();
  } // end the for loop

  MACH3LOG_DEBUG("Constructing instance of covarianceXsec");
  initParams();
  // Print
  Print();
}

// ********************************************
void covarianceXsec::InitXsecFromConfig() {
// ********************************************
  _fSystToGlobalSystIndexMap.resize(kSystTypes);

  _fDetID = std::vector<int>(_fNumPar);
  _fParamType = std::vector<SystType>(_fNumPar);
  _ParameterGroup = std::vector<std::string>(_fNumPar);

  //KS: We know at most how params we expect so reserve memory for max possible params. Later we will shrink to size to not waste memory. Reserving means slightly faster loading and possible less memory fragmentation.
  NormParams.reserve(_fNumPar);
  SplineParams.reserve(_fNumPar);

  int i = 0;
  unsigned int ParamCounter[kSystTypes] = {0};
  //ETA - read in the systematics. Would be good to add in some checks to make sure
  //that there are the correct number of entries i.e. are the _fNumPars for Names,
  //PreFitValues etc etc.
  for (auto const &param : _fYAMLDoc["Systematics"])
  {
    _fDetID[i] = (param["Systematic"]["DetID"].as<int>());
    _ParameterGroup[i] = (param["Systematic"]["ParameterGroup"].as<std::string>());

    //Fill the map to get the correlations later as well
    std::string ParamType = param["Systematic"]["Type"].as<std::string>();
    //Now load in variables for spline systematics only
    if (ParamType.find(SystType_ToString(SystType::kSpline)) != std::string::npos)
    {
      //Set param type
      _fParamType[i] = SystType::kSpline;
      // Fill Spline info
      SplineParams.push_back(GetXsecSpline(param["Systematic"]));

      if (param["Systematic"]["SplineInformation"]["SplineName"]) {
        _fSplineNames.push_back(param["Systematic"]["SplineInformation"]["SplineName"].as<std::string>());
      }

      //Insert the mapping from the spline index i.e. the length of _fSplineNames etc
      //to the Systematic index i.e. the counter for things like _fDetID and _fDetID
      _fSystToGlobalSystIndexMap[kSpline].insert(std::make_pair(ParamCounter[kSpline], i));
      ParamCounter[kSpline]++;
    } else if(param["Systematic"]["Type"].as<std::string>() == SystType_ToString(SystType::kNorm)) {
      _fParamType[i] = SystType::kNorm;
      NormParams.push_back(GetXsecNorm(param["Systematic"], i));
      _fSystToGlobalSystIndexMap[kNorm].insert(std::make_pair(ParamCounter[kNorm], i));
      ParamCounter[kNorm]++;
    }
    else if(param["Systematic"]["Type"].as<std::string>() == SystType_ToString(SystType::kFunc)){
      _fParamType[i] = SystType::kFunc;
      _fSystToGlobalSystIndexMap[kFunc].insert(std::make_pair(ParamCounter[kFunc], i));
      ParamCounter[kFunc]++;
    }
    else{
      MACH3LOG_ERROR("Given unrecognised systematic type: {}", param["Systematic"]["Type"].as<std::string>());
      std::string expectedTypes = "Expecting ";
      for (int s = 0; s < kSystTypes; ++s) {
        if (s > 0) expectedTypes += ", ";
        expectedTypes += SystType_ToString(static_cast<SystType>(s)) + "\"";
      }
      expectedTypes += ".";
      MACH3LOG_ERROR(expectedTypes);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    i++;
  }

  //Add a sanity check,
  if(_fSplineNames.size() != ParamCounter[kSpline]){
    MACH3LOG_ERROR("_fSplineNames is of size {} but found {} spline parameters", _fSplineNames.size(), ParamCounter[kSpline]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  //KS We resized them above to all params to fight memory fragmentation, now let's resize to fit only allocated memory to save RAM
  NormParams.shrink_to_fit();
  SplineParams.shrink_to_fit();

  return;
}

// ********************************************
covarianceXsec::~covarianceXsec() {
// ********************************************

}

// ********************************************
// DB Grab the Spline Names for the relevant DetID
const std::vector<std::string> covarianceXsec::GetSplineParsNamesFromDetID(const int DetID) {
// ********************************************
  std::vector<std::string> returnVec;
  for (auto &pair : _fSystToGlobalSystIndexMap[kSpline]) {
    auto &SplineIndex = pair.first;
    auto &SystIndex = pair.second;
    if ((GetParDetID(SystIndex) & DetID )){
      returnVec.push_back(_fSplineNames.at(SplineIndex));
    }
  }
  return returnVec;
}

// ********************************************
// DB Grab the Spline Modes for the relevant DetID
const std::vector< std::vector<int> > covarianceXsec::GetSplineModeVecFromDetID(const int DetID) {
// ********************************************
  std::vector< std::vector<int> > returnVec;
  //Need a counter or something to correctly get the index in _fSplineModes since it's not of length nPars
  //Should probably just make a std::map<std::string, int> for param name to FD spline index
  for (auto &pair : _fSystToGlobalSystIndexMap[kSpline]) {
    auto &SplineIndex = pair.first;
    auto &SystIndex = pair.second;
    if ((GetParDetID(SystIndex) & DetID)) { //If parameter applies to required DetID
      returnVec.push_back(SplineParams.at(SplineIndex)._fSplineModes);
    }
  }
  return returnVec;
}

// ********************************************
// Get Norm params
XsecNorms4 covarianceXsec::GetXsecNorm(const YAML::Node& param, const int Index) {
// ********************************************
  XsecNorms4 norm;
  norm.name = GetParFancyName(Index);

  // ETA Empty DummyVector can be used to specify no cut for mode, target and neutrino flavour
  // ETA Has to be of size 0 to mean apply to all
  std::vector<int> DummyModeVec;
  //Ultimately all this information ends up in the NormParams vector

  //Copy the mode information into an XsecNorms4 struct
  norm.modes = GetFromManager<std::vector<int>>(param["Mode"], DummyModeVec);
  norm.pdgs = GetFromManager<std::vector<int>>(param["NeutrinoFlavour"], DummyModeVec);
  norm.preoscpdgs = GetFromManager<std::vector<int>>(param["NeutrinoFlavourUnosc"], DummyModeVec);
  norm.targets = GetFromManager<std::vector<int>>(param["TargetNuclei"], DummyModeVec);

  //ETA - I think this can go in the norm parameters only if statement above
  int NumKinematicCuts = 0;
  if(param["KinematicCuts"]){

    NumKinematicCuts = param["KinematicCuts"].size();

    std::vector<std::string> TempKinematicStrings;
    std::vector<std::vector<double>> TempKinematicBounds;
    //First element of TempKinematicBounds is always -999, and size is then 3
    for(int KinVar_i = 0 ; KinVar_i < NumKinematicCuts ; ++KinVar_i){
      //ETA
      //This is a bit messy, Kinematic cuts is a list of maps
      for (YAML::const_iterator it = param["KinematicCuts"][KinVar_i].begin();it!=param["KinematicCuts"][KinVar_i].end();++it) {
        TempKinematicStrings.push_back(it->first.as<std::string>());
        TempKinematicBounds.push_back(it->second.as<std::vector<double>>());
        std::vector<double> bounds = it->second.as<std::vector<double>>();
      }
    }
    norm.KinematicVarStr = TempKinematicStrings;
    norm.Selection = TempKinematicBounds;
  }

  //Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
  //We set a bool to see if any bounds exist so we can short-circuit checking all of them every step
  bool HasKinBounds = false;

  if(norm.KinematicVarStr.size() > 0) HasKinBounds = true;

  norm.hasKinBounds = HasKinBounds;
  //End of kinematic bound checking

  // Set the global parameter index of the normalisation parameter
  norm.index = Index;

  return norm;
}

// ********************************************
// Grab the global syst index for the relevant DetID
// i.e. get a vector of size nSplines where each entry is filled with the global syst number
const std::vector<int> covarianceXsec::GetGlobalSystIndexFromDetID(const int DetID, const SystType Type) {
// ********************************************
  std::vector<int> returnVec;
  for (auto &pair : _fSystToGlobalSystIndexMap[Type]) {
    auto &SystIndex = pair.second;
    if ((GetParDetID(SystIndex) & DetID)) { //If parameter applies to required DetID
      returnVec.push_back(SystIndex);
    }
  }
  return returnVec;
}

// ********************************************
// Grab the global syst index for the relevant DetID
// i.e. get a vector of size nSplines where each entry is filled with the global syst number
const std::vector<int> covarianceXsec::GetSystIndexFromDetID(int DetID,  const SystType Type) {
// ********************************************
  std::vector<int> returnVec;
  for (auto &pair : _fSystToGlobalSystIndexMap[Type]) {
    auto &SplineIndex = pair.first;
    auto &SystIndex = pair.second;
    if ((GetParDetID(SystIndex) & DetID)) { //If parameter applies to required DetID
      returnVec.push_back(SplineIndex);
    }
  }
  return returnVec;
}

// ********************************************
// Get Norm params
XsecSplines1 covarianceXsec::GetXsecSpline(const YAML::Node& param) {
// ********************************************
  XsecSplines1 Spline;

  //Now get the Spline interpolation type
  if (param["SplineInformation"]["InterpolationType"]){
    for(int InterpType = 0; InterpType < kSplineInterpolations ; ++InterpType){
      if(param["SplineInformation"]["InterpolationType"].as<std::string>() == SplineInterpolation_ToString(SplineInterpolation(InterpType)))
        Spline._SplineInterpolationType = SplineInterpolation(InterpType);
    }
  } else { //KS: By default use TSpline3
    Spline._SplineInterpolationType = SplineInterpolation(kTSpline3);
  }
  Spline._SplineKnotUpBound = GetFromManager<double>(param["SplineInformation"]["SplineKnotUpBound"], 9999);
  Spline._SplineKnotLowBound = GetFromManager<double>(param["SplineInformation"]["SplineKnotLowBound"], -9999);

  //If there is no mode information given then this will be an empty vector
  Spline._fSplineModes = GetFromManager(param["SplineInformation"]["Mode"], std::vector<int>());

  return Spline;
}

// ********************************************
// DB Grab the Normalisation parameters for the relevant DetID
const std::vector<XsecNorms4> covarianceXsec::GetNormParsFromDetID(const int DetID) {
// ********************************************
  std::vector<XsecNorms4> returnVec;
  int norm_counter = 0;
  IterateOverParams(DetID,
    [&](int i) { return GetParamType(i) == kNorm; }, // Filter condition
    [&](auto) {
      XsecNorms4 Temp = NormParams[norm_counter];
      returnVec.push_back(Temp);
      norm_counter++;
    }
  );
  return returnVec;
}

// ********************************************
// DB Grab the number of parameters for the relevant DetID
int covarianceXsec::GetNumParamsFromDetID(const int DetID, const SystType Type) {
// ********************************************
  int returnVal = 0;
  IterateOverParams(DetID,
    [&](int i) { return GetParamType(i) == Type; }, // Filter condition
    [&](int) { returnVal += 1; } // Action to perform if filter passes
  );
  return returnVal;
}

// ********************************************
// DB Grab the parameter names for the relevant DetID
const std::vector<std::string> covarianceXsec::GetParsNamesFromDetID(const int DetID, const SystType Type) {
// ********************************************
  std::vector<std::string> returnVec;
  IterateOverParams(DetID,
    [&](int i) { return GetParamType(i) == Type; }, // Filter condition
    [&](int i) { returnVec.push_back(GetParFancyName(i)); } // Action to perform if filter passes
  );
  return returnVec;
}

// ********************************************
// DB DB Grab the parameter indices for the relevant DetID
const std::vector<int> covarianceXsec::GetParsIndexFromDetID(const int DetID, const SystType Type) {
// ********************************************
  std::vector<int> returnVec;
  IterateOverParams(DetID,
    [&](int i) { return GetParamType(i) == Type; }, // Filter condition
    [&](int i) { returnVec.push_back(i); } // Action to perform if filter passes
  );
  return returnVec;
}

// ********************************************
template <typename FilterFunc, typename ActionFunc>
void covarianceXsec::IterateOverParams(const int DetID, FilterFunc filter, ActionFunc action) {
// ********************************************
  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetParDetID(i) & DetID) && filter(i)) { // Common filter logic
      action(i); // Specific action for each function
    }
  }
}

// ********************************************
void covarianceXsec::initParams() {
// ********************************************
  for (int i = 0; i < _fNumPar; ++i) {
    //ETA - set the name to be xsec_% as this is what ProcessorMCMC expects
    _fNames[i] = "xsec_"+std::to_string(i);

    // Set covarianceBase parameters (Curr = current, Prop = proposed, Sigma = step)
    _fCurrVal[i] = _fPreFitValue[i];
    _fPropVal[i] = _fCurrVal[i];
  }
  //DB Set Individual Step scale for PCA parameters to the LastPCAdpar fIndivStepScale because the step scale for those parameters is set by 'eigen_values[i]' but needs an overall step scale
  //   However, individual step scale for non-PCA parameters needs to be set correctly
  if (pca) {
    for (int i = FirstPCAdpar; i <= LastPCAdpar; i++) {
      _fIndivStepScale[i] = _fIndivStepScale[LastPCAdpar-1];
    }
  }
  randomize();
  //KS: Transfer the starting parameters to the PCA basis, you don't want to start with zero..
  if (pca)
  {
    TransferToPCA();
    for (int i = 0; i < _fNumParPCA; ++i) {
      _fPreFitValue_PCA[i] = fParCurr_PCA(i);
    }
  }
}

// ********************************************
// Print everything we know about the inputs we're Getting
void covarianceXsec::Print() {
// ********************************************
  MACH3LOG_INFO("#################################################");
  MACH3LOG_INFO("Printing covarianceXsec:");

  PrintGlobablInfo();

  PrintNormParams();

  PrintSplineParams();

  PrintFunctionalParams();

  PrintParameterGroups();

  MACH3LOG_INFO("Finished");
  MACH3LOG_INFO("#################################################");

  CheckCorrectInitialisation();
} // End


// ********************************************
void covarianceXsec::PrintGlobablInfo() {
// ********************************************
  MACH3LOG_INFO("============================================================================================================================================================");
  MACH3LOG_INFO("{:<5} {:2} {:<40} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<5} {:2} {:<10}", "#", "|", "Name", "|", "Nom.", "|", "Prior", "|", "Error", "|", "Lower", "|", "Upper", "|", "StepScale", "|", "DetID", "|", "Type");
  MACH3LOG_INFO("------------------------------------------------------------------------------------------------------------------------------------------------------------");
  for (int i = 0; i < GetNumParams(); i++) {
    #ifndef USE_FPGA
      std::string ErrString = fmt::format("{:.2f}", _fError[i]);
      MACH3LOG_INFO("{:<5} {:2} {:<40} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<5} {:2} {:<10}", i, "|", GetParFancyName(i), "|", _fGenerated[i], "|", _fPreFitValue[i], "|", "+/- " + ErrString, "|", _fLowBound[i], "|", _fUpBound[i], "|", _fIndivStepScale[i], "|", _fDetID[i], "|", SystType_ToString(_fParamType[i]));
    #endif
  }
  MACH3LOG_INFO("============================================================================================================================================================");
}

// ********************************************
void covarianceXsec::PrintNormParams() {
// ********************************************
  // Output the normalisation parameters as a sanity check!
  MACH3LOG_INFO("Normalisation parameters:  {}", NormParams.size());
  if(_fSystToGlobalSystIndexMap[kNorm].size() == 0) return;

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
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┴────────────────────┴────────────────────┴────────────────────┘");
}

// ********************************************
void covarianceXsec::PrintSplineParams() {
// ********************************************
  MACH3LOG_INFO("Spline parameters: {}", _fSystToGlobalSystIndexMap[kSpline].size());
  if(_fSystToGlobalSystIndexMap[kSpline].size() == 0) return;
  MACH3LOG_INFO("=====================================================================================================================================================================");
  MACH3LOG_INFO("{:<4} {:<2} {:<40} {:<2} {:<40} {:<2} {:<20} {:<2} {:<20} {:<2} {:<20} {:<2}", "#", "|", "Name", "|", "Spline Name", "|", "Spline Interpolation", "|", "Low Knot Bound", "|", "Up Knot Bound", "|");
  MACH3LOG_INFO("---------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  for (auto &pair : _fSystToGlobalSystIndexMap[kSpline]) {
    auto &SplineIndex = pair.first;
    auto &GlobalIndex = pair.second;

    MACH3LOG_INFO("{:<4} {:<2} {:<40} {:<2} {:<40} {:<2} {:<20} {:<2} {:<20} {:<2} {:<20} {:<2}",
                  SplineIndex, "|", GetParFancyName(GlobalIndex), "|",
                  _fSplineNames[SplineIndex], "|",
                  SplineInterpolation_ToString(GetParSplineInterpolation(SplineIndex)), "|",
                  GetParSplineKnotLowerBound(SplineIndex), "|",
                  GetParSplineKnotUpperBound(SplineIndex), "|");
  }
  MACH3LOG_INFO("=====================================================================================================================================================================");
}

// ********************************************
void covarianceXsec::PrintFunctionalParams() {
// ********************************************
  MACH3LOG_INFO("Functional parameters: {}", _fSystToGlobalSystIndexMap[kFunc].size());
  if(_fSystToGlobalSystIndexMap[kFunc].size() == 0) return;
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│", "#", "Global #", "Name");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┤");
  for (auto &pair : _fSystToGlobalSystIndexMap[kFunc]) {
    auto &FuncIndex = pair.first;
    auto &GlobalIndex = pair.second;
    MACH3LOG_INFO("│{0:4}│{1:<10}│{2:40}│", std::to_string(FuncIndex), GlobalIndex, GetParFancyName(GlobalIndex));
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┘");
}

// ********************************************
void covarianceXsec::PrintParameterGroups() {
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
// KS: Check if matrix is correctly initialised
void covarianceXsec::CheckCorrectInitialisation() {
// ********************************************
  // KS: Lambda Function which simply checks if there are no duplicates in std::vector
  auto CheckForDuplicates = [](const std::vector<std::string>& names, const std::string& nameType) {
    std::unordered_map<std::string, int> seenStrings;
    for (size_t i = 0; i < names.size(); ++i) {
      const auto& name = names[i];
      if (seenStrings.find(name) != seenStrings.end()) {
        int firstIndex = seenStrings[name];
        MACH3LOG_CRITICAL("There are two systematics with the same {} '{}', first at index {}, and again at index {}", nameType, name, firstIndex, i);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      seenStrings[name] = i;
    }
  };

  // KS: Checks if there are no duplicates in fancy names etc, this can happen if we merge configs etc
  CheckForDuplicates(_fFancyNames, "_fFancyNames");
  CheckForDuplicates(_fSplineNames, "_fSplineNames");
}

// ********************************************
// Function to set to prior parameters of a given group
void covarianceXsec::SetGroupOnlyParameters(const std::string& Group) {
// ********************************************
  if(!pca) {
    for (int i = 0; i < _fNumPar; i++) {
      if(IsParFromGroup(i, Group)) _fPropVal[i] = _fPreFitValue[i];
    }
  } else {
    MACH3LOG_ERROR("SetGroupOnlyParameters not implemented for PCA");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ********************************************
// Checks if parameter belongs to a given group
bool covarianceXsec::IsParFromGroup(const int i, const std::string& Group) {
// ********************************************

  if(Group == _ParameterGroup[i]) return true;
  else return false;
}

// ********************************************
// Dump Matrix to ROOT file, useful when we need to pass matrix info to another fitting group
void covarianceXsec::DumpMatrixToFile(const std::string& Name) {
// ********************************************
  TFile* outputFile = new TFile(Name.c_str(), "RECREATE");

  TObjArray* xsec_param_names = new TObjArray();
  TObjArray* xsec_spline_interpolation = new TObjArray();
  TObjArray* xsec_spline_names = new TObjArray();

  TVectorD* xsec_param_prior = new TVectorD(_fNumPar);
  TVectorD* xsec_flat_prior = new TVectorD(_fNumPar);
  TVectorD* xsec_stepscale = new TVectorD(_fNumPar);
  TVectorD* xsec_param_nom = new TVectorD(_fNumPar);
  TVectorD* xsec_param_lb = new TVectorD(_fNumPar);
  TVectorD* xsec_param_ub = new TVectorD(_fNumPar);

  TVectorD* xsec_param_knot_weight_lb = new TVectorD(_fNumPar);
  TVectorD* xsec_param_knot_weight_ub = new TVectorD(_fNumPar);

  TVectorD* xsec_error = new TVectorD(_fNumPar);
  TVectorD* xsec_param_id = new TVectorD(_fNumPar);

  for(int i = 0; i < _fNumPar; ++i)
  {
    TObjString* nameObj = new TObjString(_fFancyNames[i].c_str());
    xsec_param_names->AddLast(nameObj);

    TObjString* splineType = new TObjString("TSpline3");
    xsec_spline_interpolation->AddLast(splineType);

    TObjString* splineName = new TObjString("");
    xsec_spline_names->AddLast(splineName);

    (*xsec_param_prior)[i] = _fPreFitValue[i];
    (*xsec_param_nom)[i] = _fGenerated[i];
    (*xsec_flat_prior)[i] = _fFlatPrior[i];
    (*xsec_stepscale)[i] = _fIndivStepScale[i];
    (*xsec_error)[i] = _fError[i];
    (*xsec_param_id)[i] = _fDetID[i];

    (*xsec_param_lb)[i] = _fLowBound[i];
    (*xsec_param_ub)[i] = _fUpBound[i];

    //Default values
    (*xsec_param_knot_weight_lb)[i] = -9999;
    (*xsec_param_knot_weight_ub)[i] = +9999;
  }

  for (auto &pair : _fSystToGlobalSystIndexMap[kSpline]) {
    auto &SplineIndex = pair.first;
    auto &SystIndex = pair.second;

    (*xsec_param_knot_weight_lb)[SystIndex] = SplineParams.at(SplineIndex)._SplineKnotLowBound;
    (*xsec_param_knot_weight_ub)[SystIndex] = SplineParams.at(SplineIndex)._SplineKnotUpBound;

    TObjString* splineType = new TObjString(SplineInterpolation_ToString(SplineParams.at(SplineIndex)._SplineInterpolationType).c_str());
    xsec_spline_interpolation->AddAt(splineType, SystIndex);

    TObjString* splineName = new TObjString(_fSplineNames[SplineIndex].c_str());
    xsec_spline_names->AddAt(splineName, SystIndex);
  }
  xsec_param_names->Write("xsec_param_names", TObject::kSingleKey);
  delete xsec_param_names;
  xsec_spline_interpolation->Write("xsec_spline_interpolation", TObject::kSingleKey);
  delete xsec_spline_interpolation;
  xsec_spline_names->Write("xsec_spline_names", TObject::kSingleKey);
  delete xsec_spline_names;

  xsec_param_prior->Write("xsec_param_prior");
  delete xsec_param_prior;
  xsec_flat_prior->Write("xsec_flat_prior");
  delete xsec_flat_prior;
  xsec_stepscale->Write("xsec_stepscale");
  delete xsec_stepscale;
  xsec_param_nom->Write("xsec_param_nom");
  delete xsec_param_nom;
  xsec_param_lb->Write("xsec_param_lb");
  delete xsec_param_lb;
  xsec_param_ub->Write("xsec_param_ub");
  delete xsec_param_ub;

  xsec_param_knot_weight_lb->Write("xsec_param_knot_weight_lb");
  delete xsec_param_knot_weight_lb;
  xsec_param_knot_weight_ub->Write("xsec_param_knot_weight_ub");
  delete xsec_param_knot_weight_ub;

  xsec_param_id->Write("xsec_param_id");
  delete xsec_param_id;
  xsec_error->Write("xsec_error");
  delete xsec_error;

  covMatrix->Write("xsec_cov");
  TH2D* CorrMatrix = GetCorrelationMatrix();
  CorrMatrix->Write("hcov");
  delete CorrMatrix;

  outputFile->Close();
  delete outputFile;

  MACH3LOG_INFO("Finished dumping covariance object");
}
