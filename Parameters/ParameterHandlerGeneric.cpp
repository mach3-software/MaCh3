#include "Parameters/ParameterHandlerGeneric.h"

// ********************************************
// ETA - YAML constructor
// this will replace the root file constructor but let's keep it in
// to do some validations
ParameterHandlerGeneric::ParameterHandlerGeneric(
    const std::vector<std::string> &YAMLFile, std::string name,
    double threshold, int FirstPCA, int LastPCA)
    : ParameterHandlerBase(YAMLFile, name, threshold, FirstPCA, LastPCA) {
  // ********************************************
  InitParametersTypeFromConfig();

  // ETA - again this really doesn't need to be here...
  for (int i = 0; i < parlist.NumParameters(); i++) {
    // Sort out the print length
    if (int(parlist.params.name[i].length()) > PrintLength)
      PrintLength = int(parlist.params.name[i].length());
  } // end the for loop

  MACH3LOG_DEBUG("Constructing instance of ParameterHandler");
  InitParams();
  // Print
  Print();
}

// ********************************************
void ParameterHandlerGeneric::InitParametersTypeFromConfig() {
  // ********************************************

  // ETA - read in the systematics. Would be good to add in some checks to make
  // sure that there are the correct number of entries i.e. are the _fNumPars
  // for Names, PreFitValues etc etc.
  for (auto const &node : _fYAMLDoc["Systematics"]) {
    auto const &pardef = node["Systematic"];

    auto group = Get<std::string>(pardef["ParameterGroup"], __FILE__, __LINE__);
    auto fancy_name =
        Get<std::string>(pardef["Names"]["FancyName"], __FILE__, __LINE__);

    // Fill the map to get the correlations later as well
    auto ParamType =
        String_ToSystType(Get<std::string>(pardef["Type"], __FILE__, __LINE__));

    // Now load in variables for spline systematics only
    switch (ParamType) {
    case kSpline: {
      SplineParams.push_back(GetSplineParameter(pardef["SplineInformation"]));
      SplineParams.back().group = group;
      SplineParams.back().fancy_name = fancy_name;
      break;
    }
    case kNorm: {
      NormParams.push_back(GetNormParameter(pardef));
      NormParams.back().group = group;
      NormParams.back().fancy_name = fancy_name;
      break;
    }
    case kFunc: {
      FuncParams.push_back(GetFunctionalParameters(pardef));
      FuncParams.back().group = group;
      FuncParams.back().fancy_name = fancy_name;
      break;
    }
    case kOsc: {
      OscParams.push_back(GetOscillationParameters(pardef));
      OscParams.back().group = group;
      OscParams.back().fancy_name = fancy_name;
      break;
    }
    case kSystTypes:
    default: {
      break;
    }
    }

  } // end loop over params

  DetermineGlobalParameterIndices();
}

void ParameterHandlerGeneric::DetermineGlobalParameterIndices() {

  GlobalParams.clear();
  for (auto &par : SplineParams) {
    par.index = parlist.FindParameterByFancyName(par.fancy_name);
    GlobalParams.push_back(&par);
  }
  for (auto &par : NormParams) {
    par.index = parlist.FindParameterByFancyName(par.fancy_name);
    GlobalParams.push_back(&par);

    if (GetLowerBound(par.index) < 0.) {
      MACH3LOG_ERROR("Normalisation Parameter {} ({}), has lower parameters "
                     "bound which can go below 0 and is equal {}",
                     GetParFancyName(par.index), par.index,
                     GetLowerBound(par.index));
      MACH3LOG_ERROR(
          "Normalisation parameters can't go bellow 0 as this is unphysical");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  for (auto &par : FuncParams) {
    par.index = parlist.FindParameterByFancyName(par.fancy_name);
    par.valuePtr = RetPointer(par.index);
    GlobalParams.push_back(&par);
  }
  for (auto &par : OscParams) {
    par.index = parlist.FindParameterByFancyName(par.fancy_name);
    par.valuePtr = RetPointer(par.index);
    GlobalParams.push_back(&par);
  }

  // LP order by index... probably doesn't help but probably doesn't hurt
  std::sort(GlobalParams.begin(), GlobalParams.end(),
            [](TypeParameterBase const *l, TypeParameterBase const *r) {
              return l->index < r->index;
            });
}

TypeParameterBase const *ParameterHandlerGeneric::GetParam(const int i) const {
  auto it = std::find_if(
      GlobalParams.begin(), GlobalParams.end(),
      [=](TypeParameterBase const *par) { return par->index == i; });
  if (it == GlobalParams.end()) {
    return nullptr;
  }
  return *it;
}

SystType ParameterHandlerGeneric::GetParamType(const int i) const {
  auto par = GetParam(i);
  return par ? par->syst_type : kSystTypes;
}

// ********************************************
ParameterHandlerGeneric::~ParameterHandlerGeneric() {
  // ********************************************
  MACH3LOG_DEBUG("Destructing ParameterHandler");
}

// ********************************************
// DB Grab the Spline Names for the relevant SampleName
const std::vector<std::string>
ParameterHandlerGeneric::GetSplineParsNamesFromSampleName(
    const std::string &SampleName) {
  // ********************************************
  std::vector<std::string> returnVec;
  for (auto &sp : SplineParams) {
    if (AppliesToSample(sp.index, SampleName)) {
      returnVec.push_back(sp.spline_name);
    }
  }
  return returnVec;
}

// ********************************************
const std::vector<SplineInterpolation>
ParameterHandlerGeneric::GetSplineInterpolationFromSampleName(
    const std::string &SampleName) {
  // ********************************************
  std::vector<SplineInterpolation> returnVec;
  for (auto &sp : SplineParams) {
    if (AppliesToSample(sp.index, SampleName)) {
      returnVec.push_back(sp.SplineInterpolationType);
    }
  }
  return returnVec;
}

// ********************************************
// DB Grab the Spline Modes for the relevant SampleName
const std::vector<std::vector<int>>
ParameterHandlerGeneric::GetSplineModeVecFromSampleName(
    const std::string &SampleName) {
  // ********************************************
  std::vector<std::vector<int>> returnVec;
  for (auto &sp : SplineParams) {
    if (AppliesToSample(sp.index, SampleName)) {
      returnVec.push_back(sp.SplineModes);
    }
  }
  return returnVec;
}

// ********************************************
// Get Norm params
NormParameter
ParameterHandlerGeneric::GetNormParameter(const YAML::Node &param) {
  // ********************************************
  NormParameter norm;

  /// ETA size 0 to mean apply to all
  /// Ultimately all this information ends up in the @NormParams vector
  norm.modes =
      GetFromManager<std::vector<int>>(param["Mode"], {}, __FILE__, __LINE__);
  norm.pdgs = GetFromManager<std::vector<int>>(param["NeutrinoFlavour"], {},
                                               __FILE__, __LINE__);
  norm.preoscpdgs = GetFromManager<std::vector<int>>(
      param["NeutrinoFlavourUnosc"], {}, __FILE__, __LINE__);
  norm.targets = GetFromManager<std::vector<int>>(param["TargetNuclei"], {},
                                                  __FILE__, __LINE__);
  int NumKinematicCuts = 0;
  if (param["KinematicCuts"]) {
    NumKinematicCuts = int(param["KinematicCuts"].size());

    std::vector<std::string> TempKinematicStrings;
    std::vector<std::vector<std::vector<double>>> TempKinematicBounds;
    // First element of TempKinematicBounds is always -999, and size is then 3
    for (int KinVar_i = 0; KinVar_i < NumKinematicCuts; ++KinVar_i) {
      // ETA: This is a bit messy, Kinematic cuts is a list of maps
      for (YAML::const_iterator it = param["KinematicCuts"][KinVar_i].begin();
           it != param["KinematicCuts"][KinVar_i].end(); ++it) {
        TempKinematicStrings.push_back(it->first.as<std::string>());
        TempKinematicBounds.push_back(Get2DBounds(it->second));
      }
      if (TempKinematicStrings.size() == 0) {
        MACH3LOG_ERROR("Received a KinematicCuts node but couldn't read the "
                       "contents (it's a list of single-element dictionaries "
                       "(python) = map of pairs (C++))");
        MACH3LOG_ERROR("For Param {}", norm.fancy_name);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    } // KinVar_i
    norm.KinematicVarStr = TempKinematicStrings;
    norm.Selection = TempKinematicBounds;
  }

  // Next ones are kinematic bounds on where normalisation parameter should
  // apply We set a bool to see if any bounds exist so we can short-circuit
  // checking all of them every step
  bool HasKinBounds = false;

  if (norm.KinematicVarStr.size() > 0)
    HasKinBounds = true;

  norm.hasKinBounds = HasKinBounds;
  // End of kinematic bound checking

  return norm;
}

// ********************************************
// Grab the global syst index for the relevant SampleName
// i.e. get a vector of size nSplines where each entry is filled with the global
// syst number
const std::vector<int>
ParameterHandlerGeneric::GetGlobalSystIndexFromSampleName(
    const std::string &SampleName, const SystType Type) {
  // ********************************************
  std::vector<int> returnVec;
  for (auto par : GlobalParams) {
    if (par->syst_type != Type) {
      continue;
    }
    if (AppliesToSample(par->index, SampleName)) {
      returnVec.push_back(par->index);
    }
  }
  return returnVec;
}

// ********************************************
// Get Norm params
SplineParameter
ParameterHandlerGeneric::GetSplineParameter(const YAML::Node &param) {
  // ********************************************

  SplineParameter Spline;

  // Now get the Spline interpolation type
  if (param["InterpolationType"]) {
    for (int InterpType = 0; InterpType < kSplineInterpolations; ++InterpType) {
      if (param["InterpolationType"].as<std::string>() ==
          SplineInterpolation_ToString(SplineInterpolation(InterpType)))
        Spline.SplineInterpolationType = SplineInterpolation(InterpType);
    }
  } else { // KS: By default use TSpline3
    Spline.SplineInterpolationType = kTSpline3;
  }

  Spline.SplineKnotUpBound = GetFromManager<double>(
      param["SplineKnotUpBound"], M3::DefSplineKnotUpBound, __FILE__, __LINE__);
  Spline.SplineKnotLowBound =
      GetFromManager<double>(param["SplineKnotLowBound"],
                             M3::DefSplineKnotLowBound, __FILE__, __LINE__);

  if (Spline.SplineKnotUpBound != M3::DefSplineKnotUpBound ||
      Spline.SplineKnotLowBound != M3::DefSplineKnotLowBound) {
    MACH3LOG_WARN(
        "Spline knot capping enabled with bounds [{}, {}]. For reliable fits, "
        "consider modifying the input generation instead.",
        Spline.SplineKnotLowBound, Spline.SplineKnotUpBound);
  }
  // If there is no mode information given then this will be an empty vector
  Spline.SplineModes =
      GetFromManager(param["Mode"], std::vector<int>(), __FILE__, __LINE__);

  if (param["SplineName"]) {
    Spline.spline_name = param["SplineName"].as<std::string>();
  }

  return Spline;
}

// ********************************************
// Get Func params
FunctionalParameter
ParameterHandlerGeneric::GetFunctionalParameters(const YAML::Node &param) {
  // ********************************************
  FunctionalParameter func;

  func.pdgs = GetFromManager<std::vector<int>>(
      param["NeutrinoFlavour"], std::vector<int>(), __FILE__, __LINE__);
  func.targets = GetFromManager<std::vector<int>>(
      param["TargetNuclei"], std::vector<int>(), __FILE__, __LINE__);
  func.modes = GetFromManager<std::vector<int>>(
      param["Mode"], std::vector<int>(), __FILE__, __LINE__);
  func.preoscpdgs = GetFromManager<std::vector<int>>(
      param["NeutrinoFlavourUnosc"], std::vector<int>(), __FILE__, __LINE__);

  // HH - Copied from GetXsecNorm
  int NumKinematicCuts = 0;
  if (param["KinematicCuts"]) {

    NumKinematicCuts = int(param["KinematicCuts"].size());

    std::vector<std::string> TempKinematicStrings;
    std::vector<std::vector<std::vector<double>>> TempKinematicBounds;
    // First element of TempKinematicBounds is always -999, and size is then 3
    for (int KinVar_i = 0; KinVar_i < NumKinematicCuts; ++KinVar_i) {
      // ETA: This is a bit messy, Kinematic cuts is a list of maps
      for (YAML::const_iterator it = param["KinematicCuts"][KinVar_i].begin();
           it != param["KinematicCuts"][KinVar_i].end(); ++it) {
        TempKinematicStrings.push_back(it->first.as<std::string>());
        TempKinematicBounds.push_back(Get2DBounds(it->second));
      }
      if (TempKinematicStrings.size() == 0) {
        MACH3LOG_ERROR("Received a KinematicCuts node but couldn't read the "
                       "contents (it's a list of single-element dictionaries "
                       "(python) = map of pairs (C++))");
        MACH3LOG_ERROR("For Param {}", func.fancy_name);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    } // KinVar_i
    func.KinematicVarStr = TempKinematicStrings;
    func.Selection = TempKinematicBounds;
  }

  return func;
}

// ********************************************
// Get Osc params
OscillationParameter
ParameterHandlerGeneric::GetOscillationParameters(const YAML::Node &) {
  // ********************************************
  return OscillationParameter();
}

// ********************************************
// HH: Grab the Functional parameters for the relevant SampleName
const std::vector<FunctionalParameter>
ParameterHandlerGeneric::GetFunctionalParametersFromSampleName(
    const std::string &SampleName) const {
  // ********************************************
  std::vector<FunctionalParameter> returnVec;
  for (auto &par : FuncParams) {
    if (AppliesToSample(par.index, SampleName)) {
      returnVec.push_back(par);
    }
  }
  return returnVec;
}

// ********************************************
// DB Grab the Normalisation parameters for the relevant SampleName
const std::vector<NormParameter>
ParameterHandlerGeneric::GetNormParsFromSampleName(
    const std::string &SampleName) const {
  // ********************************************
  std::vector<NormParameter> returnVec;
  for (auto &par : NormParams) {
    if (AppliesToSample(par.index, SampleName)) {
      returnVec.push_back(par);
    }
  }
  return returnVec;
}

// ********************************************
// KS Grab the Spline parameters for the relevant SampleName
const std::vector<SplineParameter>
ParameterHandlerGeneric::GetSplineParsFromSampleName(
    const std::string &SampleName) const {
  // ********************************************
  std::vector<SplineParameter> returnVec;
  for (auto &par : SplineParams) {
    if (AppliesToSample(par.index, SampleName)) {
      returnVec.push_back(par);
    }
  }
  return returnVec;
}

// ********************************************
// DB Grab the number of parameters for the relevant SampleName
int ParameterHandlerGeneric::GetNumParamsFromSampleName(
    const std::string &SampleName, const SystType Type) {
  // ********************************************
  int returnVal = 0;
  IterateOverParams(
      SampleName,
      [&](int i) { return GetParamType(i) == Type; }, // Filter condition
      [&](int) { returnVal += 1; } // Action to perform if filter passes
  );
  return returnVal;
}

// ********************************************
// DB Grab the parameter names for the relevant SampleName
const std::vector<std::string>
ParameterHandlerGeneric::GetParsNamesFromSampleName(
    const std::string &SampleName, const SystType Type) {
  // ********************************************
  std::vector<std::string> returnVec;
  IterateOverParams(
      SampleName,
      [&](int i) { return GetParamType(i) == Type; }, // Filter condition
      [&](int i) {
        returnVec.push_back(GetParFancyName(i));
      } // Action to perform if filter passes
  );
  return returnVec;
}

// ********************************************
// DB DB Grab the parameter indices for the relevant SampleName
const std::vector<int> ParameterHandlerGeneric::GetParsIndexFromSampleName(
    const std::string &SampleName, const SystType Type) {
  // ********************************************
  std::vector<int> returnVec;
  IterateOverParams(
      SampleName,
      [&](int i) { return GetParamType(i) == Type; }, // Filter condition
      [&](int i) {
        returnVec.push_back(i);
      } // Action to perform if filter passes
  );
  return returnVec;
}

// ********************************************
template <typename FilterFunc, typename ActionFunc>
void ParameterHandlerGeneric::IterateOverParams(const std::string &SampleName,
                                                FilterFunc const &filter,
                                                ActionFunc const &action) {
  // ********************************************
  for (auto &par : GlobalParams) {
    if ((AppliesToSample(par->index, SampleName)) &&
        filter(par->index)) { // Common filter logic
      action(par->index);     // Specific action for each function
    }
  }
}

// ********************************************
void ParameterHandlerGeneric::InitParams() {
  // ********************************************
  for (int i = 0; i < parlist.NumParameters(); ++i) {
    // ETA - set the name to be xsec_% as this is what ProcessorMCMC expects
    parlist.params.name[i] = "xsec_" + std::to_string(i);

    // KS: Plenty
    if (GetParamType(i) == kOsc) {
      parlist.params.name[i] = parlist.params.fancy_name[i];

      if (GetParam(i)->group != "Osc") {
        MACH3LOG_ERROR("Parameter {}, is of type Oscillation but doesn't "
                       "belong to Osc group",
                       parlist.params.fancy_name[i]);

        MACH3LOG_ERROR("It belongs to {} group", GetParam(i)->group);

        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    // Set ParameterHandler parameters (Curr = current, Prop = proposed, Sigma =
    // step)
    SetParCurrProp(i, GetParInit(i));
  }
}

// ********************************************
// Print everything we know about the inputs we're Getting
void ParameterHandlerGeneric::Print() {
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
void ParameterHandlerGeneric::PrintGlobablInfo() {
  // ********************************************
  MACH3LOG_INFO("=============================================================="
                "=============================================================="
                "================================");
  MACH3LOG_INFO("{:<5} {:2} {:<40} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} "
                "{:<10} {:2} {:<10} {:2} {:<20} {:2} {:<10}",
                "#", "|", "Name", "|", "Prior", "|", "Error", "|", "Lower", "|",
                "Upper", "|", "StepScale", "|", "SampleNames", "|", "Type");
  MACH3LOG_INFO("--------------------------------------------------------------"
                "--------------------------------------------------------------"
                "--------------------------------");
  for (int i = 0; i < GetNumSystematicParams(); i++) {
    std::string ErrString = fmt::format("{:.2f}", GetError(i));
    std::string SampleNameString = "";
    for (const auto &SampleName : parlist.params.samples[i]) {
      if (!SampleNameString.empty()) {
        SampleNameString += ", ";
      }
      SampleNameString += SampleName;
    }
    MACH3LOG_INFO("{:<5} {:2} {:<40} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} "
                  "{:<10} {:2} {:<10} {:2} {:<20} {:2} {:<10}",
                  i, "|", GetParFancyName(i), "|", GetParInit(i), "|",
                  "+/- " + ErrString, "|", GetLowerBound(i), "|",
                  GetUpperBound(i), "|", GetIndivStepScale(i), "|",
                  SampleNameString, "|", GetParamTypeString(i));
  }
  MACH3LOG_INFO("=============================================================="
                "=============================================================="
                "================================");
}

// ********************************************
void ParameterHandlerGeneric::PrintNormParams() {
  // ********************************************
  // Output the normalisation parameters as a sanity check!
  MACH3LOG_INFO("Normalisation parameters:  {}", NormParams.size());
  if (!NormParams.size()) {
    return;
  }

  bool have_parameter_with_kin_bounds = false;

  // KS: Consider making some class producing table..
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┬────"
                "────────────────┬────────────────────┬────────────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│{3:20}│{4:20}│{5:20}│", "#", "Global #",
                "Name", "Int. mode", "Target", "pdg");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┼────"
                "────────────────┼────────────────────┼────────────────────┤");

  for (unsigned int i = 0; i < NormParams.size(); ++i) {
    std::string intModeString;
    for (unsigned int j = 0; j < NormParams[i].modes.size(); j++) {
      intModeString += std::to_string(NormParams[i].modes[j]);
      intModeString += " ";
    }
    if (NormParams[i].modes.empty())
      intModeString += "all";

    std::string targetString;
    for (unsigned int j = 0; j < NormParams[i].targets.size(); j++) {
      targetString += std::to_string(NormParams[i].targets[j]);
      targetString += " ";
    }
    if (NormParams[i].targets.empty())
      targetString += "all";

    std::string pdgString;
    for (unsigned int j = 0; j < NormParams[i].pdgs.size(); j++) {
      pdgString += std::to_string(NormParams[i].pdgs[j]);
      pdgString += " ";
    }
    if (NormParams[i].pdgs.empty())
      pdgString += "all";

    MACH3LOG_INFO("│{: <4}│{: <10}│{: <40}│{: <20}│{: <20}│{: <20}│", i,
                  NormParams[i].index, NormParams[i].fancy_name, intModeString,
                  targetString, pdgString);

    if (NormParams[i].hasKinBounds)
      have_parameter_with_kin_bounds = true;
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┴────"
                "────────────────┴────────────────────┴────────────────────┘");

  if (have_parameter_with_kin_bounds) {
    MACH3LOG_INFO("Normalisation parameters KinematicCuts information");
    MACH3LOG_INFO(
        "┌────┬──────────┬────────────────────────────────────────┬────────────"
        "────────┬────────────────────────────────────────┐");
    MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│{3:20}│{4:40}│", "#", "Global #",
                  "Name", "KinematicCut", "Value");
    MACH3LOG_INFO(
        "├────┼──────────┼────────────────────────────────────────┼────────────"
        "────────┼────────────────────────────────────────┤");
    for (unsigned int i = 0; i < NormParams.size(); ++i) {
      // skip parameters with no KinematicCuts
      if (!NormParams[i].hasKinBounds)
        continue;

      const long unsigned int ncuts = NormParams[i].KinematicVarStr.size();
      for (long unsigned int icut = 0; icut < ncuts; icut++) {
        std::string kinematicCutValueString;
        for (const auto &value : NormParams[i].Selection[icut]) {
          for (const auto &v : value) {
            kinematicCutValueString += fmt::format("{:.2f} ", v);
          }
        }
        if (icut == 0)
          MACH3LOG_INFO("│{: <4}│{: <10}│{: <40}│{: <20}│{: <40}│", i,
                        NormParams[i].index, NormParams[i].fancy_name,
                        NormParams[i].KinematicVarStr[icut],
                        kinematicCutValueString);
        else
          MACH3LOG_INFO("│{: <4}│{: <10}│{: <40}│{: <20}│{: <40}│", "", "", "",
                        NormParams[i].KinematicVarStr[icut],
                        kinematicCutValueString);
      } // icut
    } // i
    MACH3LOG_INFO(
        "└────┴──────────┴────────────────────────────────────────┴────────────"
        "────────┴────────────────────────────────────────┘");
  } else
    MACH3LOG_INFO("No normalisation parameters have KinematicCuts defined");
}

// ********************************************
void ParameterHandlerGeneric::PrintSplineParams() {
  // ********************************************
  MACH3LOG_INFO("Spline parameters: {}", SplineParams.size());
  if (SplineParams.size() == 0) {
    return;
  }
  MACH3LOG_INFO("=============================================================="
                "=============================================================="
                "=========================================");
  MACH3LOG_INFO("{:<4} {:<2} {:<40} {:<2} {:<40} {:<2} {:<20} {:<2} {:<20} "
                "{:<2} {:<20} {:<2}",
                "#", "|", "Name", "|", "Spline Name", "|",
                "Spline Interpolation", "|", "Low Knot Bound", "|",
                "Up Knot Bound", "|");
  MACH3LOG_INFO("--------------------------------------------------------------"
                "--------------------------------------------------------------"
                "-----------------------------------------");
  for (int i = 0; i < int(SplineParams.size()); ++i) {
    MACH3LOG_INFO(
        "{:<4} {:<2} {:<40} {:<2} {:<40} {:<2} {:<20} {:<2} {:<20} {:<2} "
        "{:<20} {:<2}",
        i, "|", GetParFancyName(SplineParams[i].index), "|",
        SplineParams[i].spline_name, "|",
        SplineInterpolation_ToString(GetParSplineInterpolation(i)), "|",
        GetParSplineKnotLowerBound(i), "|", GetParSplineKnotUpperBound(i), "|");
  }
  MACH3LOG_INFO("=============================================================="
                "=============================================================="
                "=========================================");
}

// ********************************************
void ParameterHandlerGeneric::PrintFunctionalParams() {
  // ********************************************
  MACH3LOG_INFO("Functional parameters: {}", FuncParams.size());
  if (FuncParams.size() == 0) {
    return;
  }
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│", "#", "Global #", "Name");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┤");
  for (int i = 0; i < int(FuncParams.size()); ++i) {
    MACH3LOG_INFO("│{0:4}│{1:<10}│{2:40}│", std::to_string(i),
                  FuncParams[i].index, GetParFancyName(FuncParams[i].index));
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┘");
}

// ********************************************
void ParameterHandlerGeneric::PrintOscillationParams() {
  // ********************************************
  MACH3LOG_INFO("Oscillation parameters: {}", OscParams.size());
  if (OscParams.size() == 0) {
    return;
  }
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│", "#", "Global #", "Name");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┤");
  for (int i = 0; i < int(OscParams.size()); ++i) {
    MACH3LOG_INFO("│{0:4}│{1:<10}│{2:40}│", std::to_string(i),
                  OscParams[i].index, GetParFancyName(OscParams[i].index));
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┘");
}

// ********************************************
void ParameterHandlerGeneric::PrintParameterGroups() {
  // ********************************************
  // KS: Create a map to store the counts of unique strings, in principle this
  // could be in header file
  std::unordered_map<std::string, int> paramCounts;

  for (auto &par : GlobalParams) {
    paramCounts[par->group]++;
  }

  MACH3LOG_INFO("Printing parameter groups");
  // Output the counts
  for (const auto &pair : paramCounts) {
    MACH3LOG_INFO("Found {}: {} params", pair.second, pair.first);
  }
}

// ********************************************
std::vector<std::string> ParameterHandlerGeneric::GetUniqueParameterGroups() {
  // ********************************************
  std::unordered_set<std::string> uniqueGroups;

  // Fill the set with unique values
  for (auto &par : GlobalParams) {
    uniqueGroups.insert(par->group);
  }

  // Convert to vector and return
  std::vector<std::string> result(uniqueGroups.begin(), uniqueGroups.end());
  return result;
}

// ********************************************
// KS: Check if matrix is correctly initialised
void ParameterHandlerGeneric::CheckCorrectInitialisation() {
  // ********************************************
  // KS: Lambda Function which simply checks if there are no duplicates in
  // std::vector
  auto CheckForDuplicates = [](const std::vector<std::string> &names,
                               const std::string &nameType) {
    std::unordered_map<std::string, size_t> seenStrings;
    for (size_t i = 0; i < names.size(); ++i) {
      const auto &name = names[i];
      if (seenStrings.find(name) != seenStrings.end()) {
        size_t firstIndex = seenStrings[name];
        MACH3LOG_CRITICAL("There are two systematics with the same {} '{}', "
                          "first at index {}, and again at index {}",
                          nameType, name, firstIndex, i);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      seenStrings[name] = i;
    }
  };

  // KS: Checks if there are no duplicates in fancy names etc, this can happen
  // if we merge configs etc
  CheckForDuplicates(parlist.params.fancy_name, "fancy_name");

  std::vector<std::string> spline_names;
  for (auto &sp : SplineParams) {
    spline_names.push_back(sp.spline_name);
  }
  CheckForDuplicates(spline_names, "spline_names");
}

// ********************************************
// Function to set to prior parameters of a given group
void ParameterHandlerGeneric::SetGroupOnlyParameters(
    const std::vector<std::string> &Groups) {
  // ********************************************
  for (size_t i = 0; i < Groups.size(); i++) {
    SetGroupOnlyParameters(Groups[i]);
  }
}

// ********************************************
// Function to set to prior parameters of a given group
void ParameterHandlerGeneric::SetGroupOnlyParameters(
    const std::string &Group, const std::vector<double> &Pars) {
  // ********************************************
  // If empty, set the proposed to prior
  if (Pars.empty()) {
    for (int i = 0; i < parlist.NumParameters(); i++) {
      if (IsParFromGroup(i, Group)) {
        SetParProp(i, GetParInit(i));
      }
    }
  } else {
    const size_t ExpectedSize = static_cast<size_t>(GetNumParFromGroup(Group));
    if (Pars.size() != ExpectedSize) {
      MACH3LOG_ERROR("Number of param in group {} is {}, while you passed {}",
                     Group, ExpectedSize, Pars.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    int Counter = 0;
    for (int i = 0; i < parlist.NumParameters(); i++) {
      // If belongs to group set value from parsed vector, otherwise use propose
      // value
      if (IsParFromGroup(i, Group)) {
        SetParProp(i, Pars[Counter]);
        Counter++;
      }
    }
  }
}

// ********************************************
// Toggle fix/free to parameters of a given group
void ParameterHandlerGeneric::ToggleFixGroupOnlyParameters(
    const std::string &Group) {
  // ********************************************
  for (int i = 0; i < parlist.NumParameters(); ++i)
    if (IsParFromGroup(i, Group))
      ToggleFixParameter(i);
}

// ********************************************
// Toggle fix/free to parameters of several groups
void ParameterHandlerGeneric::ToggleFixGroupOnlyParameters(
    const std::vector<std::string> &Groups) {
  // ********************************************
  for (size_t i = 0; i < Groups.size(); i++)
    ToggleFixGroupOnlyParameters(Groups[i]);
}

// ********************************************
// Set parameters to be fixed in a given group
void ParameterHandlerGeneric::SetFixGroupOnlyParameters(
    const std::string &Group) {
  // ********************************************
  for (int i = 0; i < parlist.NumParameters(); ++i)
    if (IsParFromGroup(i, Group))
      SetFixParameter(i);
}

// ********************************************
// Set parameters of several groups to be fixed
void ParameterHandlerGeneric::SetFixGroupOnlyParameters(
    const std::vector<std::string> &Groups) {
  // ********************************************
  for (size_t i = 0; i < Groups.size(); i++)
    SetFixGroupOnlyParameters(Groups[i]);
}

// ********************************************
// Set parameters to be free in a given group
void ParameterHandlerGeneric::SetFreeGroupOnlyParameters(
    const std::string &Group) {
  // ********************************************
  for (int i = 0; i < parlist.NumParameters(); ++i)
    if (IsParFromGroup(i, Group))
      SetFreeParameter(i);
}

// ********************************************
// Set parameters of several groups to be fixed
void ParameterHandlerGeneric::SetFreeGroupOnlyParameters(
    const std::vector<std::string> &Groups) {
  // ********************************************
  for (size_t i = 0; i < Groups.size(); i++)
    SetFreeGroupOnlyParameters(Groups[i]);
}

// ********************************************
// Checks if parameter belongs to a given group
bool ParameterHandlerGeneric::IsParFromGroup(const int i,
                                             const std::string &Group) const {
  // ********************************************
  std::string groupLower = Group;
  auto par = GetParam(i);
  std::string paramGroupLower = par ? par->group : std::string("unmanaged");
  // KS: Convert both strings to lowercase, this way comparison will be case
  // insensitive
  std::transform(groupLower.begin(), groupLower.end(), groupLower.begin(),
                 ::tolower);
  std::transform(paramGroupLower.begin(), paramGroupLower.end(),
                 paramGroupLower.begin(), ::tolower);

  return groupLower == paramGroupLower;
}

// ********************************************
int ParameterHandlerGeneric::GetNumParFromGroup(
    const std::string &Group) const {
  // ********************************************
  int Counter = 0;
  for (int i = 0; i < parlist.NumParameters(); i++) {
    if (IsParFromGroup(i, Group))
      Counter++;
  }
  return Counter;
}

// ********************************************
// DB Grab the Normalisation parameters for the relevant sample name
std::vector<const double *> ParameterHandlerGeneric::GetOscParsFromSampleName(
    const std::string &SampleName) {
  // ********************************************
  std::vector<const double *> returnVec;
  for (auto &par : OscParams) {
    if (AppliesToSample(par.index, SampleName)) {
      returnVec.push_back(par.valuePtr);
    }
  }
  return returnVec;
}

// ********************************************
// Dump Matrix to ROOT file, useful when we need to pass matrix info to another
// fitting group
void ParameterHandlerGeneric::DumpMatrixToFile(const std::string &Name) {
  // ********************************************
  TFile *outputFile = new TFile(Name.c_str(), "RECREATE");

  TObjArray *xsec_param_names = new TObjArray();
  TObjArray *xsec_spline_interpolation = new TObjArray();
  TObjArray *xsec_spline_names = new TObjArray();

  TVectorD *xsec_param_prior = new TVectorD(parlist.NumParameters());
  TVectorD *xsec_flat_prior = new TVectorD(parlist.NumParameters());
  TVectorD *xsec_stepscale = new TVectorD(parlist.NumParameters());
  TVectorD *xsec_param_lb = new TVectorD(parlist.NumParameters());
  TVectorD *xsec_param_ub = new TVectorD(parlist.NumParameters());

  TVectorD *xsec_param_knot_weight_lb = new TVectorD(parlist.NumParameters());
  TVectorD *xsec_param_knot_weight_ub = new TVectorD(parlist.NumParameters());
  TVectorD *xsec_error = new TVectorD(parlist.NumParameters());

  for (int i = 0; i < parlist.NumParameters(); ++i) {
    TObjString *nameObj = new TObjString(GetParFancyName(i).c_str());
    xsec_param_names->AddLast(nameObj);

    TObjString *splineType = new TObjString("TSpline3");
    xsec_spline_interpolation->AddLast(splineType);

    TObjString *splineName = new TObjString("");
    xsec_spline_names->AddLast(splineName);

    (*xsec_param_prior)[i] = GetParInit(i);
    (*xsec_flat_prior)[i] = GetFlatPrior(i);
    (*xsec_stepscale)[i] = GetIndivStepScale(i);
    (*xsec_error)[i] = GetError(i);

    (*xsec_param_lb)[i] = GetLowerBound(i);
    (*xsec_param_ub)[i] = GetUpperBound(i);

    //Default values
    (*xsec_param_knot_weight_lb)[i] = -9999;
    (*xsec_param_knot_weight_ub)[i] = +9999;
  }

  for (auto &sp : SplineParams) {

    (*xsec_param_knot_weight_lb)[sp.index] = sp.SplineKnotLowBound;
    (*xsec_param_knot_weight_ub)[sp.index] = sp.SplineKnotUpBound;

    TObjString* splineType = new TObjString(SplineInterpolation_ToString(sp.SplineInterpolationType).c_str());
    xsec_spline_interpolation->AddAt(splineType, sp.index);

    TObjString* splineName = new TObjString(sp.spline_name.c_str());
    xsec_spline_names->AddAt(splineName, sp.index);
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
  xsec_param_lb->Write("xsec_param_lb");
  delete xsec_param_lb;
  xsec_param_ub->Write("xsec_param_ub");
  delete xsec_param_ub;

  xsec_param_knot_weight_lb->Write("xsec_param_knot_weight_lb");
  delete xsec_param_knot_weight_lb;
  xsec_param_knot_weight_ub->Write("xsec_param_knot_weight_ub");
  delete xsec_param_knot_weight_ub;
  xsec_error->Write("xsec_error");
  delete xsec_error;

  covMatrix->Write("xsec_cov");
  TH2D* CorrMatrix = GetCorrelationMatrix();
  CorrMatrix->Write("hcov");
  delete CorrMatrix;

  outputFile->Close();
  delete outputFile;

  MACH3LOG_INFO("Finished dumping ParameterHandler object");
}
