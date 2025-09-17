#include "Parameters/ParameterTunes.h"

// ********************************************
ParameterTunes::ParameterTunes(const YAML::Node& Settings) {
// ********************************************
  std::map<std::string, std::vector<double>> orderedTunes;

  // Loop over all systematics
  for (const auto& sysNode : Settings) {
    const auto& values = sysNode["Systematic"]["ParameterValues"];

    for (const auto& tuneNode : values) {
      std::string key = tuneNode.first.as<std::string>();
      if (key == "PreFitValue") continue;

      double val = tuneNode.second.as<double>();
      orderedTunes[key].push_back(val);
    }
  }

  size_t expectedSize = 0;
  bool hasMismatch = false;
  std::ostringstream errorMsg;
  errorMsg << "Inconsistent number of parameter values across tunes:\n";

  for (const auto& kv : orderedTunes) {
    const auto& vals = kv.second;
    if (expectedSize == 0) {
      expectedSize = vals.size();
    } else if (vals.size() != expectedSize) {
      hasMismatch = true;
    }
    errorMsg << " - Tune '" << kv.first << "': " << vals.size() << " values\n";
  }

  if (hasMismatch) {
    errorMsg << "Expected all tunes to have " << expectedSize << " values.";
    MACH3LOG_ERROR("{}", errorMsg.str());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  int idx = 0;
  for (const auto& kv : orderedTunes) {
    TuneNames.push_back(kv.first);
    TuneValues.push_back(kv.second);
    TuneMap[kv.first] = idx++;
  }

  // KS: Log final summary of found tunes
  MACH3LOG_INFO("Found {} tunes:", TuneNames.size());
  for (size_t i = 0; i < TuneNames.size(); ++i) {
    MACH3LOG_INFO("  Tune {} {}", i, TuneNames[i]);
  }
  #ifdef DEBUG
  PrintTunes();
  #endif
}

// ********************************************
ParameterTunes::~ParameterTunes() {
// ********************************************

}

// ********************************************
std::vector<double> ParameterTunes::GetTune(const int TuneNumber) const {
// ********************************************
  if (TuneNumber >= static_cast<int>(TuneValues.size())) {
    MACH3LOG_ERROR("I only have {} tune(s), but you asked for index {}", TuneValues.size(), TuneNumber);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return TuneValues[TuneNumber];
}

// ********************************************
std::vector<double> ParameterTunes::GetTune(const std::string& TuneName) const {
// ********************************************
  auto it = TuneMap.find(TuneName);
  if (it == TuneMap.end()) {
    MACH3LOG_ERROR("Tune '{}' not found. Available tunes: {}", TuneName, fmt::join(TuneNames, ", "));
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return GetTune(it->second);
}

// ********************************************
void ParameterTunes::PrintTunes() const {
// ********************************************
  MACH3LOG_INFO("Available tunes:");
  for (size_t i = 0; i < TuneNames.size(); ++i) {
    MACH3LOG_INFO("  Tune {} {}", i, TuneNames[i]);
      for (size_t j = 0; j < TuneValues[i].size(); ++j) {
        MACH3LOG_INFO("    Value {} = {}", j, TuneValues[i][j]);
    }
  }
}
