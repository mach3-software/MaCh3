//MaCh3 includes
#include "MaCh3Modes.h"

// *******************
MaCh3Modes::MaCh3Modes(std::string const &filename) {
// *******************

  // Load config
  YAML::Node config = YAML::LoadFile(filename);

  std::string GetMaCh3ModeName(const int Index);
  NModes = 0;

  Title = config["Title"].as<std::string>();
  Generator = config["GeneratorName"].as<std::string>();

  std::vector<std::string> names = config["MaCh3Modes"].as<std::vector<std::string>>();

  for(size_t i = 0; i < names.size(); i++)
  {
    DeclareNewMode(names[i],
                    config[names[i]]["Name"].as<std::string>(),
                    config[names[i]]["GeneratorMaping"].as<std::vector<int>>());
  }
  // Add unknown category, it's better to have garbage category where all undefined modes will go rather than get random crashes
  DeclareNewMode("UNKNOWN_BAD",
                 "UNKNOWN_BAD",
                 {});
  // This is hack to not have bad mode
  NModes--;

  PrepareMap();

  Print();
}

// *******************
void MaCh3Modes::Print() {
// *******************

  MACH3LOG_INFO("Printing MaCh3 Modes Called: {}", Title);

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("{:<5} {:2} {:<20} {:2} {:<20} {:2} {:<30}", "#", "|", "Name", "|", "FancyName", "|", std::string(Generator+" Modes"));
  MACH3LOG_INFO("------------------------------------------------------------------------");
  for(int i = 0; i < NModes; ++i) {
    auto Name = fMode[i].Name;
    auto FancyName = fMode[i].FancyName;
    auto Values = fMode[i].GeneratorMaping;

    std::string generatorModes;
    for (const int& element : Values) {
      generatorModes += std::to_string(element) + " ";
    }
    MACH3LOG_INFO("{:<5} {:2} {:<20} {:2} {:<20} {:2} {:<30}", i, "|", Name, "|", FancyName, "|", generatorModes);
  }
  MACH3LOG_INFO("========================================================================");

  MACH3LOG_INFO("==========================");
  MACH3LOG_INFO("{:<10} {:2} {:<30}", std::string(Generator + " Modes"), "|", "Name");
  MACH3LOG_INFO("--------------------------");
  for (size_t i = 0; i < ModeMap.size(); ++i) {
    MACH3LOG_INFO("{:<10} {:2} {:<30}", i, "|", fMode[ModeMap[i]].Name);
  }
  MACH3LOG_INFO("==========================");
}

// *******************
MaCh3Modes_t MaCh3Modes::EnsureModeNameRegistered(std::string const &name) {
// *******************
  if (Mode.count(name)) {
    return Mode[name];
  }
  MaCh3Modes_t index = Mode.size();
  Mode[name] = index;
  return index;
}

// *******************
void MaCh3Modes::DeclareNewMode(std::string const &name,
                    std::string const &fancyname,
                    std::vector<int> const &GenMap) {
// *******************
  MaCh3ModeInfo newinfo;
  newinfo.GeneratorMaping = GenMap;
  newinfo.FancyName = fancyname;
  newinfo.Name = name;

  MaCh3Modes_t index = EnsureModeNameRegistered(name);

  fMode.emplace(index, newinfo);
  NModes++;
}

// *******************
void MaCh3Modes::PrepareMap() {
// *******************

  int maxElement = 0;
  for(int i = 0; i < NModes; ++i) {
    for(size_t j = 0; j < fMode[i].GeneratorMaping.size(); j++) {
      maxElement = std::max(maxElement, fMode[i].GeneratorMaping[j]);
    }
  }
  ModeMap.resize(maxElement+1, Mode["UNKNOWN_BAD"]);

  for(int m = 0; m < NModes; ++m) {
    for(size_t i = 0; i < fMode[m].GeneratorMaping.size(); ++i) {
      if(fMode[m].GeneratorMaping[i] < 0) {
        MACH3LOG_ERROR("Negative value of Mode {} for mode {}", fMode[m].GeneratorMaping[i], fMode[m].Name);
        throw;
      }
      if(ModeMap[fMode[m].GeneratorMaping[i]] != Mode["UNKNOWN_BAD"])
      {
        MACH3LOG_ERROR("Generator mode {} already defined in mode {} can't do this for mode {}", fMode[m].GeneratorMaping[i], fMode[ModeMap[fMode[m].GeneratorMaping[i]]].Name, fMode[m].Name);
        throw;
      }
      ModeMap[fMode[m].GeneratorMaping[i]] = m;
    }
  }
}

// *******************
std::string MaCh3Modes::GetMaCh3ModeName(const int Index) {
// *******************
  if(Index < 0)
    MACH3LOG_CRITICAL("Mode you look for is smaller than 0 and equal to {}", Index);

  // return UNKNOWN_BAD if out of boundary
  if(Index > NModes)
  {
    MACH3LOG_DEBUG("Asking for mode {}, while I only have {}, returning {} mode", Index, NModes, fMode[NModes].Name);
    return fMode[NModes].Name;
  }
  return fMode[Index].Name;
}

// *******************
std::string MaCh3Modes::GetMaCh3ModeFancyName(const int Index) {
// *******************
  // return UNKNOWN_BAD if out of boundary
  if(Index < 0)
    MACH3LOG_CRITICAL("Mode you look for is smaller than 0 and equal to {}", Index);

  if(Index > NModes)
  {
    MACH3LOG_DEBUG("Asking for mode {}, while I only have {}, returning {} mode", Index, NModes, fMode[NModes].Name);
    return fMode[NModes].FancyName;
  }
  return fMode[Index].FancyName;
}

// *******************
MaCh3Modes_t MaCh3Modes::GetMode(std::string name) {
// *******************
  if (Mode.count(name)) {
    return Mode[name];
  }
  MACH3LOG_DEBUG("Asking for mode {}, while I only have {}, returning {} mode", Index, NModes, fMode[NModes].Name);

  // return UNKNOWN_BAD
  return NModes;
}

// *******************
MaCh3Modes_t MaCh3Modes::GetModeFromGenerator(const int Index) {
// *******************
  if(Index < 0)
    MACH3LOG_CRITICAL("Mode you look for is smaller than 0 and equal to {}", Index);

  if(Index > NModes)
  {
    MACH3LOG_DEBUG("Asking for mode {}, while I only have {}, returning {} mode", Index, NModes, fMode[NModes].Name);
    return NModes;
  }
  return ModeMap[Index];
}
