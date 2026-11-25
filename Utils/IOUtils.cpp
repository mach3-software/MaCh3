#include "Utils/IOUtils.h"

_MaCh3_Safe_Include_Start_ //{
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RConfig.hxx>
_MaCh3_Safe_Include_End_ //}

namespace M3 {
// **************************************************************************
TFile* Open(const std::string& Name, const std::string& Type, const std::string& File, const int Line) {
// **************************************************************************
  TFile* OutFile = new TFile(Name.c_str(), Type.c_str());

  // Check if the file is successfully opened and usable
  if (OutFile->IsZombie()) {
    MACH3LOG_ERROR("Failed to open file: {}", Name);
    std::string lowerType = Type;
    std::transform(lowerType.begin(), lowerType.end(), lowerType.begin(), ::tolower);
    if (lowerType == "recreate") {
      MACH3LOG_ERROR("Check if directory exist");
    }
    delete OutFile;
    throw MaCh3Exception(File, Line);
  }
  return OutFile;
}

// ***************************************************************************
void EnableROOTMultiThreading() {
// ***************************************************************************
  #ifdef MULTITHREAD
  auto nThreads = GetNThreads();
  const bool IMTWasEnabled = ROOT::IsImplicitMTEnabled();
  if(!IMTWasEnabled) ROOT::EnableImplicitMT(nThreads);
  #endif
}

// ***************************************************************************
void DisableROOTMultiThreading() {
// ***************************************************************************
  const bool IMTWasEnabled = ROOT::IsImplicitMTEnabled();
  if(!IMTWasEnabled) ROOT::DisableImplicitMT();
}

} //end M3 namespace
