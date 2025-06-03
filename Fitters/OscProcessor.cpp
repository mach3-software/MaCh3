#include "OscProcessor.h"

// ****************************
OscProcessor::OscProcessor(const std::string &InputFile) : MCMCProcessor(InputFile) {
// ****************************
  //KS: WARNING this only work when you project from Chain, will nor work when you try SetBranchAddress etc. Turn it on only if you know how to use it
  PlotJarlskog = false;

/// @todo Here where we should add all unitarity triangles, fancy Jarlskog studies and other hacky things that only make sense for oscitations

}

// ****************************
// The destructor
OscProcessor::~OscProcessor() {
// ****************************

}

// ***************
// Read the Osc cov file and get the input central values and errors
void OscProcessor::ReadXSecFile() {
// ***************
  // Call base class function
  MCMCProcessor::ReadXSecFile();

  // KS: Check if OscParams were enabled, in future we will also get
  for(size_t i = 0; i < ParameterGroup.size(); i++) {
    if(ParameterGroup[i] == "Osc"){
      OscEnabled = true;
      break;
    }
  }

  if(PlotJarlskog && OscEnabled)
  {
    Chain->SetAlias("J_cp", "TMath::Sqrt(sin2th_13)*TMath::Sqrt(1.-sin2th_13)*TMath::Sqrt(1.-sin2th_13)*TMath::Sqrt(sin2th_12)*TMath::Sqrt(1.-sin2th_12)*TMath::Sqrt(sin2th_23)*TMath::Sqrt(1.-sin2th_23)*TMath::Sin(delta_cp)");
    BranchNames.push_back("J_cp");
    ParamType.push_back(kXSecPar);
    nParam[kXSecPar]++;
    nDraw++;

    /// @todo we should actually calculate central value and prior error but leave it for now...
    ParamNom[kXSecPar].push_back( 0. );
    ParamCentral[kXSecPar].push_back( 0. );
    ParamErrors[kXSecPar].push_back( 1. );
    // Push back the name
    ParamNames[kXSecPar].push_back("J_cp");
    ParamFlat[kXSecPar].push_back( false );
  } else if(PlotJarlskog && !OscEnabled) {
    MACH3LOG_ERROR("Trying to enable Jarlskog without oscillations");
    throw MaCh3Exception(__FILE__,__LINE__);
  }
}
