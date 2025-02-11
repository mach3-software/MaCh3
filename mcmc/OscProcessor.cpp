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
void OscProcessor::ReadOSCFile() {
// ***************
  // Call base class constructor
  MCMCProcessor::ReadOSCFile();

  if(PlotJarlskog)
  {
    Chain->SetAlias("J_cp", "TMath::Sqrt(sin2th_13)*TMath::Sqrt(1.-sin2th_13)*TMath::Sqrt(1.-sin2th_13)*TMath::Sqrt(sin2th_12)*TMath::Sqrt(1.-sin2th_12)*TMath::Sqrt(sin2th_23)*TMath::Sqrt(1.-sin2th_23)*TMath::Sin(delta_cp)");
    BranchNames.push_back("J_cp");
    ParamType.push_back(kOSCPar);
    nParam[kOSCPar]++;
    nDraw++;

    /// @todo we should actually calculate central value and prior error but leave it for now...
    ParamNom[kOSCPar].push_back( 0. );
    ParamCentral[kOSCPar].push_back( 0. );
    ParamErrors[kOSCPar].push_back( 1. );
    // Push back the name
    ParamNames[kOSCPar].push_back("J_cp");
      ParamFlat[kOSCPar].push_back( false );
  }
}
