#include "covarianceXsec.h"

// ********************************************
// ETA - YAML constructor
// this will replace the root file constructor but let's keep it in
// to do some validations
covarianceXsec::covarianceXsec(const std::vector<std::string>& YAMLFile, const char *name, double threshold, int FirstPCA, int LastPCA)
               : covarianceBase(YAMLFile, name, threshold, FirstPCA, LastPCA){
// ********************************************

  InitXsecFromConfig();
  SetupNormPars();

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

  _fDetID = std::vector<int>(_fNumPar);
  _fParamType = std::vector<SystType>(_fNumPar);
  //_fDetString = std::vector<std::string>(_fNumPar);
  isFlux.resize(_fNumPar);
  //Vector of vectors of strings to contain potentially multiple variables that
  //might be cut on
  _fKinematicPars = std::vector<std::vector<std::string>>(_fNumPar);
  //Vector of vector of ints to contain the lower and upper bounds of a cut
  //for a particular kinematic variables
  _fKinematicBounds = std::vector<std::vector<std::vector<double>>>(_fNumPar);

  //KS: We know at most how params we expect so reserve memory for max possible params. Later we will shrink to size to not waste memory. Reserving means slightly faster loading and possible less memory fragmentation.
  _fSplineInterpolationType.reserve(_fNumPar);
  _fTargetNuclei.reserve(_fNumPar);
  _fNeutrinoFlavour.reserve(_fNumPar);
  _fNeutrinoFlavourUnosc.reserve(_fNumPar);
  _fNormModes.reserve(_fNumPar);
  _fSplineKnotUpBound.reserve(_fNumPar);
  _fSplineKnotLowBound.reserve(_fNumPar);

  int i = 0;

  //ETA - read in the systematics. Would be good to add in some checks to make sure
  //that there are the correct number of entries i.e. are the _fNumPars for Names,
  //PreFitValues etc etc.
  for (auto const &param : _fYAMLDoc["Systematics"])
  {
    _fDetID[i] = (param["Systematic"]["DetID"].as<int>());

	 //Fill the map to get the correlations later as well
	 std::string ParamType = param["Systematic"]["Type"].as<std::string>();
	 int nFDSplines;
	 //Now load in variables for spline systematics only
	 if (ParamType.find(SystType_ToString(SystType::kSpline)) != std::string::npos)
     {
       _fParamType[i] = SystType::kSpline;
	   if (param["Systematic"]["SplineInformation"]["FDSplineName"]) {
		 _fFDSplineNames.push_back(param["Systematic"]["SplineInformation"]["FDSplineName"].as<std::string>());
		 nFDSplines++;
	   }

	   if (param["Systematic"]["SplineInformation"]["FDMode"]) {
		 _fFDSplineModes.push_back(param["Systematic"]["SplineInformation"]["FDMode"].as<std::vector<int>>());
	   }

      if (param["Systematic"]["SplineInformation"]["NDSplineName"]) {
        _fNDSplineNames.push_back(param["Systematic"]["SplineInformation"]["NDSplineName"].as<std::string>());
      }

      //Now get the Spline interpolation type
      if (param["Systematic"]["SplineInformation"]["InterpolationType"]){
        for(int InterpType = 0; InterpType < kSplineInterpolations ; ++InterpType){
          if(param["Systematic"]["SplineInformation"]["InterpolationType"].as<std::string>() == SplineInterpolation_ToString(SplineInterpolation(InterpType)))
            _fSplineInterpolationType.push_back(SplineInterpolation(InterpType));
        }
      } else { //KS: By default use TSpline3
        _fSplineInterpolationType.push_back(SplineInterpolation(kTSpline3));
      }

      _fSplineKnotUpBound.push_back(GetFromManager<double>(param["Systematic"]["SplineInformation"]["SplineKnotUpBound"], 9999));
      _fSplineKnotLowBound.push_back(GetFromManager<double>(param["Systematic"]["SplineInformation"]["SplineKnotLowBound"], -9999));

      } else if(param["Systematic"]["Type"].as<std::string>() == SystType_ToString(SystType::kNorm)) {
       _fParamType[i] = SystType::kNorm;
        // ETA Empty DummyVector can be used to specify no cut for mode, target and neutrino flavour
       // ETA Has to be of size 0 to mean apply to all
        std::vector<int> DummyModeVec;
        //Ultimately all this information ends up in the NormParams vector

        // Set the target of the normalisation parameter
       _fTargetNuclei.push_back(GetFromManager<std::vector<int>>(param["Systematic"]["TargetNuclei"], DummyModeVec));
       _fNeutrinoFlavour.push_back(GetFromManager<std::vector<int>>(param["Systematic"]["NeutrinoFlavour"], DummyModeVec));
       _fNeutrinoFlavourUnosc.push_back(GetFromManager<std::vector<int>>(param["Systematic"]["NeutrinoFlavourUnosc"], DummyModeVec));
       _fNormModes.push_back(GetFromManager<std::vector<int>>(param["Systematic"]["Mode"], DummyModeVec));
      }
     else if(param["Systematic"]["Type"].as<std::string>() == SystType_ToString(SystType::kFunc)){
       _fParamType[i] = SystType::kFunc;
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
       throw;
      }

	 //ETA - I think this can go in the norm parameters only if statement above
	 int NumKinematicCuts = 0;
	 if(param["Systematic"]["KinematicCuts"]){

	   NumKinematicCuts = param["Systematic"]["KinematicCuts"].size();

	   std::vector<std::string> TempKinematicStrings;
	   std::vector<std::vector<double>> TempKinematicBounds;
	   //First element of TempKinematicBounds is always -999, and size is then 3
	   for(int KinVar_i = 0 ; KinVar_i < NumKinematicCuts ; ++KinVar_i){ 
		 //ETA
		 //This is a bit messy, Kinematic cuts is a list of maps
		 //The for loop h
		 for (YAML::const_iterator it = param["Systematic"]["KinematicCuts"][KinVar_i].begin();it!=param["Systematic"]["KinematicCuts"][KinVar_i].end();++it) {
		   TempKinematicStrings.push_back(it->first.as<std::string>());
		   TempKinematicBounds.push_back(it->second.as<std::vector<double>>());
		   std::vector<double> bounds = it->second.as<std::vector<double>>();
		 }
	   }
        _fKinematicPars.at(i) = TempKinematicStrings;
        _fKinematicBounds.at(i) = TempKinematicBounds;
      }
	 if(_fFancyNames[i].find("b_")==0) isFlux[i] = true;
	 else isFlux[i] = false;
	 i++;
  }
  _fSplineInterpolationType.shrink_to_fit();
  _fTargetNuclei.shrink_to_fit();
  _fNeutrinoFlavour.shrink_to_fit();
  _fNeutrinoFlavourUnosc.shrink_to_fit();
  _fNormModes.shrink_to_fit();
  _fSplineKnotUpBound.shrink_to_fit();
  _fSplineKnotLowBound.shrink_to_fit();
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
  returnVec.reserve(_fNumPar);
  int nd_counter = 0;
  //int counter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
      if (GetParamType(i) == kSpline) { //If parameter is implemented as a spline

        //ETA - Horrible hard-code becuase ND and FD spline names are different...
        // this is only true at T2K and needs to change!
        if(DetID == 1){
          returnVec.push_back(_fNDSplineNames[nd_counter]);
          nd_counter++;
        }
        else{
          returnVec.push_back(GetParName(i));
        }
      }
    }
  }
  returnVec.shrink_to_fit();

  return returnVec;
}

// ********************************************
// ETA - this is a complete fudge for now and is only here because on
// T2K the splines at ND280 and SK have different names... 
// this is completely temporary and will be removed in the future
const std::vector<std::string> covarianceXsec::GetFDSplineFileParsNamesFromDetID(const int DetID) {
// ********************************************
  std::vector<std::string> returnVec;
  int FDSplineCounter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
	//std::cout << " Param i has DetID " << GetParDetID(i) << " from yaml" << std::endl;
	if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
	  if (GetParamType(i) == kSpline) { //If parameter is implemented as a spline
		//This is horrible but has to be here whilst the ND and FD splines are named
		//differently on T2K. This is easy to fix but isn't currently available.
		returnVec.push_back(_fFDSplineNames[FDSplineCounter]); //Append spline name
		FDSplineCounter++;
	  }
	}
  }
  return returnVec;
}

// ********************************************
// ETA - this is another fudge for now and is only here because on
// T2K the splines at ND280 and SK have different names... 
// this is completely temporary and will be removed in the future
const std::vector<std::string> covarianceXsec::GetNDSplineFileParsNamesFromDetID(const int DetID) {
// ********************************************
  std::vector<std::string> returnVec;
  int NDSplineCounter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
	if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
	  if (GetParamType(i) == kSpline) { //If parameter is implemented as a spline
		//This is horrible but has to be here whilst the ND and ND splines are named
		//differently on T2K. This is easy to fix but isn't currently available.
		std::cout << "Getting ND spline name " << _fNDSplineNames[NDSplineCounter] << std::endl;
		returnVec.push_back(_fNDSplineNames[NDSplineCounter]); //Append spline name
	  }
	  NDSplineCounter++;
	}
  }
  return returnVec;
}

// ********************************************
// DB Grab the Spline Names for the relevant DetID
const std::vector<std::string> covarianceXsec::GetSplineFileParsNamesFromDetID(const int DetID) {
// ********************************************
  std::vector<std::string> returnVec;

  int FDSplineCounter = 0;
  int NDSplineCounter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
      if (GetParamType(i) == kSpline) { //If parameter is implemented as a spline
		//This is horrible but has to be here whilst the ND and FD splines are named
		//differently on T2K. This is easy to fix but isn't currently available.
		if((GetParDetID(i) & 1) == 1){
		  returnVec.push_back(_fNDSplineNames[NDSplineCounter]);
		  NDSplineCounter++;
		} else{
		  returnVec.push_back(_fFDSplineNames[FDSplineCounter]);
		  FDSplineCounter++;
		}
      }
    }
  }

  return returnVec;
}

// ********************************************
// DB Grab the Spline Modes for the relevant DetID
const std::vector< std::vector<int> > covarianceXsec::GetSplineModeVecFromDetID(const int DetID) {
// ********************************************
  std::vector< std::vector<int> > returnVec;

  //Need a counter or something to correctly get the index in _fFDSplineModes since it's not of length nPars
  //Should probably just make a std::map<std::string, int> for param name to FD spline index
  int nFDSplineCounter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
	if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
	  if (GetParamType(i) == kSpline) { //If parameter is implemented as a spline
		returnVec.push_back(_fFDSplineModes[nFDSplineCounter]);	
		nFDSplineCounter++;
	  }
	}
  }
  return returnVec;
}

// ********************************************
// Get Norm params
XsecNorms4 covarianceXsec::GetXsecNorm(const int i, const int norm_counter) {
// ********************************************
  XsecNorms4 norm;
  norm.name = GetParFancyName(i);

  //Copy the mode information into an XsecNorms4 struct
  norm.modes = _fNormModes[norm_counter];
  norm.pdgs = _fNeutrinoFlavour[norm_counter];
  norm.preoscpdgs = _fNeutrinoFlavourUnosc[norm_counter];
  norm.targets = _fTargetNuclei[norm_counter];

  //Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
  //We set a bool to see if any bounds exist so we can short-circuit checking all of them every step
  bool HasKinBounds = false;

  if(_fKinematicPars.at(i).size() > 0) HasKinBounds = true;

  for(unsigned int KinematicCut_i = 0; KinematicCut_i < _fKinematicPars[i].size(); ++KinematicCut_i) {
    //Push back with the string for the kinematic cut
    norm.KinematicVarStr.push_back(_fKinematicPars.at(i).at(KinematicCut_i));
    //Push back with the bounds for the kinematic cut
    norm.Selection.push_back(_fKinematicBounds.at(i).at(KinematicCut_i));
  }
  norm.hasKinBounds = HasKinBounds;
  //End of kinematic bound checking

  // Set the global parameter index of the normalisation parameter
  norm.index = i;

  return norm;
}

// ********************************************
// DB Grab the Normalisation parameters
// I have changed this because it is quite nice for the covariance object not to care
// which samplePDF a parameter should affect or not.
void covarianceXsec::SetupNormPars(){
// ********************************************
  //ETA - in case NormParams already is filled
  NormParams.clear();

  int norm_counter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
    if (GetParamType(i) == kNorm) { //If parameter is implemented as a normalisation
      //Add this parameter to the vector of parameters
      NormParams.push_back(GetXsecNorm(i, norm_counter));
      norm_counter++;
    }
  }
  return; 
}

// ********************************************
// DB Grab the Normalisation parameters for the relevant DetID
const std::vector<XsecNorms4> covarianceXsec::GetNormParsFromDetID(const int DetID) {
// ********************************************
  std::vector<XsecNorms4> returnVec;
  int norm_counter = 0;

  for (int i = 0; i < _fNumPar; ++i) {
    if (GetParamType(i) == kNorm) { //If parameter is implemented as a normalisation
      if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
        //Add this parameter to the vector of parameters
          returnVec.push_back(GetXsecNorm(i, norm_counter));
      }
      norm_counter++;
    }
  }
  return returnVec;
}


// ********************************************
// DB Grab the number of parameters for the relevant DetID
int covarianceXsec::GetNumParamsFromDetID(const int DetID, const SystType Type) {
// ********************************************
  int returnVal = 0;

  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
      if (GetParamType(i) == Type) { //If parameter is implemented as a functional parameter
        returnVal += 1;
      }
    }
  }
  return returnVal;
}

// ********************************************
// DB Grab the parameter names for the relevant DetID
const std::vector<std::string> covarianceXsec::GetParsNamesFromDetID(const int DetID, const SystType Type) {
// ********************************************
  std::vector<std::string> returnVec;

  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
      if (GetParamType(i) == Type) { //If parameter is implemented as a functional param
        returnVec.push_back(GetParFancyName(i));
      }
    }
  }
  return returnVec;
}

// ********************************************
// DB DB Grab the parameter indices for the relevant DetID
const std::vector<int> covarianceXsec::GetParsIndexFromDetID(const int DetID, const SystType Type) {
// ********************************************
  std::vector<int> returnVec;

  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
      if (GetParamType(i) == Type) { //If parameter is implemented as a functional param
        returnVec.push_back(i);
      }
    }
  }
  return returnVec;
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
  CorrelateSteps();
}

// ********************************************
// Print everything we know about the inputs we're Getting
void covarianceXsec::Print() {
// ********************************************
  MACH3LOG_INFO("#################################################");
  MACH3LOG_INFO("Printing covarianceXsec:");

  MACH3LOG_INFO("============================================================================================================================================================");
  MACH3LOG_INFO("{:<5} {:2} {:<40} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<5} {:2} {:<10}", "#", "|", "Name", "|", "Nom.", "|", "Prior", "|", "Error", "|", "Lower", "|", "Upper", "|", "StepScale", "|", "DetID", "|", "Type");
  MACH3LOG_INFO("------------------------------------------------------------------------------------------------------------------------------------------------------------");
  for (int i = 0; i < GetNumParams(); i++) {
    std::string ErrString = fmt::format("{:.2f}", _fError[i]);
    MACH3LOG_INFO("{:<5} {:2} {:<40} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<10} {:2} {:<5} {:2} {:<10}", i, "|", GetParFancyName(i), "|", _fGenerated[i], "|", _fPreFitValue[i], "|", "+/- " + ErrString, "|", _fLowBound[i], "|", _fUpBound[i], "|", _fIndivStepScale[i], "|", _fDetID[i], "|", SystType_ToString(_fParamType[i]));
  }
  MACH3LOG_INFO("============================================================================================================================================================");

  // Output the normalisation parameters as a sanity check!
  MACH3LOG_INFO("Normalisation parameters:  {}", NormParams.size());

  //KS: Consider making some class producing table..
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┬───────────────┬───────────────┬───────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│{3:15}│{4:15}│{5:15}│", "#", "Global #", "Name", "Int. mode", "Target", "pdg");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┼───────────────┼───────────────┼───────────────┤");

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

    MACH3LOG_INFO("│{: <4}│{: <10}│{: <40}│{: <15}│{: <15}│{: <15}│", i, NormParams[i].index, NormParams[i].name, intModeString, targetString, pdgString);
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┴───────────────┴───────────────┴───────────────┘");

  std::vector<int> SplineParsIndex;
  for (int i = 0; i < _fNumPar; ++i)
  {
    if (GetParamType(i) == kSpline) { SplineParsIndex.push_back(i); }
  }

  MACH3LOG_INFO("Spline parameters: {}", SplineParsIndex.size());

  MACH3LOG_INFO("=========================================================================================================================");
  MACH3LOG_INFO("{:<4} {:<2} {:<40} {:<2} {:<20} {:<2} {:<20} {:<2} {:<20} {:<2}", "#", "|", "Name", "|", "Spline Interpolation", "|", "Low Knot Bound", "|", "Up Knot Bound", "|");
  MACH3LOG_INFO("-------------------------------------------------------------------------------------------------------------------------");
  for (unsigned int i = 0; i < SplineParsIndex.size(); ++i) {
    MACH3LOG_INFO("{:<4} {:<2} {:<40} {:<2} {:<20} {:<2} {:<20} {:<2} {:<20} {:<2}", i, "|", GetParFancyName(SplineParsIndex[i]), "|",
                  SplineInterpolation_ToString(_fSplineInterpolationType[i]), "|",
                  _fSplineKnotLowBound[i], "|", _fSplineKnotUpBound[i], "|");
  }
  MACH3LOG_INFO("=========================================================================================================================");

  std::vector<int> FuncParsIndex;
  for (int i = 0; i < _fNumPar; ++i)
  {
    if (GetParamType(i) == kFunc) { FuncParsIndex.push_back(i); }
  }

  MACH3LOG_INFO("Functional parameters: {}", FuncParsIndex.size());
  MACH3LOG_INFO("┌────┬──────────┬────────────────────────────────────────┐");
  MACH3LOG_INFO("│{0:4}│{1:10}│{2:40}│", "#", "Global #", "Name");
  MACH3LOG_INFO("├────┼──────────┼────────────────────────────────────────┤");
  for (unsigned int i = 0; i < FuncParsIndex.size(); ++i) {
      MACH3LOG_INFO("│{0:4}│{1:<10}│{2:40}│", std::to_string(i), FuncParsIndex[i], GetParFancyName(FuncParsIndex[i]));
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────────────────────────┘");

  CheckCorrectInitialisation();
} // End

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
  CheckForDuplicates(_fNDSplineNames, "_fNDSplineNames");
  CheckForDuplicates(_fFDSplineNames, "_fFDSplineNames");
}

// ********************************************
// Sets the proposed Flux parameters to the prior values
void covarianceXsec::setFluxOnlyParameters() {
// ********************************************
  if(!pca) {
    for (int i = 0; i < _fNumPar; i++) {
      if(isFlux[i]) _fPropVal[i] = _fPreFitValue[i];
    }
  } else {
    MACH3LOG_ERROR("setFluxOnlyParameters not implemented for PCA");
    throw;
  }
}

// ********************************************
// Sets the proposed Flux parameters to the prior values
void covarianceXsec::setXsecOnlyParameters() {
// ********************************************
  if(!pca) {
    for (int i = 0; i < _fNumPar; i++) {
      if(!isFlux[i]) _fPropVal[i] = _fPreFitValue[i];
    }
  } else {
    MACH3LOG_ERROR("setXsecOnlyParameters not implemented for PCA");
    throw;
  }
}
