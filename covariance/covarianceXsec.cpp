#include "covarianceXsec.h"

// ********************************************
// ETA - YAML constructor
// this will replace the root file constructor but let's keep it in
// to do some validations
covarianceXsec::covarianceXsec(std::vector<std::string> YAMLFile, double threshold, int FirstPCAdpar, int LastPCAdpar)
               : covarianceBase(YAMLFile, threshold, FirstPCAdpar, LastPCAdpar){
// ********************************************

  setName("xsec_cov");
  InitXsecFromConfig();
  SetupNormPars();

  //ETA - again this really doesn't need to be hear...
  for (int i = 0; i < _fNumPar; i++)
  {
    // Sort out the print length
    if(_fNames[i].length() > PrintLength) PrintLength = _fNames[i].length();
  } // end the for loop

  MACH3LOG_INFO("Constructing instance of covarianceXsec");
  initParams(0.001);
  // Print
  Print();
}

// ********************************************
void covarianceXsec::InitXsecFromConfig() {
// ********************************************

  _fDetID = std::vector<int>(_fNumPar);
  _fParamType = std::vector<std::string>(_fNumPar);
  //_fDetString = std::vector<std::string>(_fNumPar);
  isFlux.resize(_fNumPar);
  //Vector of vectors of strings to contain potentially multiple variables that
  //might be cut on
  _fKinematicPars = std::vector<std::vector<std::string>>(_fNumPar);
  //Vector of vector of ints to contain the lower and upper bounds of a cut
  //for a particular kinematic variables
  _fKinematicBounds = std::vector<std::vector<std::vector<double>>>(_fNumPar);

  int i = 0;

  //ETA - read in the systematics. Would be good to add in some checks to make sure
  //that there are the correct number of entries i.e. are the _fNumPars for Names,
  //PreFitValues etc etc.
  for (auto const &param : _fYAMLDoc["Systematics"])
  {
    _fParamType[i] = (param["Systematic"]["Type"].as<std::string>());
    _fDetID[i] = (param["Systematic"]["DetID"].as<int>());

	 //Fill the map to get the correlations later as well
	 std::string ParamType = param["Systematic"]["Type"].as<std::string>();
	 //Now load in varaibles for spline systematics only
	 int SplineCounter = 0;
	 if (ParamType.find("Spline") != std::string::npos) {
	   //ETA - do some checks whether various Spline-related information exists
	   if(CheckNodeExists(param["Systematic"], "SplineInformation")){ 

		 if (param["Systematic"]["SplineInformation"]["SplineName"]) {
		   _fSplineNames.push_back(param["Systematic"]["SplineInformation"]["SplineName"].as<std::string>());
		 }

		 _fSplineModes.push_back(GetFromManager(param["Systematic"]["SplineInformation"]["Mode"], std::vector<int>()));

		 if (param["Systematic"]["SplineInformation"]["InterpolationType"]){

		   for(int InterpType = 0 ; InterpType < kSplineInterpolations ; ++InterpType){
			 if(param["Systematic"]["SplineInformation"]["InterpolationType"].as<std::string>() == SplineInterpolation_ToString(SplineInterpolation(InterpType)))
			 {
			   _fSplineInterpolationType.push_back(SplineInterpolation(InterpType));
			 }
		   }
		 }else{
		   //KS: By default use TSpline3
		   _fSplineInterpolationType.push_back(SplineInterpolation(kTSpline3));
		 }
	   }
	   else{
         MACH3LOG_ERROR("Spline Information does not exist for spline {}", param["Systematic"]["FancyName"].as<std::string>());
	   }

	   //Insert the mapping from the spline index i.e. the length of _fSplineNames etc
	   //to the Systematic index i.e. the counter for things like _fDetID and _fDetID
	   _fSplineToSystIndexMap.insert(std::pair{SplineCounter, i}); 
       SplineCounter++;
	 } else if(param["Systematic"]["Type"].as<std::string>() == "Norm") {
 
	   //Empty DummyVector can be used to specify no cut for mode, target and neutrino flavour
	   std::vector<int> DummyModeVec;
	   //Ultimately all thsi information ends up in the NormParams vector

	   // Set the target of the normalisation parameter
	   if(param["Systematic"]["TargetNuclei"]){
		 _fTargetNuclei.push_back(param["Systematic"]["TargetNuclei"].as<std::vector<int>>());
	   } else{
		 //Has to be of size 0 to mean apply to all
		 _fTargetNuclei.push_back(DummyModeVec);
	   }

	   // Set the neutrino flavours to apply to 
	   if(param["Systematic"]["NeutrinoFlavour"]){
		 _fNeutrinoFlavour.push_back(param["Systematic"]["NeutrinoFlavour"].as<std::vector<int>>());
	   } else{
		 //Has to be of size 0 to mean apply to all
		 _fNeutrinoFlavour.push_back(DummyModeVec);
	   }

	   // Set the unoscillated neutrino flavours to apply to which is often used for flux systs 
	   if(param["Systematic"]["NeutrinoFlavourUnosc"]){
		 _fNeutrinoFlavourUnosc.push_back(param["Systematic"]["NeutrinoFlavourUnosc"].as<std::vector<int>>());
	   } else{
		 //Has to be of size 0 to mean apply to all
		 _fNeutrinoFlavourUnosc.push_back(DummyModeVec);
	   }

	   //First check to see if we have specified a mode
	   //std::cout << "Found a norm parameter at " << i << std::endl;
	   if(param["Systematic"]["Mode"]){
		 _fNormModes.push_back(param["Systematic"]["Mode"].as<std::vector<int>>());
	   } else{
		 //Has to be of size 0 to mean apply to all
		 _fNormModes.push_back(DummyModeVec);
	   }
	 }
     else if(param["Systematic"]["Type"].as<std::string>() == "Functional"){
	   //std::cout << "Found a functional parameter!!" << std::endl;
	    
	 }
     else{
       MACH3LOG_ERROR("Given unrecognised systematic type: {}", param["Systematic"]["Type"].as<std::string>());
       MACH3LOG_ERROR("Expecting \"Norm\", \"Spline\" or \"Functional\"");
       throw;
	 }

	 //ETA - I think this can go in the norm parameters only if statement above
	 int NumKinematicCuts = 0;
	 if(param["Systematic"]["KinematicCuts"]){

	   NumKinematicCuts = param["Systematic"]["KinematicCuts"].size();
	   //std::cout << "Number of Kinematic cuts is " << NumKinematicCuts << std::endl;

	   std::vector<std::string> TempKinematicStrings;
	   std::vector<std::vector<double>> TempKinematicBounds;
	   //First element of TempKinematicBounds is always -999, and size is then 3

	   for(int KinVar_i = 0 ; KinVar_i < NumKinematicCuts ; ++KinVar_i){ 
		 //ETA
		 //This is a bit messy, Kinematic cuts is a list of maps
		 //The for loop h
		 for (YAML::const_iterator it=param["Systematic"]["KinematicCuts"][KinVar_i].begin();it!=param["Systematic"]["KinematicCuts"][KinVar_i].end();++it) {
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

  return;
}

// ********************************************
covarianceXsec::~covarianceXsec() {
// ********************************************


}

// ********************************************
// DB Grab the Number of splines for the relevant DetID
int covarianceXsec::GetNumSplineParamsFromDetID(int DetID) {
  int returnVal = 0; 

  for (const auto &[SplineIndex, SystIndex] : _fSplineToSystIndexMap){
    if ((GetParDetID(SystIndex) & DetID)) { //If parameter applies to required DetID
		returnVal += 1;
    }
  }
  return returnVal;
}
// ********************************************

// ********************************************
// DB Grab the Spline Names for the relevant DetID
const std::vector<std::string> covarianceXsec::GetSplineParsNamesFromDetID(int DetID) {

  std::vector<std::string> returnVec;
  for (const auto &[SplineIndex, SystIndex] : _fSplineToSystIndexMap){
    if ((GetParDetID(SystIndex) & DetID )){
      returnVec.push_back(_fSplineNames.at(SplineIndex));
	}
  }

  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the Spline Modes for the relevant DetID
// ETA - This funciton will only work for the far detector.
const std::vector< std::vector<int> > covarianceXsec::GetSplineModeVecFromDetID(int DetID) {
  std::vector< std::vector<int> > returnVec;

  //Need a counter or something to correctly get the index in _fSplineModes since it's not of length nPars
  //Should probably just make a std::map<std::string, int> for param name to FD spline index
  for (const auto &[SplineIndex, SystIndex] : _fSplineToSystIndexMap){
	if ((GetParDetID(SystIndex) & DetID)) { //If parameter applies to required DetID
		returnVec.push_back(_fSplineModes.at(SplineIndex));	
	  }	
  }

  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the Spline Indices for the relevant DetID
const std::vector<int> covarianceXsec::GetSplineParsIndexFromDetID(int DetID) {
  std::vector<int> returnVec;

  for (const auto &[SplineIndex, SystIndex] : _fSplineToSystIndexMap){
    if ((GetParDetID(SystIndex) & DetID)) { //If parameter applies to required DetID
		returnVec.push_back(SplineIndex);
      } 
  }
  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the Normalisation parameters
// ETA - this should be the same as GetNormParsFromDetID but not dependent on the DetID
// I have changed this because it is quite nice for the covariance object not to care
// which samplePDF a parameter should affect or not.
void covarianceXsec::SetupNormPars(){

  //ETA - in case NormParams already is filled
  NormParams.clear();

  int norm_counter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
	if (strcmp(GetParamType(i), "Norm") == 0) { //If parameter is implemented as a normalisation
	
		std::vector<int> temp;
		XsecNorms4 norm;
		norm.name = GetParFancyName(i);

		//Copy the mode information into an XsecNorms4 struct
		norm.modes = _fNormModes[norm_counter];	
		norm.pdgs = _fNeutrinoFlavour[norm_counter];
		norm.preoscpdgs = _fNeutrinoFlavourUnosc[norm_counter];
		norm.targets = _fTargetNuclei[norm_counter];
		
		//Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
		//We set a bool to see if any bounds exist so we can short-circuit checking all of them every step
		bool HasKinBounds=false;

		////////////////////
		//New generic cuts things 
		////////////////////
		if(_fKinematicPars.at(i).size() > 0){
		  HasKinBounds = true;
		}

		for(unsigned int KinematicCut_i = 0 ; KinematicCut_i < _fKinematicPars[i].size() ; ++KinematicCut_i){
		  //Push back with the string for the kinematic cut
		  norm.KinematicVarStr.push_back(_fKinematicPars.at(i).at(KinematicCut_i));
		  //Push back with the bounds for the kinematic cut
          norm.Selection.push_back(_fKinematicBounds.at(i).at(KinematicCut_i));
		}
		norm.hasKinBounds=HasKinBounds;
		//End of kinematic bound checking

		// Set the global parameter index of the normalisation parameter
		norm.index=i;
		//Add this parameter to the vector of parameters
		NormParams.push_back(norm);
		norm_counter++;
	  }
	}

  return; 
}
// ********************************************


// ********************************************
// DB Grab the Normalisation parameters for the relevant DetID
// ETA - I think this doesn't need to be the same as scanParameters, haven't we already got this info??
const std::vector<XsecNorms4> covarianceXsec::GetNormParsFromDetID(int DetID) {
  std::vector<XsecNorms4> returnVec;
  int norm_counter = 0;

  for (int i = 0; i < _fNumPar; ++i) {
	if (strcmp(GetParamType(i), "Norm") == 0) { //If parameter is implemented as a normalisation

	  if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID

		std::vector<int> temp;

		XsecNorms4 norm;
		norm.name=GetParFancyName(i);

		//Copy the mode information into an XsecNorms4 struct
		norm.modes = _fNormModes[norm_counter];	
		norm.pdgs = _fNeutrinoFlavour[norm_counter];
		norm.preoscpdgs = _fNeutrinoFlavourUnosc[norm_counter];
		norm.targets = _fTargetNuclei[norm_counter];
		
		//Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
		//We set a bool to see if any bounds exist so we can shortcircuit checking all of them every step
		bool HasKinBounds=false;

		////////////////////
		//New generic cuts things 
		////////////////////
		if(_fKinematicPars.at(i).size() > 0){
		  HasKinBounds = true;
		}

		for(unsigned int KinematicCut_i = 0 ; KinematicCut_i < _fKinematicPars[i].size() ; ++KinematicCut_i){
		  //Push back with the string for the kinematic cut
		  //std::cout << "----------------------" << std::endl;
		  //std::cout << "Will apply a cut on " << _fKinematicPars.at(i).at(KinematicCut_i) << std::endl;
		  norm.KinematicVarStr.push_back(_fKinematicPars.at(i).at(KinematicCut_i));
		  //std::cout << "With bounds " << _fKinematicBounds.at(i).at(KinematicCut_i).at(0) << " to " << _fKinematicBounds.at(i).at(KinematicCut_i).at(1) << std::endl;
		  //Push back with the bounds for the kinematic cut
          norm.Selection.push_back(_fKinematicBounds.at(i).at(KinematicCut_i));
		}

		norm.hasKinBounds=HasKinBounds;
		//End of kinematic bound checking

		// Set the global parameter index of the normalisation parameter
		norm.index=i;
		//Add this parameter to the vector of parameters
		returnVec.push_back(norm);
	  }
	  norm_counter++;
	}
  }

  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the number of Normalisation parameters for the relevant DetID
int covarianceXsec::GetNumFuncParamsFromDetID(int DetID) {
  int returnVal = 0;

  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
      if (strcmp(GetParamType(i), "Functional") == 0) { //If parameter is implemented as a functional parameter
		returnVal += 1;
      }
    }
  }

  return returnVal;
}
// ********************************************

// ********************************************
// DB Grab the Functional parameter names for the relevant DetID
const std::vector<std::string> covarianceXsec::GetFuncParsNamesFromDetID(int DetID) {
  std::vector<std::string> returnVec;

  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
      if (strcmp(GetParamType(i), "Functional") == 0) { //If parameter is implemented as a functional param
		returnVec.push_back(GetParFancyName(i));
      }
    }
  }

  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the Functional parameter indices for the relevant DetID
const std::vector<int> covarianceXsec::GetFuncParsIndexFromDetID(int DetID) {
  std::vector<int> returnVec;

  for (int i = 0; i < _fNumPar; ++i) {
	//std::cout << "TRYING TO SETUP FUNCTIONAL PARAMETER for " << i << " which is of type " << GetParamType(i) << std::endl;
    if ((GetParDetID(i) & DetID)) { //If parameter applies to required DetID
      if (strcmp(GetParamType(i), "Functional") == 0) { //If parameter is implemented as a functional param
		//std::cout << "Found Functional parameter" << std::endl;
		returnVec.push_back(i);
      }
    }
  }

  return returnVec;
}
// ********************************************

// ********************************************
// Scan the parameters, e.g. number of spline parameters, functional parameters, Far only parameters, and so on
//void covarianceXsec::ScanParameters() {
//  // ********************************************
//  
//  // Should be able to count the normalisation parameters from the covarianceXsec class
//  nNearSplineParams     = 0;
//  nSplineParamsUniq = 0; //!!decide what to do with this, doesn't currently have an Far analogue
//  nNearNormParams       = 0;
//  nNearFuncParams       = 0;
//  // Parameters that apply to Far only
//  
//  nFarSplineParams = 0;
//  nFarNormParams = 0;
//  nFaronlyNormParams = 0;
//  nFarFuncParams = 0;
//
//  int norm_counter = 0;
//
//  for (int i = 0; i < _fNumPar; ++i) {
//
//    //ETA - need to rethink this and just check against strings
//    /*bool isValidDetID = false;
//    int DetIDCounter = 0;
//    int ParamaDetID = GetXSecParamID(i,1);
//    //DB Loop over all supported DetIDs to ensure Root/XML inputs are familiar
//    for (int iKnownDetID=0;iKnownDetID<MaCh3Utils::nKnownDetIDs;iKnownDetID++) {
//      if ((ParamaDetID & MaCh3Utils::KnownDetIDsMap[iKnownDetID]) == MaCh3Utils::KnownDetIDsMap[iKnownDetID]) {
//	    isValidDetID = true;
//	    //DetIDCounter += MaCh3Utils::KnownDetIDsMap[iKnownDetID];
//      }
//    }
//    //DB Throw if Param DetID is unsupported. Also check that only supported DetIDs are contained in the param DetID
//    if (!isValidDetID || ((ParamaDetID - DetIDCounter)!=0)) {
//      std::cerr << "Unknown DetID:" << GetXSecParamID(i,1) << std::endl;
//      std::cerr << "DetIDCounter:" << DetIDCounter << std::endl;
//      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
//      throw;
//    }
//	*/
//    
//
//	///////////////
//	// ETA: for norm parameters the ND and FD treatment is the same so just need to check
//	// against the type of parameter. For splines and functional there are some
//	// different checks and hacks that need to be done for ND or FD.	
//	///////////////
//	
//	//ETA - adding in a counter for the number of norms as xsec_norm_kinematic_type is only of length of the number of norm parameters
//	//Not sure how this wasn't broken already??
//    //I think ScanParameters should just get the counters and set names? Even the setting names I don't like,
//    // as "Scan" implies it is just counting and printing stuff	
//	if(strcmp(GetParamType(i), "Norm") == 0){
//	  
//	  
//	  //XsecNorms4 tmp_xsec;
//	  //tmp_xsec.name=GetParName(i);
//
//	  //tmp_xsec.modes=_fNormModes[norm_counter];
//      //
//      
//
//	  //Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
//	  //We set a bool to see if any bounds exist so we can shortcircuit checking all of them every step
//	  bool haskinbounds=false;
//
//	  ////////////////////
//	  //New generic cuts things 
//	  ////////////////////
//	  //Only consider the kinematic string and the boundaries if you've actually given it a string to use...
//	  /*if( _fKinematicPars[i].size() > 0 ){
//		haskinbounds = true;
//
//		//ETA - This can be a vector :) can provide different kinematic variables to cut on
//		tmp_xsec.KinematicVarStr = _fKinematicPars[i];
//
//		//ETA - This can be a vector :) can provide different kinematic variables to cut on
//		std::vector< std::vector<double> > Selections(_fKinematicPars[i].size());
//
//		//ETA - push back kinematic type with dummy -999 since this needs to be converted into an enum for a kinematic type within
//		//a samplePDFFD daughter class
//		for(unsigned int KinVar_i = 0 ; KinVar_i < _fKinematicPars[i].size() ; ++KinVar_i) {
//		  Selections[KinVar_i].push_back(-999.9);
//		  Selections[KinVar_i].push_back(_fKinematicBounds[i][KinVar_i][0]);
//		  Selections[KinVar_i].push_back(_fKinematicBounds[i][KinVar_i][1]);
//		  //std::cout << "  - " << _fKinematicPars[i][KinVar_i] << " from " << Selections[KinVar_i][1] << " to " << Selections[KinVar_i][2] << std::endl;
//		}
//		tmp_xsec.Selection = Selections;
//	  }
//
//	  tmp_xsec.hasKinBounds=haskinbounds;
//	  //End of kinematic bound checking	
//	  // Set the global parameter index of the normalisation parameter
//	  tmp_xsec.index=i;
//	  */
//	  if((GetParDetID(i) & 24)){
//		//Add this parameter to the vector of parameters
//		//FarNormParams.push_back(tmp_xsec);
//		nFarNormParams++;
//
//		//non-Near affecting Far affecting norm parameters
//		if(!((GetParDetID(i) & 1))) {
//		  nFaronlyNormParams++;
//		}
//	  }
//	  
//	  if((GetParDetID(i) & 1)){
//		//NearNormParams.push_back(tmp_xsec);
//		nNearNormParams++;
//	  }
//	  else{std::cout << "Found DetId of " << GetParDetID(i) << std::endl;}
//	  norm_counter++;
//
//	}//End FarNormPars
//
//	/////////
//	//ETA:
//	//For splines and functional parameter there are some differences with the ND 
//	//and the FD still... so need to explicitly check against a hard-coded DetID
//	////////////
//
//    // Also make some helper arrays with all the Far parameters, so we have Near parameters, Far parameters, and all parameters
//    // These are parameters that have 24 or 25 in ID(i,1), so can be fancy and do a bitwise comparison to 01000
//	if ((GetParDetID(i) & 24) == 24) {//Far pars
//
//	  // Now check if it's a spline parameter or not
//	  //This needs to be updated to check against a string
//	  if (strcmp(GetParamType(i), "Spline") == 0) {//FarSplinePars
//		FarSplineParsNames.push_back(GetParName(i));
//		FarSplineParsIndex.push_back(i);
//
//		FarSplineModes.push_back(_fFDSplineModes[nFarSplineParams]);
//		nFarSplineParams++;
//
//		// Or a normalisation parameter
//	  } //End FarSplinePars
//	  else if (strcmp(GetParamType(i), "Functional") == 0){//Far functional parameter
//		nFarFuncParams++;
//		FarFuncParsNames.push_back(GetParName(i));
//		FarFuncParsIndex.push_back(i);
//	  }//End Far funcpars
//	  else if(!strcmp(GetParamType(i), "Norm") == 0){
//		std::cerr << "Found a parameter in covarianceXsec which wasn't Functional, Spline or Norm!" << std::endl;
//		std::cerr << "This is undefined behaviour currently, and implementation should change" << std::endl;
//		std::cerr << "Param " << GetParName(i) << " (param " << i << ") = " << GetParamType(i) << std::endl;
//		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
//		throw;
//	  }
//	}//End Far pars
//
//    //Near affecting parameters
//	/*if((GetParDetID(i) & 1) == 1) {//Near affecting parameters
//
//	  if (strcmp(GetParamType(i), "Functional") == 0) {//Near affecting func pars
//		NearfuncParsNames.push_back(GetParName(i));
//		NearfuncParsIndex.push_back(i);
//		nNearFuncParams++;
//	  }//End Near affecting func pars
//	  else if (strcmp(GetParamType(i), "Spline") == 0) {
//		//std::cout << "FOUND AN ND SPLINE PARAMETER!!" << std::endl;
//		//NearsplineParsNames.push_back(GetParName(i));
//		//NearsplineParsIndex.push_back(i);
//		//nNearSplineParams++;
//	  }//End Near affecting spline pars
//	  else if(!strcmp(GetParamType(i), "Norm") == 0){
//		std::cerr << "Found a parameter in covarianceXsec which wasn't Functional, Spline or Norm!" << std::endl;
//		std::cerr << "This is undefined behaviour currently, and implementation should change" << std::endl;
//		std::cerr << "Param " << GetParName(i) << " (param " << i << ") = " << GetParamType(i) << std::endl;
//		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
//		throw;
//	  }
//	  */
//	}//End Near affecting parameters
//  }
//  
//  // Now count the repeated parameters and save their indices and names
//  // This is useful because we only need to store 1 spline for e.g. neutrino and anti-neutrino interactions, but we still want to have separate neutrino/anti-neutrino parameters in covarianceXsec class
//  // This then feeds into samplePDFND2014
//  std::string paramName = "empty";
//  std::string nextName = "empty";
//  // Counter for unique splines (not repeated)
//  nSplineParamsUniq = 0;
//  // Counter for shared parameter (identified here by presence of "_nubar" in end of name)
//  // Maybe this should be worked on for the future, perhaps by having the input matrix specify instead which parameters share splines with each other; would require whoever makes the splines to talk to whoever makes the covariance matrix
//  nSplineParamsShare = 0;
//  int itCnt = 0;
//  for (std::vector<std::string>::iterator it = NearsplineParsNames.begin(); it != NearsplineParsNames.end(); it++) {
//
//    // The name of the current parameter
//    std::string nextName = *it;
//    std::string nextName_cut = nextName.substr(0, nextName.find_last_of("_"));
//
//
//    paramName = nextName_cut;
//    splineParsUniqIndex.push_back(NearsplineParsIndex[itCnt]);
//    splineParsUniqNames.push_back(*it);
//    splineParsShareToUniq.push_back(nSplineParamsUniq);
//
//    nSplineParamsUniq++;
//    itCnt++;
//  } // end for loop
//  return;
//} // end ScanParameters

// ********************************************
void covarianceXsec::initParams(double fScale) {
  // ********************************************

  for (int i = 0; i < _fNumPar; ++i) {
    //ETA - set the name to be xsec_% as this is what ProcessorMCMC expects
    _fNames[i] = "xsec_"+std::to_string(i);

    // Set covarianceBase parameters (Curr = current, Prop = proposed, Sigma = step)
    _fCurrVal[i] = _fPreFitValue[i];
    _fPropVal[i] = _fCurrVal[i];

    // Any param with nom == 1 should be > 0
    if (_fGenerated[i] == 1) {
      // If the _fPreFitValue is negative we should try to throw it above this
      // We don't really do this ever...
      while (_fPreFitValue[i] <= 0) {
        _fPreFitValue[i] = random_number[0]->Gaus(_fPreFitValue[i], fScale*TMath::Sqrt( (*covMatrix)(i,i) ));
      }
    }
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
  if (pca) TransferToPCA();
  CorrelateSteps();
}

// ********************************************
// Print everything we know about the inputs we're Getting
void covarianceXsec::Print() {
// ********************************************

  // Get the precision so we can go back to normal for cout later
  //std::streamsize ss = std::cout.precision();

  MACH3LOG_INFO("#################################################");
  MACH3LOG_INFO("Printing covarianceXsec:");

  std::cout<<"======================================================================================================================"<<std::endl;
  std::cout << std::left << std::setw(5) << "#" << std::setw(2) << "|" << std::setw(25) << "Name" << std::setw(2) << "|" << std::setw(10) << "Nom." << std::setw(2) << "|" << std::setw(10) << "Prior" << std::setw(2) << "|" << std::setw(15) << "Error" << std::setw(2) << "|" << std::setw(10) << "Lower" << std::setw(2) << "|" << std::setw(10) << "Upper" << "|" << std::setw(10) << "StepScale" << "|" << std::setw(5) << "DetID" << std::endl;;
  std::cout<<"----------------------------------------------------------------------------------------------------------------------"<<std::endl;

  for (int i = 0; i < GetNumParams(); i++) {
    std::cout << std::left << std::setprecision(3) << std::setw(5) << i << std::setw(2) << "|" << std::setw(25) << GetParFancyName(i) << std::setw(2) << "|" << std::setw(10) << _fGenerated[i] << std::setw(2) << "|" << std::setw(10) << _fPreFitValue[i] << std::setw(2) << "|" << "+/- " << std::setw(11) << _fError[i] << std::setw(2) << "|" << std::setw(10) << _fLowBound[i] << std::setw(2) << "|" << std::setw(10) << _fUpBound[i] << "|" << std::setw(10) << _fIndivStepScale[i] << "|" << _fDetID[i] << std::endl;
  }
  std::cout<<"======================================================================================================================"<<std::endl;

  // Output the normalisation parameters as a sanity check!
  MACH3LOG_INFO("Normalisation parameters:  {}", NormParams.size());

  //KS: Consider making some class producing table..
  MACH3LOG_INFO("┌────┬──────────┬────────────────────┬──────────┬──────────┬──────────┐");
  MACH3LOG_INFO("│#   │Global #  │Name                │Int. mode │Target    │pdg       │");
  MACH3LOG_INFO("├────┼──────────┼────────────────────┼──────────┼──────────┼──────────┤");

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

    MACH3LOG_INFO("│{: <4}│{: <10}│{: <20}│{: <10}│{: <10}│{: <10}│", i, NormParams[i].index, NormParams[i].name, intModeString, targetString, pdgString);
  }
  MACH3LOG_INFO("└────┴──────────┴────────────────────┴──────────┴──────────┴──────────┘");

  std::vector<int> SplineParsIndex;
  for (int i = 0; i < _fNumPar; ++i)
  {
    if (strcmp(GetParamType(i), "Spline") == 0) { SplineParsIndex.push_back(i); }
  }

  MACH3LOG_INFO("Spline parameters: {}", SplineParsIndex.size());


  MACH3LOG_INFO("=====================================================");
  MACH3LOG_INFO("{:<4} {:<2} {:<10} {:<2} {:<30} {:<2}", "#", "|", "Name", "|", "Spline Interpolation", "|");
  MACH3LOG_INFO("-----------------------------------------------------");
  for (auto index = SplineParsIndex.begin() ; index != SplineParsIndex.end() ; ++index){
	MACH3LOG_INFO("{:<4} {:<2} {:<10} {:<2} {:<30} {:<2}", *index, "|", GetParFancyName(*index), "|", SplineInterpolation_ToString(GetParSplineInterpolation(*index)), "|");
  }
  MACH3LOG_INFO("=====================================================");

  std::vector<int> FuncParsIndex;
  for (int i = 0; i < _fNumPar; ++i)
  {
    if (strcmp(GetParamType(i), "Functional") == 0) { FuncParsIndex.push_back(i); }
  }

  MACH3LOG_INFO("Functional parameters: {}", FuncParsIndex.size());
  MACH3LOG_INFO("=================================");
  MACH3LOG_INFO("{0:4} {1:2} {2:10}", "#", "|", "Name");
  MACH3LOG_INFO("---------------------------------");
  for (unsigned int i = 0; i < FuncParsIndex.size(); ++i) {
    MACH3LOG_INFO("{0:4} {1:2} {2:10}", std::to_string(i), "|", GetParFancyName(FuncParsIndex[i]));
  }
  MACH3LOG_INFO("=================================");

} // End

// ********************************************
// Sets the proposed Flux parameters to the prior values
void covarianceXsec::setFluxOnlyParameters() {
// ********************************************
  if(!pca)
  {
    for (int i = 0; i < _fNumPar; i++)
    {
      if(isFlux[i]) _fPropVal[i] = _fPreFitValue[i];
    }
  }
  else
  {
    MACH3LOG_ERROR("setFluxOnlyParameters not implemented for PCA");
    throw;
  }

}

// ********************************************
// Sets the proposed Flux parameters to the prior values
void covarianceXsec::setXsecOnlyParameters() {
// ********************************************
  if(!pca)
  {
    for (int i = 0; i < _fNumPar; i++)
    {
      if(!isFlux[i]) _fPropVal[i] = _fPreFitValue[i];
    }
  }
  else
  {
    MACH3LOG_ERROR("setXsecOnlyParameters not implemented for PCA");
    throw;
  }
}
