#include "covarianceXsec.h"

// ********************************************
covarianceXsec::covarianceXsec(const char *name, const char *file,
   	                           double threshold,int FirstPCAdpar,
							   int LastPCAdpar)
  : covarianceBase(name, file,0,threshold,FirstPCAdpar,LastPCAdpar) {
// ********************************************
 
  int ncols = _fNumPar;

  // Leave these as arrays for backwards compatibility
  // Probably should be made into vectors but at least size isn't hard-coded ^,^
  // Load up the parameters id (nPars x nCols big, contains information about if these parameters are spline, normalisation or some other parameter)
  //xsec_param_id_a   = new int*[_fNumPar];
  //for (int i = 0; i < _fNumPar; ++i) {
  //  xsec_param_id_a[i] = new int[ncols];
  //}
  xsec_param_nom_a  = new double[_fNumPar];
  xsec_param_lb_a   = new double[_fNumPar];
  xsec_param_ub_a   = new double[_fNumPar];
  xsec_param_prior_a = new double[_fNumPar];

  //DB Resize step scale vector to expected size
  xsec_stepscale_vec.resize(_fNumPar);

  //ETA - We really don't need this loop and can get rid of it shortly
  for (int i = 0; i < _fNumPar; i++) {
    // Fill the prior central value
	xsec_param_nom_a[i] = _fPreFitValue[i];

    // Fill the lower bound
    xsec_param_lb_a[i]  = _fLowBound[i];

    // Fill the upper bound
    xsec_param_ub_a[i]  = _fUpBound[i];
    // Fill the prior uncertainty
	xsec_param_prior_a[i]= _fError[i];
    // DB Fill the stepscales vector
    xsec_stepscale_vec[i] = _fIndivStepScale[i];

	//for (int j = 0; j < ncols; j++) {
    //  xsec_param_id_a[i][j] = (*xsec_param_id)(i,j);
    //}
  } // end the for loop

  //infile->Close();
  //delete infile;
  // Scan through the input parameters and find which are normalisation, which are splines, and so on
  ScanParameters();

  std::cout << "Constructing instance of covarianceXsec" << std::endl;
  initParams(0.001);
  // Print
  Print();

  // Set the parameters
  // setParameters; why is this done in testFGD2 and not here?
}

// ********************************************
// ETA - YAML constructor
// this will replace the root file constructor but let's keep it in
// to do some validations
covarianceXsec::covarianceXsec(const char *YAMLFile)
  : covarianceBase(YAMLFile) {

  std::cout << "About to ParseYaml" << std::endl;
  ParseYAML(YAMLFile);

  int ncols = _fNumPar;
  // Leave these as arrays for backwards compatibility
  // Probably should be made into vectors but at least size isn't hard-coded ^,^
  // Load up the parameters id (nPars x nCols big, contains information about if these parameters are spline, normalisation or some other parameter)
  //xsec_param_id_a   = new int*[_fNumPar];
  //for (int i = 0; i < _fNumPar; ++i) {
  //  xsec_param_id_a[i] = new int[ncols];
  //}
  xsec_param_nom_a  = new double[_fNumPar];
  xsec_param_lb_a   = new double[_fNumPar];
  xsec_param_ub_a   = new double[_fNumPar];
  xsec_param_prior_a = new double[_fNumPar];

  //DB Resize step scale vector to expected size
  xsec_stepscale_vec.resize(_fNumPar);

  //ETA - again this really doesn't need to be hear...
  for (int i = 0; i < _fNumPar; i++) {
    // Fill the nominal
	xsec_param_nom_a[i] = _fPreFitValue[i];

    // Fill the lower bound
    xsec_param_lb_a[i]  = _fLowBound[i];

    // Fill the upper bound
    xsec_param_ub_a[i]  = _fUpBound[i];

    // Fill the prior
	xsec_param_prior_a[i]= _fError[i];

	// Sort out the print length
    if(_fNames[i].length() > PrintLength) PrintLength = _fNames[i].length();

    // DB Fill the stepscales vector
    xsec_stepscale_vec[i] = _fIndivStepScale[i];
  } // end the for loop

  //infile->Close();
  //delete infile;
  // Scan through the input parameters and find which are normalisation, which are splines, and so on
  ScanParameters();

  std::cout << "Constructing instance of covarianceXsec" << std::endl;
  //initParams(0.001);
  // Print
  Print();

  // Set the parameters
  // setParameters; why is this done in testFGD2 and not here?
}


void covarianceXsec::ParseYAML(const char* FileName)
{
  std::cout << "Let's read some YAML!" << std::endl;

  _fYAMLDoc = YAML::LoadFile(FileName);
  _fNumPar = _fYAMLDoc["Systematics"].size();
  _fNames = std::vector<std::string>(_fNumPar);
  _fGenerated = std::vector<double>(_fNumPar);
  _fPreFitValue = std::vector<double>(_fNumPar);
  _fError = std::vector<double>(_fNumPar);
  _fLowBound = std::vector<double>(_fNumPar);
  _fUpBound = std::vector<double>(_fNumPar);
  _fIndivStepScale = std::vector<double>(_fNumPar);
  _fFlatPrior = std::vector<bool>(_fNumPar);
  _fDetID = std::vector<int>(_fNumPar);
  _fCovMatrix = new TMatrixDSym(_fNumPar);
  _fParamType = std::vector<std::string>(_fNumPar);
  _fDetString = std::vector<std::string>(_fNumPar);

  std::cout << "Found " << _fNumPar << " systematics in yaml" << std::endl;

  //Vector of vectors of strings to contain potentially multiple variables that
  //might be cut on
  _fKinematicPars = std::vector<std::vector<std::string>>(_fNumPar);
  //Vector of vector of ints to contain the lower and upper bounds of a cut
  //for a particular kinematic variables
  _fKinematicBounds = std::vector<std::vector<std::vector<double>>>(_fNumPar);

  int i=0;

  std::vector<std::map<std::string,double>> Correlations(_fNumPar);
  std::map<std::string, int> CorrNamesMap;

  //ETA - read in the systematics. Would be good to add in some checks to make sure
  //that there are the correct number of entries i.e. are the _fNumPars for Names,
  //PreFitValues etc etc.
  for (auto const &param : _fYAMLDoc["Systematics"]) {
     //std::cout << param["Systematic"]["Names"]["ParameterName"].as<std::string>() << std::endl;

     _fNames[i] = (param["Systematic"]["Names"]["ParameterName"].as<std::string>());
     _fPreFitValue[i] = (param["Systematic"]["ParameterValues"]["PreFitValue"].as<double>());
     _fGenerated[i] = (param["Systematic"]["ParameterValues"]["Generated"].as<double>());
     _fIndivStepScale[i] = (param["Systematic"]["StepScale"]["MCMC"].as<double>());
     _fDetID[i] = (param["Systematic"]["DetID"].as<int>());
     _fError[i] = (param["Systematic"]["Error"].as<double>());
	 _fParamType[i] = (param["Systematic"]["Type"].as<std::string>());

	 //ETA - a bit of a fudge but works
	 std::vector<double> TempBoundsVec = param["Systematic"]["ParameterBounds"].as<std::vector<double>>();
     _fLowBound[i] = TempBoundsVec[0];
     _fUpBound[i] = TempBoundsVec[1];

	 //ETA - now for parameters which are optional and have default values
	 if (param["Systematic"]["FlatPrior"]) {
	   _fFlatPrior[i] = param["Systematic"]["FlatPrior"].as<bool>();
	 } else {
	   _fFlatPrior[i] = false;
	 }

	 //Fill the map to get the correlations later as well
     CorrNamesMap[param["Systematic"]["Names"]["ParameterName"].as<std::string>()]=i;
	 std::string ParamType = param["Systematic"]["Type"].as<std::string>();
	 int nFDSplines;
	 //Now load in varaibles for spline systematics only
	 if (ParamType.find("Spline") != std::string::npos) {
	   //std::cout << "Reading in a Spline Parameter" << std::endl;

	   if (param["Systematic"]["SplineInformation"]["FDSplineName"]) {
		 _fFDSplineNames.push_back(param["Systematic"]["SplineInformation"]["FDSplineName"].as<std::string>());
		 nFDSplines++;
	   }

	   if (param["Systematic"]["SplineInformation"]["FDMode"]) {
		 //std::cout << "Pushing back _fFDSplineModes for param " << i << std::endl;
		 _fFDSplineModes.push_back(param["Systematic"]["SplineInformation"]["FDMode"].as<std::vector<int>>());
		 std::cout << "_fFDSplineModes is of size " << _fFDSplineModes.size() << std::endl;
		 std::cout << "_fFDSplineModes[0] is of size " << _fFDSplineModes[0].size() << std::endl;
	   }

	   if (param["Systematic"]["SplineInformation"]["NDSplineName"]) {
		 _fNDSplineNames.push_back(param["Systematic"]["SplineInformation"]["NDSplineName"].as<std::string>());
		 //MASSIVE HACKK!!! This only works because we only need the interpolation type at the ND
		 //Now get the Spline interpolation type
		 if (param["Systematic"]["SplineInformation"]["InterpolationType"]){
		   for(int InterpType = 0; InterpType < kSplineInterpolations ; InterpType++){
			 if(param["Systematic"]["SplineInformation"]["InterpolationType"].as<std::string>() == SplineInterpolation_ToString(SplineInterpolation(InterpType)))
			 {
			   _fSplineInterpolationType.push_back(SplineInterpolation(InterpType));
			 }
		   }
		 }
	   }

	   else{std::cout << "PROBLEM!!! No SPLINEINFORMATION FOR " << param["Systematic"]["Names"]["ParameterName"].as<std::string>().c_str() << std::endl;}



	 } else if(param["Systematic"]["Type"].as<std::string>() == "Norm") {
	   //std::cout << "Norm parameter" << std::endl;
 
	   //Empty DummyVector can be used to specify no cut for mode, target and neutrino flavour
	   std::vector<int> DummyModeVec;

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
	   std::cout << "Found a functional parameter!!" << std::endl;
	    
	 }
     else{
       std::cerr << "[ERROR] Given unrecognised systematic type: " << param["Systematic"]["Type"].as<std::string>().c_str() << std::endl; 
       std::cerr << "[ERROR] Expecting \"Norm\", \"Spline\" or \"Functional\"" << std::endl;
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
       
	 //Also loop through the correlations
	 if(param["Systematic"]["Correlations"]) {
	   //std::cout << "FOUND " << param["Systematic"]["Correlations"].size() << " CORRELATIONS!" << std::endl;
	   for(unsigned int Corr_i = 0 ; Corr_i < param["Systematic"]["Correlations"].size()  ; ++Corr_i){
		 for (YAML::const_iterator it=param["Systematic"]["Correlations"][Corr_i].begin();it!=param["Systematic"]["Correlations"][Corr_i].end();++it) {
		   //std::cout << "Correlation with " << it->first.as<std::string>() << " of " << it->second.as<double>() << std::endl;
           Correlations[i][it->first.as<std::string>()] = it->second.as<double>();
		 }
	   }
	 }
	 i++;
  }
  
	 
  //ETA
  //Now that we've been through all systematic let's fill the covmatrix
  //This makes the root TCov from YAML
  for(int i=0; i < _fNumPar; i++) {
	(*_fCovMatrix)(i,i)=_fError[i]*_fError[i];
	//Get the map of parameter name to correlation fomr the Correlations object
	for (auto const& [key, val] : Correlations[i]) {
	  int index = -1;

	  //If you found the parameter name then get the index
	  if (CorrNamesMap.find(key) != CorrNamesMap.end()) {
		index=CorrNamesMap[key];
	  }
	  else {
		std::cout << "Parameter " << key << " not in list! Check your spelling?" << std::endl;
		exit(5);
	  }

	  //
	  double Corr1 = val;
	  double Corr2 = 0;
	  if(Correlations[index].find(_fNames[i]) != Correlations[index].end()) {
		Corr2 = Correlations[index][_fNames[i]];
		//Do they agree to better than float precision?
		if(std::abs(Corr2 - Corr1) > FLT_EPSILON) {
		  std::cout << "Correlations are not equal between " << _fNames[i] << " and " << key << std::endl;
		  std::cout << "Got : " << Corr2  << " and " << Corr1 << std::endl;
		  exit(5);
		}
	  } else {
		std::cout << "Correlation does not appear reciprocally between " << _fNames[i] << " and " << key << std::endl;
		exit(5);
	  }
	  (*_fCovMatrix)(i,index)= (*_fCovMatrix)(index,i) = Corr1*_fError[i]*_fError[index];
	}
  } 

  //Now make positive definite
  MakePosDef(_fCovMatrix);

  setCovMatrix(_fCovMatrix);
   
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
  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
		returnVal += 1;
      }
    }
  }

  return returnVal;
}
// ********************************************

// ********************************************
// DB Grab the Spline Names for the relevant DetID
const std::vector<std::string> covarianceXsec::GetSplineParsNamesFromDetID(int DetID) {

  std::vector<std::string> returnVec;
  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
        returnVec.push_back(GetParName(i));
      }
    }
  }

  return returnVec;
}
// ********************************************
//

// ETA - this is a complete fudge for now and is only here because on
// T2K the splines at ND280 and SK have different names... 
// this is completely temporary and will be removed in the future
const std::vector<std::string> covarianceXsec::GetFDSplineFileParsNamesFromDetID(int DetID) {

  std::vector<std::string> returnVec;
  int FDSplineCounter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
	//std::cout << " Param i has DetID " << GetXsecParamDetID(i) << " from yaml" << std::endl;
	if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
	  if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
		//This is horrible but has to be here whilst the ND and FD splines are named
		//differently on T2K. This is easy to fix but isn't currently available.
		//std::cout << "FDSplineIndex is " << std::endl;
		//std::cout << "Getting FD spline name " << _fFDSplineNames[FDSplineCounter] << " compared DetID " << DetID << " with " << GetXsecParamDetID(i) << std::endl;
		returnVec.push_back(_fFDSplineNames[FDSplineCounter]); //Append spline name
		FDSplineCounter++;
	  }
	}
  }
  return returnVec;
}

// ETA - this is another fudge for now and is only here because on
// T2K the splines at ND280 and SK have different names... 
// this is completely temporary and will be removed in the future
const std::vector<std::string> covarianceXsec::GetNDSplineFileParsNamesFromDetID(int DetID) {

  std::vector<std::string> returnVec;
  int NDSplineCounter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
	if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
	  if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
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
const std::vector<std::string> covarianceXsec::GetSplineFileParsNamesFromDetID(int DetID) {
  std::vector<std::string> returnVec;

  int FDSplineCounter = 0;
  int NDSplineCounter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
		//This is horrible but has to be here whilst the ND and FD splines are named
		//differently on T2K. This is easy to fix but isn't currently available.
		if(DetID == 1){
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

// ********************************************
// DB Grab the Spline Modes for the relevant DetID
const std::vector< std::vector<int> > covarianceXsec::GetSplineModeVecFromDetID(int DetID) {
  std::vector< std::vector<int> > returnVec;

  //Need a counter or something to correctly get the index in _fFDSplineModes since it's not of length nPars
  //Should probably just make a std::map<std::string, int> for param name to FD spline index
  int nFDSplineCounter = 0;
  for (int i = 0; i < _fNumPar; ++i) {
	if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
	  if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
		returnVec.push_back(_fFDSplineModes[nFDSplineCounter]);	
		nFDSplineCounter++;
	  }
	}
  }

  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the Spline Indices for the relevant DetID
const std::vector<int> covarianceXsec::GetSplineParsIndexFromDetID(int DetID) {
  std::vector<int> returnVec;

  for (int i = 0; i < _fNumPar; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
		returnVec.push_back(i);
      }
    }
  }

  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the Normalisation parameters for the relevant DetID
// ETA - I think this doesn't need to be the same as scanParameters, haven't we already got this info??
const std::vector<XsecNorms4> covarianceXsec::GetNormParsFromDetID(int DetID) {
  std::vector<XsecNorms4> returnVec;
  int norm_counter = 0;

  for (int i = 0; i < _fNumPar; ++i) {
	if (strcmp(GetXsecParamType(i), "Norm") == 0) { //If parameter is implemented as a normalisation

	  if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID   

		std::vector<int> temp;

		XsecNorms4 norm;
		norm.name=GetParName(i);

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
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Functional") == 0) { //If parameter is implemented as a functional parameter
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
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Functional") == 0) { //If parameter is implemented as a functional param
		returnVec.push_back(GetParName(i));
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
	//std::cout << "TRYING TO SETUP FUNCTIONAL PARAMETER for " << i << " which is of type " << GetXsecParamType(i) << std::endl;
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Functional") == 0) { //If parameter is implemented as a functional param
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
void covarianceXsec::ScanParameters() {
  // ********************************************
  
  // Should be able to count the normalisation parameters from the covarianceXsec class
  nNearSplineParams     = 0;
  nSplineParamsUniq = 0; //!!decide what to do with this, doesn't currently have an Far analogue
  nNearNormParams       = 0;
  nNearFuncParams       = 0;
  // Parameters that apply to Far only
  
  nFarSplineParams = 0;
  nFarNormParams = 0;
  nFaronlyNormParams = 0;
  nFarFuncParams = 0;

  int norm_counter = -1;

  for (int i = 0; i < _fNumPar; ++i) {

    //ETA - need to rethink this and just check against strings
    /*bool isValidDetID = false;
    int DetIDCounter = 0;
    int ParamaDetID = GetXSecParamID(i,1);
    //DB Loop over all supported DetIDs to ensure Root/XML inputs are familiar
    for (int iKnownDetID=0;iKnownDetID<MaCh3Utils::nKnownDetIDs;iKnownDetID++) {
      if ((ParamaDetID & MaCh3Utils::KnownDetIDsMap[iKnownDetID]) == MaCh3Utils::KnownDetIDsMap[iKnownDetID]) {
	    isValidDetID = true;
	    //DetIDCounter += MaCh3Utils::KnownDetIDsMap[iKnownDetID];
      }
    }
    //DB Throw if Param DetID is unsupported. Also check that only supported DetIDs are contained in the param DetID
    if (!isValidDetID || ((ParamaDetID - DetIDCounter)!=0)) {
      std::cerr << "Unknown DetID:" << GetXSecParamID(i,1) << std::endl;
      std::cerr << "DetIDCounter:" << DetIDCounter << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
	*/
    
    if(GetParName(i).find("b_")==0)isFlux.push_back(true);
    else isFlux.push_back(false);

	///////////////
	// ETA: for norm parameters the ND and FD treatment is the same so just need to check
	// against the type of parameter. For splines and functional there are some
	// different checks and hacks that need to be done for ND or FD.	
	///////////////
	
	//ETA - adding in a counter for the number of norms as xsec_norm_kinematic_type is only of length of the number of norm parameters
	//Not sure how this wasn't broken already??
	//This needs to be updated to check against a string
	if(strcmp(GetXsecParamType(i), "Norm") == 0){
	  
	  XsecNorms4 tmp_xsec;
	  tmp_xsec.name=GetParName(i);

	  tmp_xsec.modes=_fNormModes[nFarNormParams];

	  //Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
	  //We set a bool to see if any bounds exist so we can shortcircuit checking all of them every step
	  bool haskinbounds=false;

	  ////////////////////
	  //New generic cuts things 
	  ////////////////////
	  //Only consider the kinematic string and the boundaries if you've actually given it a string to use...
	  if( _fKinematicPars[i].size() > 0 ){
		haskinbounds = true;

		//ETA - This can be a vector :) can provide different kinematic variables to cut on
		tmp_xsec.KinematicVarStr = _fKinematicPars[i];

		//ETA - This can be a vector :) can provide different kinematic variables to cut on
		std::vector< std::vector<double> > Selections(_fKinematicPars[i].size());

		//ETA - push back kinematic type with dummy -999 since this needs to be converted into an enum for a kinematic type within
		//a samplePDFFD daughter class
		for(unsigned int KinVar_i = 0 ; KinVar_i < _fKinematicPars[i].size() ; ++KinVar_i) {
		  Selections[KinVar_i].push_back(-999.9);
		  Selections[KinVar_i].push_back(_fKinematicBounds[i][KinVar_i][0]);
		  Selections[KinVar_i].push_back(_fKinematicBounds[i][KinVar_i][1]);
		  //std::cout << "  - " << _fKinematicPars[i][KinVar_i] << " from " << Selections[KinVar_i][1] << " to " << Selections[KinVar_i][2] << std::endl;
		}
		tmp_xsec.Selection = Selections;
	  }

	  tmp_xsec.hasKinBounds=haskinbounds;
	  //End of kinematic bound checking	

	  // Set the global parameter index of the normalisation parameter
	  tmp_xsec.index=i;
	  if(GetXsecParamDetID(i) & 24){
		//Add this parameter to the vector of parameters
		FarNormParams.push_back(tmp_xsec);
		nFarNormParams++;

		//non-Near affecting Far affecting norm parameters
		if (!(GetXsecParamDetID(i) & 1)) {
		  nFaronlyNormParams++;
		}
	  }
	  else if(GetXsecParamDetID(i) & 1){
		NearNormParams.push_back(tmp_xsec);
		nNearNormParams++;
	  }
	  norm_counter++;

	}//End FarNormPars

	/////////
	//ETA:
	//For splines and functional parameter there are some differences with the ND 
	//and the FD still... so need to explicitly check against a hard-coded DetID
	////////////

    // Also make some helper arrays with all the Far parameters, so we have Near parameters, Far parameters, and all parameters
    // These are parameters that have 24 or 25 in ID(i,1), so can be fancy and do a bitwise comparison to 01000
	if ((GetXsecParamDetID(i) & 24) == 24) {//Far pars

	  // Now check if it's a spline parameter or not
	  //This needs to be updated to check against a string
	  if (strcmp(GetXsecParamType(i), "Spline") == 0) {//FarSplinePars
		FarSplineParsNames.push_back(GetParName(i));
		FarSplineParsIndex.push_back(i);

		FarSplineModes.push_back(_fFDSplineModes[nFarSplineParams]);
		nFarSplineParams++;

		// Or a normalisation parameter
	  } //End FarSplinePars
	  else if (strcmp(GetXsecParamType(i), "Functional") == 0){//Far functional parameter
		nFarFuncParams++;
		FarFuncParsNames.push_back(GetParName(i));
		FarFuncParsIndex.push_back(i);
	  }//End Far funcpars
	  else if(!strcmp(GetXsecParamType(i), "Norm") == 0){
		std::cerr << "Found a parameter in covarianceXsec which wasn't Functional, Spline or Norm!" << std::endl;
		std::cerr << "This is undefined behaviour currently, and implementation should change" << std::endl;
		std::cerr << "Param " << GetParName(i) << " (param " << i << ") = " << GetXsecParamType(i) << std::endl;
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		throw;
	  }
	}//End Far pars

    //Near affecting parameters
	if((GetXsecParamDetID(i) & 1) == 1) {//Near affecting parameters

	  if (strcmp(GetXsecParamType(i), "Functional") == 0) {//Near affecting func pars	
		NearfuncParsNames.push_back(GetParName(i));
		NearfuncParsIndex.push_back(i);
		nNearFuncParams++;
	  }//End Near affecting func pars
	  else if (strcmp(GetXsecParamType(i), "Spline") == 0) {
		//std::cout << "FOUND AN ND SPLINE PARAMETER!!" << std::endl;
		NearsplineParsNames.push_back(GetParName(i));
		NearsplineParsIndex.push_back(i);
		nNearSplineParams++;
	  }//End Near affecting spline pars
	  else if(!strcmp(GetXsecParamType(i), "Norm") == 0){
		std::cerr << "Found a parameter in covarianceXsec which wasn't Functional, Spline or Norm!" << std::endl;
		std::cerr << "This is undefined behaviour currently, and implementation should change" << std::endl;
		std::cerr << "Param " << GetParName(i) << " (param " << i << ") = " << GetXsecParamType(i) << std::endl;
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		throw;
	  }
	}//End Near affecting parameters
  }
  
  // Now count the repeated parameters and save their indices and names
  // This is useful because we only need to store 1 spline for e.g. neutrino and anti-neutrino interactions, but we still want to have separate neutrino/anti-neutrino parameters in covarianceXsec class
  // This then feeds into samplePDFND2014
  std::string paramName = "empty";
  std::string nextName = "empty";
  // Counter for unique splines (not repeated)
  nSplineParamsUniq = 0;
  // Counter for shared parameter (identified here by presence of "_nubar" in end of name)
  // Maybe this should be worked on for the future, perhaps by having the input matrix specify instead which parameters share splines with each other; would require whoever makes the splines to talk to whoever makes the covariance matrix
  nSplineParamsShare = 0;
  int itCnt = 0;
  for (std::vector<std::string>::iterator it = NearsplineParsNames.begin(); it != NearsplineParsNames.end(); it++) {

    // The name of the current parameter
    std::string nextName = *it;
    std::string nextName_cut = nextName.substr(0, nextName.find_last_of("_"));

    //DEPRECATED
    //KS: Shared splines are no longer used and looking for "_nubar" leads to confusion
    //and even to seg faults, right now none of splines will be assgined as "shared"
    //we can easily revert this change
    //
    // Anti-neutrino parameters will share the same spline name with the neutrino counter-parts
    //if (nextName.find("_nubar") != std::string::npos) {
    if (false) { 
      splineParsShareIndex.push_back(NearsplineParsIndex[itCnt]);
      splineParsShareNames.push_back(*it);
      // Which uniq spline does the parameter share with
      splineParsShareToUniq.push_back(splineParsShareToUniq.back());
      nSplineParamsShare++;

      itCnt++;
      continue;
    }

    paramName = nextName_cut;
    splineParsUniqIndex.push_back(NearsplineParsIndex[itCnt]);
    splineParsUniqNames.push_back(*it);
    splineParsShareToUniq.push_back(nSplineParamsUniq);

    nSplineParamsUniq++;
    itCnt++;
  } // end for loop
  return;
} // end ScanParameters

// ********************************************
void covarianceXsec::initParams(double fScale) {
  // ********************************************

  for (int i = 0; i < size; ++i) {
    char pname[10];
    sprintf(pname, "xsec_%i", i);
	_fNames[i] = std::string(pname);

    // This just assign the initial parameters to the prior
    _fPreFitValue[i] = _fPreFitValue[i];
    // Any param with nom == 1 should be > 0
    if (_fPreFitValue[i] == 1) {
      // If the _fPreFitValue is negative we should try to throw it above this
      // We don't really do this ever...
      while (_fPreFitValue[i] <= 0) {
        _fPreFitValue[i] = random_number[0]->Gaus(_fPreFitValue[i], 0.1*fScale*TMath::Sqrt( (*covMatrix)(i,i) ));
      }
    }

    // Set covarianceBase parameters (Curr = current, Prop = proposed, Sigma = step)
    _fCurrVal[i] = _fPreFitValue[i];
    _fPropVal[i] = _fCurrVal[i];
    _fError[i] = fScale;
  }

  //DB Set Individual Step scale for PCA parameters to the lastpcadpar fIndivStepScale because the step scale for those parameters is set by 'eigen_values[i]' but needs an overall step scale
  //   However, individual step scale for non-PCA parameters needs to be set correctly
  if (pca) {
    for (int i=FirstPCAdpar;i<=LastPCAdpar;i++) {
      xsec_stepscale_vec[i] = xsec_stepscale_vec[LastPCAdpar-1];
    }
  }

  setIndivStepScale(xsec_stepscale_vec);
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
  std::streamsize ss = std::cout.precision();

  std::cout << "#################################################" << std::endl;
  std::cout << "Printing covarianceXsec:" << std::endl;

    std::cout << "Number of parameters: " << GetNumParams() << std::endl;
  std::cout << std::endl;

  std::cout << "Global parameter map:" << std::endl;
  std::cout << std::left << std::setw(5) << "#" << std::setw(2) << "|" << std::setw(25) << "Name" << std::setw(2) << "|" << std::setw(10) << "Nom." << std::setw(2) << "|" << std::setw(10) << "Prior" << std::setw(2) << "|" << std::setw(15) << "Error" << std::setw(2) << "|" << std::setw(10) << "Lower" << std::setw(2) << "|" << std::setw(10) << "Upper" << "|" << std::setw(15) << "IndivStepScale" << std::endl;;

  for (int i = 0; i < GetNumParams(); i++) {
    std::cout << std::left << std::setprecision(3) << std::setw(5) << i << std::setw(2) << "|" << std::setw(25) << GetParName(i) << std::setw(2) << "|" << std::setw(10) << _fGenerated[i] << std::setw(2) << "|" << std::setw(10) << _fPreFitValue[i] << std::setw(2) << "|" << "+/- " << std::setw(11) << _fError[i] << std::setw(2) << "|" << std::setw(10) << _fLowBound[i] << std::setw(2) << "|" << std::setw(10) << _fUpBound[i] << "|" << std::setw(15) << _fIndivStepScale[i] << std::endl;
  }

  std::cout << std::endl;

  // Output the normalisation parameters as a sanity check!
  std::cout << "Near detector normalisation parameters:" << nNearNormParams << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(2) << "|" << std::setw(10) << "Global #" << std::setw(2) << "|" << std::setw(20) << "Name" << std::setw(2) << "|" << std::setw(10) << "Int. mode" << std::setw(2) << "|" << std::setw(10) << "Target" << std::setw(2) << "|" << std::setw(10) << "Type" << std::endl;
  for (int i = 0; i < nNearNormParams; ++i) {
    std::cout << std::setw(4) << i << std::setw(2) << "|" << std::setw(10) << NearNormParams.at(i).index << std::setw(2) << "|" << std::setw(20) << NearNormParams.at(i).name << std::setw(2) << "|" << std::setw(10);
    for(int j = 0; j < int((NearNormParams.at(i).modes).size()); j++){
      std::cout<< NearNormParams.at(i).modes.at(j) <<" ";
    }
    std::cout<< std::setw(2) << "|" << " ";
    for (int j = 0; j < int((NearNormParams.at(i).targets).size()); j++) {
      std::cout << (NearNormParams.at(i).targets).at(j) << " ";
    }

    std::cout << std::setw(2) << "|" << " ";

    for (int j = 0; j < int((NearNormParams.at(i).pdgs).size()); j++) {
      std::cout << (NearNormParams.at(i).pdgs).at(j) << " ";
    }

    std::cout << std::endl;
  }

  std::cout << std::endl;

  // Start Far printing

  std::cout << std::endl;

  std::cout << "Far detector spline parameters: " << nFarSplineParams << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(2) << "|" << std::setw(10) << "Global #" << std::setw(2) << "|" << std::setw(20) << "Name" << std::endl;
  for (int i = 0; i < nFarSplineParams; ++i) {
    std::cout << std::setw(4) << i << std::setw(2) << "|" << std::setw(10) << FarSplineParsIndex.at(i) << std::setw(2) << "|" << std::setw(20) << FarSplineParsNames.at(i) << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Far detector functional parameters: " << nFarFuncParams << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(2) << "|" << std::setw(10) << "Global #" << std::setw(2) << "|" << std::setw(20) << "Name" << std::endl;
  for (int i = 0; i < nFarFuncParams; ++i) {
    std::cout << std::setw(4) << i << std::setw(2) << "|" << std::setw(10) << FarFuncParsIndex.at(i) << std::setw(2) << "|" << std::setw(20) << FarFuncParsNames.at(i) << std::endl;
  }
  std::cout << std::endl;

  // Output the Far normalisation parameters as a sanity check!
  std::cout << "Far detector normalisation parameters: " << nFarNormParams << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(2) << "|" << std::setw(10) << "Global #" << std::setw(2) << "|" << std::setw(20) << "Name" << std::setw(2) << "|" << std::setw(10) << "Int. mode" << std::setw(2) << "|" << std::setw(10) << "Target" << std::setw(2) << "|" << std::setw(10) << "Type" << std::endl;
  for (int i = 0; i < nFarNormParams; ++i) {
    std::cout << std::setw(4) << i << std::setw(2) << "|" << std::setw(10) << FarNormParams[i].index << std::setw(2) << "|" << std::setw(20) << FarNormParams[i].name << std::setw(2) << "|" << std::setw(10);
    for(unsigned j = 0; j <FarNormParams[i].modes.size(); j++){
      std::cout<< FarNormParams.at(i).modes.at(j) << " ";
    }
    std::cout << std::setw(2) << "|" << " ";
    for (unsigned j = 0; j < FarNormParams[i].targets.size(); j++) {
      std::cout << FarNormParams.at(i).targets.at(j) << " ";
    }

    std::cout << std::setw(2) << "|" << " ";

    for (unsigned j = 0; j < FarNormParams[i].pdgs.size(); j++) {
      std::cout << FarNormParams.at(i).pdgs.at(j) << " ";
    }

    std::cout << std::endl;
  }
  std::cout << std::endl;

  // FINISH Far priting


  std::cout << "Near detetor affecting spline parameters: " << nNearSplineParams << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(2) << "|" << std::setw(10) << "Global #" << std::setw(2) << "|" << std::setw(20) << "Name" << std::endl;
  for (int i = 0; i < nNearSplineParams; ++i) {
    std::cout << std::setw(4) << i << std::setw(2) << "|" << std::setw(10) << NearsplineParsIndex.at(i) << std::setw(2) << "|" << std::setw(20) << NearsplineParsNames.at(i) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "Unique spline parameters: " << nSplineParamsUniq << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(2) << "|" << std::setw(10) << "Global #" << std::setw(2) << "|" << std::setw(20) << "Name" << std::endl;
  for (int i = 0; i < nSplineParamsUniq; ++i) {
    std::cout << std::setw(4) << i << std::setw(2) << "|" << std::setw(10) << splineParsUniqIndex.at(i) << std::setw(2) << "|" << std::setw(20) << splineParsUniqNames.at(i) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "What spline parameter uses what spline?" << std::endl;
  std::cout << std::setw(10) << "Spline #" << std::setw(15) << "Uses unique #" << std::setw(2) << "|" << std::setw(20) << "Param Name" << std::setw(20) << "Uses Name" << std::endl;
  for (int i = 0; i < nNearSplineParams; ++i) {
    std::cout << std::setw(10) << i << std::setw(15) << splineParsShareToUniq.at(i) << std::setw(2) << "|" << std::setw(20) << NearsplineParsNames.at(i) << std::setw(20) << splineParsUniqNames.at(splineParsShareToUniq.at(i)) << std::endl;
  }

  std::cout << std::endl;

  std::cout << "Near affecting Functional parameters: " << nNearFuncParams << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(2) << "|" << std::setw(10) << "Global #" << std::setw(2) << "|" << std::setw(20) << "Name" << std::endl;
  for (int i = 0; i < nNearFuncParams; ++i) {
    std::cout << std::setw(4) << i << std::setw(2) << "|" << std::setw(10) << NearfuncParsIndex.at(i) << std::setw(2) << "|" << std::setw(20) << NearfuncParsNames.at(i) << std::endl;
  }

  std::cout << "\nDone printing covarianceXsec" << std::endl;
  std::cout << "#################################################" << std::endl;

  // Set precision of cout to what we started with
  std::cout.precision(ss);

} // End

// ********************************************
// Sets the proposed Flux parameters to the prior values
void covarianceXsec::setFluxOnlyParameters() {
// ********************************************
    for (int i = 0; i < size; i++) 
    {
      if(isFlux[i]) _fPropVal[i] = _fPreFitValue[i];
    }
}
// ********************************************
// Sets the proposed Flux parameters to the prior values
void covarianceXsec::setXsecOnlyParameters() {
// ********************************************
    for (int i = 0; i < size; i++) 
    {
      if(!isFlux[i]) _fPropVal[i] = _fPreFitValue[i];
    }
}
