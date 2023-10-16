#include "covarianceXsec.h"

// ********************************************
covarianceXsec::covarianceXsec(const char *name, const char *file,
   	                           double threshold,int FirstPCAdpar,
							   int LastPCAdpar)
  : covarianceBase(name, file,0,threshold,FirstPCAdpar,LastPCAdpar) {
// ********************************************

  /*
   * ETA - this is all redundent now
  TFile *infile = new TFile(file, "READ");

  xsec_param_norm_modes = NULL;
  xsec_param_norm_horncurrents = NULL;
  xsec_param_norm_elem = NULL;
  xsec_param_norm_nupdg = NULL;

  // Now to the special objects that only 2016a and above have
  if(!(xsec_param_norm_modes = (TObjArray*)(infile->Get("xsec_norm_modes")))){
	std::cerr<<"Can't find xec_norm_modes in xseccov"<<std::endl;
	throw;
  }
  if(!(xsec_param_norm_horncurrents = (TObjArray*)(infile->Get("xsec_norm_horncurrents")))){
	std::cerr<<"Can't find xec_norm_horncurrents in xseccov"<<std::endl;
	throw;
  }
  if(!(xsec_param_norm_elem  = (TObjArray*)(infile->Get("xsec_norm_elements")))){std::cerr<<"Can't find xec_norm_elements in xseccov"<<std::endl;throw;}
  if(!(xsec_param_norm_nupdg = (TObjArray*)(infile->Get("xsec_norm_nupdg")))){std::cerr<<"Can't find xec_norm_nupdg in xseccov"<<std::endl;throw;}
  if(!(xsec_param_norm_preoscnupdg = (TObjArray*)(infile->Get("xsec_norm_prod_nupdg")))){std::cerr<<"Can't find xec_norm_prod_nupdg in xseccov"<<std::endl;throw;}
  //ETA - adding in string which will then be parsed in samplePDF class into a kinematic variable to cut on
  //LW - Changing to if no kinematic type, make empty array
  //if(!(xsec_kinematic_type = (TObjArray*)(infile->Get("xsec_norm_kinematic_type")))){std::cerr<< "[ERROR]::" << __FILE__ << ":" << __LINE__ << " cannot find xsec_kinematic_type in xseccov" << std::endl; throw;}
  if(!(xsec_param_fd_spline_modes = (TObjArray*)(infile->Get("fd_spline_modes")))){std::cerr<<"Can't find fd_spline_modes in xseccov"<<std::endl;throw;}
  if(!(xsec_param_fd_spline_names = (TObjArray*)(infile->Get("fd_spline_names")))){std::cerr<<"Can't find fd_spline_names in xseccov"<<std::endl;throw;}
  //if(!(xsec_param_nd_spline_names = (TObjArray*)(infile->Get("nd_spline_names")))){std::cerr<<"Can't find nd_spline_names in xseccov"<<std::endl;throw;}
  
  if(!(xsec_kinematic_type = (TObjArray*)(infile->Get("xsec_norm_kinematic_type")))){xsec_kinematic_type = new TObjArray();}
  if(!(xsec_param_nd_spline_names = (TObjArray*)(infile->Get("nd_spline_names")))){xsec_param_nd_spline_names= new TObjArray();}


  // Check that the size of all the arrays are good
  if (xsec_param_norm_modes->GetEntries() != xsec_param_norm_elem->GetEntries() || 
      xsec_param_norm_modes->GetEntries() != xsec_param_norm_nupdg->GetEntries() ||
      xsec_param_norm_modes->GetEntries() != xsec_param_norm_horncurrents->GetEntries() ||
      xsec_param_norm_modes->GetEntries() != xsec_param_norm_preoscnupdg->GetEntries() *//* ||
	   xsec_param_norm_modes->GetEntries() != xsec_kinematic_type->GetEntries() *//* ){	
    std::cerr << "Number of entries in input matrix normalisation parameters is wrong!" << std::endl;
    std::cerr << "Element GetEntries =    " << xsec_param_norm_elem->GetEntries() << std::endl;
    std::cerr << "Modes GetEntries =      " << xsec_param_norm_modes->GetEntries() << std::endl;
    std::cerr << "Horncurrents GetEntries =      " << xsec_param_norm_horncurrents->GetEntries() << std::endl;
    std::cerr << "Nupdgs GetEntries =      " << xsec_param_norm_nupdg->GetEntries() << std::endl;
    std::cerr << "Prodnupdgs GetEntries =      " << xsec_param_norm_preoscnupdg->GetEntries() << std::endl;
	std::cerr << "Kinematic types GetEntries =       " << xsec_kinematic_type->GetEntries() << std::endl;
    throw;
  }

  //ETA- don't need this explicit difference.
  // I think we just keep one vector for the spline modes. This could get a bit
  // weird in the case that there are different modes across systs that you
  // want to correlate. But I think this will be fine.
  if(xsec_param_fd_spline_modes->GetEntries() != xsec_param_fd_spline_names->GetEntries() ){
    std::cerr <<"Number of entries in input matrix fd spline parameters is wrong!"<<std::endl;
    std::cerr <<"Far spline modes GetEntries = "<<xsec_param_fd_spline_modes->GetEntries()<<std::endl;
    std::cerr <<"Far spline names GetEntries = "<<xsec_param_fd_spline_names->GetEntries()<<std::endl;
    throw;
  }

  // set extra xsec objects
  xsec_param_nom    = (TVectorD*)(infile->Get("xsec_param_nom"));
  xsec_param_prior  = (TVectorD*)(infile->Get("xsec_param_prior"));
  xsec_param_id     = (TMatrixD*)(infile->Get("xsec_param_id"));
  xsec_param_lb     = (TVectorD*)(infile->Get("xsec_param_lb"));
  xsec_param_ub     = (TVectorD*)(infile->Get("xsec_param_ub"));
  xsec_stepscale    = (TVectorD*)(infile->Get("xsec_stepscale"));
  TVectorD* flat_prior = (TVectorD*)(infile->Get("xsec_flat_prior"));
  xsec_kinematic_ub = (TVectorD*)(infile->Get("xsec_norm_kinematic_ub"));
  xsec_kinematic_lb = (TVectorD*)(infile->Get("xsec_norm_kinematic_lb"));
  TObjArray* objarr_name = (TObjArray*)(infile->Get("xsec_param_names"));

  nPars = xsec_param_prior->GetNrows();

  int ncols = xsec_param_id->GetNcols();

  // Check that the input matrix all makes sense in terms of size
  if (xsec_param_nom->GetNrows() != nPars) {
    std::cerr << "Reading " << file << ":" << name << std::endl;
    std::cerr << "Number of rows of nominal != number of rows of prior" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  */

  int ncols = _fNumPar;
  //ParseYAML(name);

  nPars = _fNumPar;
  // Leave these as arrays for backwards compatibility
  // Probably should be made into vectors but at least size isn't hard-coded ^,^
  // Load up the parameters id (nPars x nCols big, contains information about if these parameters are spline, normalisation or some other parameter)
  xsec_param_id_a   = new int*[nPars];
  for (int i = 0; i < nPars; ++i) {
    xsec_param_id_a[i] = new int[ncols];
  }
  xsec_param_nom_a  = new double[nPars];
  xsec_param_lb_a   = new double[nPars];
  xsec_param_ub_a   = new double[nPars];
  xsec_param_prior_a = new double[nPars];

  //DB Resize step scale vector to expected size
  xsec_stepscale_vec.resize(nPars);

  for (int i = 0; i < nPars; i++) {
    // Fill the prior central value
	xsec_param_nom_a[i] = _fPreFitValue[i];

    // Fill the lower bound
    xsec_param_lb_a[i]  = _fLowBound[i];

    // Fill the upper bound
    xsec_param_ub_a[i]  = _fUpBound[i];
    // Fill the prior uncertainty
	xsec_param_prior_a[i]= _fError[i];
    // Fill the names
    xsec_param_names.push_back(_fNames[i]);

    if(xsec_param_names[i].length() > PrintLength) PrintLength = xsec_param_names.back().length();
    // DB Fill the stepscales vector
    xsec_stepscale_vec[i] = _fIndivStepScale[i];

    if(_fFlatPrior[i]){setEvalLikelihood(i, true);}
    
    for (int j = 0; j < ncols; j++) {
      xsec_param_id_a[i][j] = (*xsec_param_id)(i,j);
    }

    // Also write the covarianceBase ones (grr, this is a bit of mess, try to clean it up)
    _fUpBound[i] = xsec_param_ub_a[i];
    // Also write the covarianceBase ones (grr, this is a bit of mess, try to clean it up)
    _fLowBound[i] = xsec_param_lb_a[i];
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

  //int ncols = _fNumPar;
  ParseYAML(YAMLFile);

  int ncols = _fNumPar;
  nPars = _fNumPar;
  // Leave these as arrays for backwards compatibility
  // Probably should be made into vectors but at least size isn't hard-coded ^,^
  // Load up the parameters id (nPars x nCols big, contains information about if these parameters are spline, normalisation or some other parameter)
  xsec_param_id_a   = new int*[nPars];
  for (int i = 0; i < nPars; ++i) {
    xsec_param_id_a[i] = new int[ncols];
  }
  xsec_param_nom_a  = new double[nPars];
  xsec_param_lb_a   = new double[nPars];
  xsec_param_ub_a   = new double[nPars];
  xsec_param_prior_a = new double[nPars];

  //DB Resize step scale vector to expected size
  xsec_stepscale_vec.resize(nPars);

  for (int i = 0; i < nPars; i++) {
    // Fill the nominal
	xsec_param_nom_a[i] = _fPreFitValue[i];

    // Fill the lower bound
    xsec_param_lb_a[i]  = _fLowBound[i];

    // Fill the upper bound
    xsec_param_ub_a[i]  = _fUpBound[i];

    // Fill the prior
	xsec_param_prior_a[i]= _fError[i];

    // Fill the names
    xsec_param_names.push_back(_fNames[i]);

    if(xsec_param_names[i].length() > PrintLength) PrintLength = xsec_param_names.back().length();

    // DB Fill the stepscales vector
    xsec_stepscale_vec[i] = _fIndivStepScale[i];

    //if((*flat_prior)(i)) setEvalLikelihood(i, false);
    if(!_fFlatPrior[i]){setEvalLikelihood(i, false);}
    
	//xsec_param_id_a[i][0] is the "type" of systematic, spline or otherwise
	//xsec_param_id_a[i][1] is the DetId
    //for (int j = 0; j < ncols; j++) {
    //  xsec_param_id_a[i][j] = (*xsec_param_id)(i,j);
    //}	

	//ETA - does this need to be here still?
    // Also write the covarianceBase ones (grr, this is a bit of mess, try to clean it up)
    _fUpBound[i] = xsec_param_ub_a[i];
    // Also write the covarianceBase ones (grr, this is a bit of mess, try to clean it up)
    _fLowBound[i] = xsec_param_lb_a[i];
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

  for (auto const &param : _fYAMLDoc["Systematics"]) {
     std::cout << param["Systematic"]["Names"]["ParameterName"].as<std::string>() << std::endl;

     _fNames[i] = (param["Systematic"]["Names"]["ParameterName"].as<std::string>());
     _fPreFitValue[i] = (param["Systematic"]["ParameterValues"]["PreFitValue"].as<double>());
     _fGenerated[i] = (param["Systematic"]["ParameterValues"]["Generated"].as<double>());

	 //ETA - a bit of a fudge but works
	 std::vector<double> TempBoundsVec = param["Systematic"]["ParameterBounds"].as<std::vector<double>>();
     _fLowBound[i] = TempBoundsVec[0];
     _fUpBound[i] = TempBoundsVec[1];


     _fIndivStepScale[i] = (param["Systematic"]["StepScale"]["MCMC"].as<double>());
     _fDetID[i] = (param["Systematic"]["DetID"].as<int>());
     _fError[i] = (param["Systematic"]["Error"].as<double>());
	 _fParamType[i] = (param["Systematic"]["Type"].as<std::string>());

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
	   std::cout << "Reading in a Spline Parameter" << std::endl;

	   if (param["Systematic"]["SplineInformation"]["FDSplineName"]) {
		 _fFDSplineNames.push_back(param["Systematic"]["SplineInformation"]["FDSplineName"].as<std::string>());
		 nFDSplines++;
	   }

	   if (param["Systematic"]["SplineInformation"]["FDMode"]) {
		 std::cout << "Pushing back _fFDSplineModes for param " << i << std::endl;
		 _fFDSplineModes.push_back(param["Systematic"]["SplineInformation"]["FDMode"].as<std::vector<int>>());
		 std::cout << "_fFDSplineModes is of size " << _fFDSplineModes.size() << std::endl;
		 std::cout << "_fFDSplineModes[0] is of size " << _fFDSplineModes[0].size() << std::endl;
	   }

	   if (param["Systematic"]["NDSplineName"]) {
		 _fNDSplineNames.push_back(param["Systematic"]["NDSplineName"].as<std::string>());
	   }

	 } else if(param["Systematic"]["Type"].as<std::string>() == "Norm") {
	   std::cout << "Norm parameter" << std::endl;

	   //First check to see if we have specified a mode
	   std::vector<int> DummyModeVec;
	   if(param["Systematic"]["Mode"]){
		 _fNormModes.push_back(param["Systematic"]["Mode"].as<std::vector<int>>());
	   } else{
		 //Has to be of size 0 to mean apply to all
		 _fNormModes.push_back(DummyModeVec);
	   }
	 }

	 int NumKinematicCuts = 0;
	 if(param["Systematic"]["KinematicCuts"]){

	   NumKinematicCuts = param["Systematic"]["KinematicCuts"].size();
	   std::cout << "Number of Kinematic cuts is " << NumKinematicCuts << std::endl;

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
	   std::cout << "FOUND " << param["Systematic"]["Correlations"].size() << " CORRELATIONS!" << std::endl;
	   for(int Corr_i = 0 ; Corr_i < param["Systematic"]["Correlations"].size()  ; ++Corr_i){
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
const int covarianceXsec::GetNumSplineParamsFromDetID(int DetID) {
  int returnVal = 0; 
  for (int i = 0; i < nPars; ++i) {
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
  for (int i = 0; i < nPars; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
        returnVec.push_back(GetParameterName(i));
      }
    }
  }

  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the Spline Names for the relevant DetID
const std::vector<std::string> covarianceXsec::GetSplineFileParsNamesFromDetID(int DetID) {
  std::vector<std::string> returnVec;

  int FDSplineCounter = 0;
  for (int i = 0; i < nPars; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (strcmp(GetXsecParamType(i), "Spline") == 0) { //If parameter is implemented as a spline
        returnVec.push_back(_fFDSplineNames[FDSplineCounter]); //Append spline name
		FDSplineCounter++;
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
  for (int i = 0; i < nPars; ++i) {
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

  for (int i = 0; i < nPars; ++i) {
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

  for (int i = 0; i < nPars; ++i) {
	if (strcmp(GetXsecParamType(i), "Norm") == 0) { //If parameter is implemented as a normalisation

	  if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID   

		std::vector<int> temp;

		XsecNorms4 norm;
		norm.name=GetParameterName(i);

		// Set the mode of the normalisation parameter
		//TVectorD* tempVector = (TVectorD*)(xsec_param_norm_modes->At(i));
		//for (int j = 0; j < tempVector->GetNrows(); ++j) {
		//  temp.push_back(tempVector[0][j]);
		//}
		//norm.modes=temp;
		//temp.clear();

		//Copy the mode information into an XsecNorms4 struct
		norm.modes = _fNormModes[norm_counter];

		// Set the target of the normalisation parameter
		// ETA - we can remove these from XsecNorms4 now I think
		/*
		tempVector = (TVectorD*)(xsec_param_norm_horncurrents->At(i));
		for (int j = 0; j < tempVector->GetNrows(); ++j) {
		  temp.push_back(tempVector[0][j]);
		}
		norm.horncurrents=temp;
		temp.clear();

		// Set the target of the normalisation parameter
		tempVector = (TVectorD*)(xsec_param_norm_elem->At(i));
		for (int j = 0; j < tempVector->GetNrows(); ++j) {
		  temp.push_back(tempVector[0][j]);
		}
		norm.targets=temp;
		temp.clear();

		// Set the pdg of the normalisation parameter
		tempVector = (TVectorD*)(xsec_param_norm_nupdg->At(i));
		for (int j = 0; j < tempVector->GetNrows(); ++j) {
		  temp.push_back(tempVector[0][j]);
		}
		norm.pdgs=temp;
		temp.clear();

		// Set the preoscillation neutrino pdg of the normalisation parameter
		tempVector = (TVectorD*)(xsec_param_norm_preoscnupdg->At(i));
		for (int j = 0; j < tempVector->GetNrows(); ++j) {
		  temp.push_back(tempVector[0][j]);
		}
		norm.preoscpdgs=temp;
		temp.clear();
		*/

		//Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
		//We set a bool to see if any bounds exist so we can shortcircuit checking all of them every step
		bool HasKinBounds=false;

		////////////////////
		//New generic cuts things 
		////////////////////

		///
		if(_fKinematicPars.at(i).size() > 0){
		  HasKinBounds = true;
		}

		for(int KinematicCut_i = 0 ; KinematicCut_i < _fKinematicPars[i].size() ; ++KinematicCut_i){
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
const int covarianceXsec::GetNumFuncParamsFromDetID(int DetID) {
  int returnVal = 0;

  for (int i = 0; i < nPars; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (GetXSecParamID(i, 0) == -2) { //If parameter is implemented as a functional parameter
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

  for (int i = 0; i < nPars; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (GetXSecParamID(i, 0) == -2) { //If parameter is implemented as a functional parameter
		returnVec.push_back(GetParameterName(i));
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

  for (int i = 0; i < nPars; ++i) {
    if ((GetXsecParamDetID(i) & DetID) == DetID) { //If parameter applies to required DetID
      if (GetXSecParamID(i, 0) == -2) { //If parameter is implemented as a functional parameter
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

  for (int i = 0; i < nPars; ++i) {

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
    
    if(GetParameterName(i).find("b_")==0)isFlux.push_back(true);
    else isFlux.push_back(false);

	//ETA - adding in a counter for the number of norms as xsec_norm_kinematic_type is only of length of the number of norm parameters
	//Not sure how this wasn't broken already??
	//This needs to be updated to check against a string
	if(strcmp(GetXsecParamType(i), "Norm") == 0){norm_counter++;} 

    // Also make some helper arrays with all the Far parameters, so we have Near parameters, Far parameters, and all parameters
    // These are parameters that have 24 or 25 in ID(i,1), so can be fancy and do a bitwise comparison to 01000
	if ((GetXsecParamDetID(i) & 24) == 24) {//Far pars

	  // Now check if it's a spline parameter or not
	  //This needs to be updated to check against a string
	  if (strcmp(GetXsecParamType(i), "Spline") == 0) {//FarSplinePars
		std::cout << GetXsecParamType(i) << std::endl;
		FarSplineParsNames.push_back(GetParameterName(i));
		FarSplineParsIndex.push_back(i);

		for(int Mode_i = 0 ; Mode_i < _fFDSplineModes[nFarSplineParams].size() ; ++Mode_i){
		  std::cout << "Mode to apply to is " << _fFDSplineModes[nFarSplineParams][Mode_i] << std::endl;
		}

		FarSplineModes.push_back(_fFDSplineModes[nFarSplineParams]);
		nFarSplineParams++;

		// Or a normalisation parameter
	  } //End FarSplinePars
	  else if (strcmp(GetXsecParamType(i), "Norm") == 0) {//FarNormPars
		XsecNorms4 tmp_xsec;
		tmp_xsec.name=GetParameterName(i);

		tmp_xsec.modes=_fNormModes[i];

		/*
		// Set the mode of the normalisation parameter
		TVectorD* tempVector = (TVectorD*)(xsec_param_norm_modes->At(i));
		std::vector<int> temp;
		for (int j = 0; j < tempVector->GetNrows(); ++j) {
		  temp.push_back(tempVector[0][j]);
		}
		tmp_xsec.modes=temp;
		temp.clear();
		*/

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
		  for(int KinVar_i = 0 ; KinVar_i < _fKinematicPars[i].size() ; ++KinVar_i) {
			Selections[KinVar_i].push_back(-999.9);
			Selections[KinVar_i].push_back(_fKinematicBounds[i][KinVar_i][0]);
			Selections[KinVar_i].push_back(_fKinematicBounds[i][KinVar_i][1]);
			std::cout << "  - " << _fKinematicPars[i][KinVar_i] << " from " << Selections[KinVar_i][1] << " to " << Selections[KinVar_i][2] << std::endl;
		  }

		  tmp_xsec.Selection = Selections;
		}

		tmp_xsec.hasKinBounds=haskinbounds;
		//End of kinematic bound checking	

		// Set the global parameter index of the normalisation parameter
		tmp_xsec.index=i;
		//Add this parameter to the vector of parameters
		FarNormParams.push_back(tmp_xsec);
		nFarNormParams++;

		if (!(GetXsecParamDetID(i) & 1)) {//non-Near affecting Far affecting norm parameters
		  nFaronlyNormParams++;
		}
	  }//End FarNormPars
	  else if (strcmp(GetXsecParamType(i), "functional") != 0){//Far functional parameter
		nFarFuncParams++;
		FarFuncParsNames.push_back(GetParameterName(i));
		FarFuncParsIndex.push_back(i);
	  }//End Far funcpars
	}//End Far pars

    //Near affecting parameters
	if ((GetXsecParamDetID(i) & 1)) {//Near affecting parameters
	  //Check if the parameter is normalisation parameter
	  if (strcmp(GetXsecParamType(i), "Norm") == 0){//Near affecting norm pars
		XsecNorms4 tmp_xsec;

		tmp_xsec.name=GetParameterName(i);

		tmp_xsec.modes=_fNormModes[i];
		//temp.clear();

		//ETA - these are no longer needed as stored as kinematic cuts
		// need to think about how to add these into XsecNorms4 though.
		// We could do it here as we have the kinematic cut information
		// already stored.
		// Set the target of the normalisation parameter
		//tempVector = (TVectorD*)(xsec_param_norm_horncurrents->At(i));
		//for (int j = 0; j < tempVector->GetNrows(); ++j) {
		//  temp.push_back(tempVector[0][j]);
		//}
		//tmp_xsec.horncurrents=temp;
		//temp.clear();

		// Set the target of the normalisation parameter
		//tempVector = (TVectorD*)(xsec_param_norm_elem->At(i));
		//for (int j = 0; j < tempVector->GetNrows(); ++j) {
		//  temp.push_back(tempVector[0][j]);
		//}
		//tmp_xsec.targets=temp;
		//temp.clear();

		// Set the pdg of the normalisation parameter
		//tempVector = (TVectorD*)(xsec_param_norm_nupdg->At(i));
		//for (int j = 0; j < tempVector->GetNrows(); ++j) {
		//  temp.push_back(tempVector[0][j]);
		//}
		//tmp_xsec.pdgs=temp;
		//temp.clear();

		//// Set the preoscillation neutrino pdg of the normalisation parameter
		//tempVector = (TVectorD*)(xsec_param_norm_preoscnupdg->At(i));
		//for (int j = 0; j < tempVector->GetNrows(); ++j) {
		//  temp.push_back(tempVector[0][j]);
		//}
		//tmp_xsec.preoscpdgs=temp;
		//temp.clear();

		//Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
		//We set a bool to see if any bounds exist so we can shortcircuit checking all of them every step
		bool haskinbounds=false;

		////////////////////
		//New generic cuts things 
		// just the same as for the Far affecting params above
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
		  for(int KinVar_i = 0 ; KinVar_i < _fKinematicPars[i].size() ; ++KinVar_i) {
			  Selections[KinVar_i].push_back(-999.9);
			  Selections[KinVar_i].push_back(_fKinematicBounds[i][KinVar_i][0]);
			  Selections[KinVar_i].push_back(_fKinematicBounds[i][KinVar_i][1]);
		  }
		  tmp_xsec.Selection = Selections;
		}
		else{std::cout << "Did not find kinematic bound for param " << i << std::endl;}

		//IF YOU HAVE NEW KINEMATIC BOUNDS TO LOAD IN PUT THEM HERE

		tmp_xsec.hasKinBounds=haskinbounds;
		//End of kinematic bound checking

		// Set the global parameter index of the normalisation parameter
		tmp_xsec.index=i;

		//IF YOU HAVE NEW KINEMATIC BOUNDS TO LOAD IN PUT THEM HERE

		// Set the global parameter index of the normalisation parameter
		tmp_xsec.index=i;
		//Add this parameter to the vector of parameters
		NearNormParams.push_back(tmp_xsec);

		// Count how many normalisation parameters we have
		nNearNormParams++;

		// Function parameters, such as BeRPA
	  }//End Near affecting normpars
	  else if (strcmp(GetXsecParamType(i), "functional") != 0) {//Near affecting func pars

		NearfuncParsNames.push_back(GetParameterName(i));
		NearfuncParsIndex.push_back(i);
		nNearFuncParams++;

		// If they are greater >= 0 it's a spline parameter
	  }//End Near affecting func pars
	  else if (strcmp(GetXsecParamType(i), "Spline") != 0) {
		NearsplineParsNames.push_back(GetParameterName(i));
		NearsplineParsIndex.push_back(i);
		nNearSplineParams++;

		//Fill the name of the Far spline objects in the spline files
		//NearSplineFileParsNames.push_back(std::string(((TObjString*)xsec_param_nd_spline_names->At(i))->GetString()));
	  }//End Near affecting spline pars
	  else {
		std::cerr << "Found a parameter in covarianceXsec which wasn't -2, -1 or above 0!" << std::endl;
		std::cerr << "This is undefined behaviour currently, and implementation should change" << std::endl;
		std::cerr << "Param " << GetParameterName(i) << " (param " << i << ") = " << GetXsecParamType(i) << std::endl;
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		//throw;
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
    _fPreFitValue[i] = _fPreFitValue[i];//(*xsec_param_prior)[i];
    //_fPreFitValue[i] = (*xsec_param_nom)[i];

    // Any param with nom == 1 should be > 0
    if (_fPreFitValue[i] == 1) {
      // If the _fPreFitValue is negative we should try to throw it above this
      // We don't really do this ever...
      while (_fPreFitValue[i] <= 0) {
        _fPreFitValue[i] = random_number[0]->Gaus((*xsec_param_prior)[i], 0.1*fScale*TMath::Sqrt( (*covMatrix)(i,i) ));
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
void covarianceXsec::setEvalLikelihood(int i, bool eL) {
  // ********************************************

  if (i > size) {
    std::cerr << "Can't setEvalLikelihood for xsec_" << i << " (" << GetParameterName(i) << ")" << " because size of covarianceXsec = " << size << std::endl;
    std::cerr << "Fix this in your config file please!" << std::endl;
    throw;
  } else {
    std::cout << "Setting " << GetParameterName(i) << " (parameter " << i << ") to flat prior" << std::endl;
    _fFlatPrior[i] = eL;
  }
}

// ********************************************
void covarianceXsec::toggleFixParameter(int i) {
  // ********************************************

  if(!pca){
	if (i > size) {
	  std::cerr << "Can't toggleFixParameter for parameter " << i << " because size of covariance =" << size << std::endl;
	  std::cerr << "Fix this in your config file please!" << std::endl;
	  exit(-1);
	} else {
	  _fError[i] *= -1.0;
	  std::cout << "Setting " << GetParameterName(i) << " (parameter " << i << ") to fixed at " << _fCurrVal[i] << std::endl;

	}
  } else {
	//KS: Find xsec parameter  in PCA base
	int isDecom = -1;
	for (int im = 0; im < npars; ++im) { if(isDecomposed_PCA[im] == i) isDecom = im; }

	if(isDecom < 0) {
	  std::cerr << "Parameter " << GetParameterName(i) << " is PCA decomposed can't fix this" << std::endl;
	  //throw; 
	} else {
	  fParSigma_PCA[isDecom] *= -1.0;
	  std::cout << "Setting un-decomposed " << getParName(i) << "(parameter " << i <<"/"<< isDecom<< " in PCA base) to fixed at " << _fCurrVal[i] << std::endl;
	}
  }
  return;
}

// ********************************************
void covarianceXsec::setXsecParNames() {
  // ********************************************
  // This feels sort of silly but fine...
  // Shouldn't really be called because a lot of post-processing depends on having name xsec_i
  // Have made covarianceXsec->GetParameterName(i) which returns the cross-section parameter name which is read from the ROOT input file (e.g. xsec_covariance_2015v0.root)
  for (int i = 0; i < size; ++i) {
	_fNames[i] = xsec_param_names[i]; 
  }
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
    std::cout << std::left << std::setprecision(3) << std::setw(5) << i << std::setw(2) << "|" << std::setw(25) << GetParameterName(i) << std::setw(2) << "|" << std::setw(10) << _fPreFitValue[i] << std::setw(2) << "|" << std::setw(10) << _fPreFitValue[i] << std::setw(2) << "|" << "+/- " << std::setw(11) << _fError[i] << std::setw(2) << "|" << std::setw(10) << _fLowBound[i] << std::setw(2) << "|" << std::setw(10) << _fUpBound[i] << "|" << std::setw(15) << _fIndivStepScale[i] << std::endl;
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
