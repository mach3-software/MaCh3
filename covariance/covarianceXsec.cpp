#include "covarianceXsec.h"

// ********************************************
covarianceXsec::covarianceXsec(const char *name, const char *file, double threshold,int firstpcapar,int lastpcapar)
  : covarianceBase(name, file,0,threshold,firstpcapar,lastpcapar) {
// ********************************************

  MakePosDef();
  TFile *infile = new TFile(file, "READ");


  xsec_param_norm_modes = NULL;
  xsec_param_norm_horncurrents = NULL;
  xsec_param_norm_elem = NULL;
  xsec_param_norm_nupdg = NULL;


  // Now to the special objects that only 2016a and above have
  if(!(xsec_param_norm_modes = (TObjArray*)(infile->Get("xsec_norm_modes")))){std::cerr<<"Can't find xec_norm_modes in xseccov"<<std::endl;throw;}
  if(!(xsec_param_norm_horncurrents = (TObjArray*)(infile->Get("xsec_norm_horncurrents")))){std::cerr<<"Can't find xec_norm_horncurrents in xseccov"<<std::endl;throw;}
  if(!(xsec_param_norm_elem  = (TObjArray*)(infile->Get("xsec_norm_elements")))){std::cerr<<"Can't find xec_norm_elements in xseccov"<<std::endl;throw;}
  if(!(xsec_param_norm_nupdg = (TObjArray*)(infile->Get("xsec_norm_nupdg")))){std::cerr<<"Can't find xec_norm_nupdg in xseccov"<<std::endl;throw;}
  if(!(xsec_param_norm_preoscnupdg = (TObjArray*)(infile->Get("xsec_norm_prod_nupdg")))){std::cerr<<"Can't find xec_norm_prod_nupdg in xseccov"<<std::endl;throw;}
  //ETA - adding in string which will then be parsed in samplePDF class into a kinematic variable to cut on
  if(!(xsec_kinematic_type = (TObjArray*)(infile->Get("xsec_norm_kinematic_type")))){std::cerr<< "[ERROR]::" << __FILE__ << ":" << __LINE__ << " cannot find xsec_kinematic_type in xseccov" << std::endl; throw;}
  if(!(xsec_param_fd_spline_modes = (TObjArray*)(infile->Get("fd_spline_modes")))){std::cerr<<"Can't find fd_spline_modes in xseccov"<<std::endl;throw;}
  if(!(xsec_param_fd_spline_names = (TObjArray*)(infile->Get("fd_spline_names")))){std::cerr<<"Can't find fd_spline_names in xseccov"<<std::endl;throw;}
  if(!(xsec_param_nd_spline_names = (TObjArray*)(infile->Get("nd_spline_names")))){std::cerr<<"Can't find nd_spline_names in xseccov"<<std::endl;throw;}
  
  // Check that the size of all the arrays are good
  if (xsec_param_norm_modes->GetEntries() != xsec_param_norm_elem->GetEntries() || 
      xsec_param_norm_modes->GetEntries() != xsec_param_norm_nupdg->GetEntries() ||
      xsec_param_norm_modes->GetEntries() != xsec_param_norm_horncurrents->GetEntries() ||
      xsec_param_norm_modes->GetEntries() != xsec_param_norm_preoscnupdg->GetEntries() ||
	  xsec_param_norm_modes->GetEntries() != xsec_kinematic_type->GetEntries() ){	
    std::cerr << "Number of entries in input matrix normalisation parameters is wrong!" << std::endl;
    std::cerr << "Element GetEntries =    " << xsec_param_norm_elem->GetEntries() << std::endl;
    std::cerr << "Modes GetEntries =      " << xsec_param_norm_modes->GetEntries() << std::endl;
    std::cerr << "Horncurrents GetEntries =      " << xsec_param_norm_horncurrents->GetEntries() << std::endl;
    std::cerr << "Nupdgs GetEntries =      " << xsec_param_norm_nupdg->GetEntries() << std::endl;
    std::cerr << "Prodnupdgs GetEntries =      " << xsec_param_norm_preoscnupdg->GetEntries() << std::endl;
	std::cerr << "Kinematic types GetEntries =       " << xsec_kinematic_type->GetEntries() << std::endl;
    throw;
  }

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
    xsec_param_nom_a[i] = (*xsec_param_nom)(i);
    // Fill the lower bound
    xsec_param_lb_a[i]  = (*xsec_param_lb)(i);
    // Fill the upper bound
    xsec_param_ub_a[i]  = (*xsec_param_ub)(i);
    // Fill the prior
    xsec_param_prior_a[i]=(*xsec_param_prior)(i);
    // Fill the names
    xsec_param_names.push_back(std::string(((TObjString*)objarr_name->At(i))->GetString()));
    // DB Fill the stepscales vector
    xsec_stepscale_vec[i] = (*xsec_stepscale)(i);    

    for (int j = 0; j < ncols; j++) {
      xsec_param_id_a[i][j] = (*xsec_param_id)(i,j);
    }

    // Also write the covarianceBase ones (grr, this is a bit of mess, try to clean it up)
    fParHiLimit[i] = xsec_param_ub_a[i];
    // Also write the covarianceBase ones (grr, this is a bit of mess, try to clean it up)
    fParLoLimit[i] = xsec_param_lb_a[i];
  } // end the for loop


  infile->Close();
  delete infile;

  // Scan through the input parameters and find which are normalisation, which are splines, and so on
  scanParameters();

  initParams(0.001);

  // Print
  Print();

  // Set the parameters
  // setParameters; why is this done in testFGD2 and not here?
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
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID
      if (GetXSecParamID(i, 0) >= 0) { //If parameter is implemented as a spline
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
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID
      if (GetXSecParamID(i, 0) >= 0) { //If parameter is implemented as a spline
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

  for (int i = 0; i < nPars; ++i) {
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID
      if (GetXSecParamID(i, 0) >= 0) { //If parameter is implemented as a spline
        returnVec.push_back(std::string(((TObjString*)xsec_param_fd_spline_names->At(i))->GetString())); //Append spline name
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

  for (int i = 0; i < nPars; ++i) {
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID
      if (GetXSecParamID(i, 0) >= 0) { //If parameter is implemented as a spline
        
	TVectorD* tempVector = (TVectorD*)(xsec_param_fd_spline_modes->At(i));
	std::vector<int> temp;
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
	returnVec.push_back(temp);	

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
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID
      if (GetXSecParamID(i, 0) >= 0) { //If parameter is implemented as a spline
	returnVec.push_back(i);
      }
    }
  }

  return returnVec;
}
// ********************************************

// ********************************************
// DB Grab the Normalisation parameters for the relevant DetID
// Eta - I think this doesn't need to be the same as scanParameters, haven't we already got this info??
const std::vector<XsecNorms4> covarianceXsec::GetNormParsFromDetID(int DetID) {
  std::vector<XsecNorms4> returnVec;
  int norm_counter = 0;

  for (int i = 0; i < nPars; ++i) {
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID   
      if (GetXSecParamID(i, 0) == -1) { //If parameter is implemented as a normalisation
	std::vector<int> temp;


	
	XsecNorms4 norm;
        norm.name=GetParameterName(i);

        // Set the mode of the normalisation parameter
        TVectorD* tempVector = (TVectorD*)(xsec_param_norm_modes->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
	  temp.push_back(tempVector[0][j]);
	}
        norm.modes=temp;
	temp.clear();

        // Set the target of the normalisation parameter
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

	//Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
	//We set a bool to see if any bounds exist so we can shortcircuit checking all of them every step
	bool haskinbounds=false;

	//IF YOU HAVE NEW KINEMATIC BOUNDS TO LOAD IN PUT THEM HERE
	// ETA - why do we have to do this twice??
	
	////////////////////
	//New generic cuts things 
	////////////////////
	
	//Only consider the kinematic string and the boundaries if you've actually given it a string to use...
	if( ((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString().Length() > 0 ){
	  std::cout << "Found Kinematic bound for parameter " << i << std::endl;
	  std::cout << "Will apply a cut on " << std::string(((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString()) << std::endl;
	  std::cout << "With lower bound " << (*xsec_kinematic_lb)(i) << " and upper bound " << (*xsec_kinematic_ub)(i) << std::endl; 
	  haskinbounds = true;

	  norm.KinematicVarStr.push_back(std::string(((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString()));

	  //ETA - This can be a vector :) can provide different kinematic variables to cut on
	  std::vector<double> single_selec;

	  //ETA - push back kinematic type with dummy -999 since this needs to be converted into an enum for a kinematic type within
	  //a samplePDFFD daughter class
	  single_selec.push_back(-999);
	  single_selec.push_back((*xsec_kinematic_lb)(i));
	  single_selec.push_back((*xsec_kinematic_ub)(i));

	  norm.Selection.push_back(single_selec);
	}

	norm.hasKinBounds=haskinbounds;
	//End of kinematic bound checking

        // Set the global parameter index of the normalisation parameter
        norm.index=i;
	//Add this parameter to the vector of parameters
	returnVec.push_back(norm);
	    norm_counter++;
      }
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
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID
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
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID
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
    if ((GetXSecParamID(i, 1) & DetID) == DetID) { //If parameter applies to required DetID
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
void covarianceXsec::scanParameters() {
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



    bool isValidDetID = false;
    int DetIDCounter = 0;
    int ParamaDetID = GetXSecParamID(i,1);
    //DB Loop over all supported DetIDs to ensure Root/XML inputs are familiar
    for (int iKnownDetID=0;iKnownDetID<MaCh3Utils::nKnownDetIDs;iKnownDetID++) {
      if ((ParamaDetID & MaCh3Utils::KnownDetIDsMap[iKnownDetID]) == MaCh3Utils::KnownDetIDsMap[iKnownDetID]) {
	isValidDetID = true;
	DetIDCounter += MaCh3Utils::KnownDetIDsMap[iKnownDetID];
      }
    }
    //DB Throw if Param DetID is unsupported. Also check that only supported DetIDs are contained in the param DetID
    if (!isValidDetID || ((ParamaDetID - DetIDCounter)!=0)) {
      std::cerr << "Unknown DetID:" << GetXSecParamID(i,1) << std::endl;
      std::cerr << "DetIDCounter:" << DetIDCounter << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
    
    if(GetParameterName(i).find("b_")==0)isFlux.push_back(true);
    else isFlux.push_back(false);

	//ETA - adding in a counter for the number of norms as xsec_norm_kinematic_type is only of length of the number of norm parameters
	//Not sure how this wasn't broken already??
	if(GetXSecParamID(i, 0) == -1){norm_counter++;} 

    // Also make some helper arrays with all the Far parameters, so we have Near parameters, Far parameters, and all parameters
    // These are parameters that have 24 or 25 in ID(i,1), so can be fancy and do a bitwise comparison to 01000
    if ((GetXSecParamID(i, 1) & 24) == 24) {//Far pars
      
      // Now check if it's a spline parameter or not
      if (GetXSecParamID(i, 0) >= 0) {//FarSplinePars
        FarSplineParsNames.push_back(GetParameterName(i));
	
	//Fill the name of the Far spline objects in the spline files
	FarSplineFileParsNames.push_back(std::string(((TObjString*)xsec_param_fd_spline_names->At(i))->GetString()));


        FarSplineParsIndex.push_back(i);
        TVectorD* tempVector = (TVectorD*)(xsec_param_fd_spline_modes->At(i));
        std::vector<int> temp;
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
	  temp.push_back(tempVector[0][j]);
	}
	FarSplineModes.push_back(temp);
	temp.clear();
	
        nFarSplineParams++;

        // Or a normalisation parameter
      } //End FarSplinePars
      else if (GetXSecParamID(i, 0) == -1) {//FarNormPars
		std::cout << "Par " << i << " is a normalisation parameter" << std::endl;
	XsecNorms4 tmp_xsec;
        tmp_xsec.name=GetParameterName(i);

        // Set the mode of the normalisation parameter
        TVectorD* tempVector = (TVectorD*)(xsec_param_norm_modes->At(i));
        std::vector<int> temp;
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
	  temp.push_back(tempVector[0][j]);
	}
        tmp_xsec.modes=temp;
	temp.clear();

        // Set the target of the normalisation parameter
        tempVector = (TVectorD*)(xsec_param_norm_horncurrents->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
        tmp_xsec.horncurrents=temp;
        temp.clear();

        // Set the target of the normalisation parameter
        tempVector = (TVectorD*)(xsec_param_norm_elem->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
        tmp_xsec.targets=temp;
        temp.clear();

	
        // Set the pdg of the normalisation parameter
        tempVector = (TVectorD*)(xsec_param_norm_nupdg->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
        tmp_xsec.pdgs=temp;
	temp.clear();


        // Set the preoscillation neutrino pdg of the normalisation parameter
        tempVector = (TVectorD*)(xsec_param_norm_preoscnupdg->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
        tmp_xsec.preoscpdgs=temp;
	temp.clear();


	//Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
	//We set a bool to see if any bounds exist so we can shortcircuit checking all of them every step
	bool haskinbounds=false;


	////////////////////
	//New generic cuts things 
	////////////////////
	

	/*
	std::cout << "############~~" << std::endl;
	std::cout << "lb is of size " << xsec_kinematic_lb->GetNrows() << std::endl;
	std::cout << "############~~" << std::endl;

	std::cout << "############~~" << std::endl;
	std::cout << "xsec_kinematic_type is of size " << xsec_kinematic_type->GetEntries() << std::endl;
	std::cout << "############~~" << std::endl;

	std::cout << "##############" << std::endl;
	std::cout << "xsec_param_norm_modes " << xsec_param_norm_modes->GetEntries() << " and we're on i: " << i << std::endl;
	std::cout << "##############" << std::endl;

	std::cout << "############~~" << std::endl;
	std::cout << "xsec_param_norm_preoscnupdf is of size " << xsec_param_norm_preoscnupdg->GetEntries() << " and we're on i: " << i << std::endl;
	std::cout << "############~~" << std::endl;
	*/


	
	//Only consider the kinematic string and the boundaries if you've actually given it a string to use...
	if( ((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString().Length() > 0){
	  std::cout << "Because xsec_kinematic_type was " << ((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString() << std::endl;
	  haskinbounds = true;
	  //ETA - This can be a vector :) can provide different kinematic variables to cut on
	  tmp_xsec.KinematicVarStr.push_back(std::string(((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString()));

	  //ETA - This can be a vector :) can provide different kinematic variables to cut on
	  std::vector<double> single_selec;

	  //ETA - push back kinematic type with dummy -999 since this needs to be converted into an enum for a kinematic type within
	  //a samplePDFFD daughter class
	  single_selec.push_back(-999);
	  single_selec.push_back((*xsec_kinematic_lb)(i));
	  single_selec.push_back((*xsec_kinematic_ub)(i));

	  tmp_xsec.Selection.push_back(single_selec);

	  std::cout << "Found kinematic bound on " << ((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString() << " from " << single_selec[1] << " to " << single_selec[2] << " for parameter " << i << std::endl;
	}
	else{
	  std::cout << "~~~" << std::endl;
	  std::cout << "Did not find kinematic bound for param " << i << std::endl;
	  std::cout << "~~~" << std::endl;
	}



	//IF YOU HAVE NEW KINEMATIC BOUNDS TO LOAD IN PUT THEM HERE

	tmp_xsec.hasKinBounds=haskinbounds;
	//End of kinematic bound checking
	

        // Set the global parameter index of the normalisation parameter
        tmp_xsec.index=i;
	//Add this parameter to the vector of parameters
	FarNormParams.push_back(tmp_xsec);

	std::cout << "FarNormParams[i].Selections has " << FarNormParams[nFarNormParams].Selection.size() << std::endl;

        nFarNormParams++;
	if (!(GetXSecParamID(i, 1) & 1)) {//non-Near affecting Far affecting norm parameters
	  nFaronlyNormParams++;
	}
      }//End FarNormPars
      else if (GetXSecParamID(i, 0) == -2) {//Far functional parameter
        nFarFuncParams++;
        FarFuncParsNames.push_back(GetParameterName(i));
        FarFuncParsIndex.push_back(i);
      }//End Far funcpars
    }//End Far pars

    //Near affecting parameters
    if ((GetXSecParamID(i, 1) & 1)) {//Near affecting parameters
      if (GetXSecParamID(i, 0) == -1) {//Near affecting norm pars
	XsecNorms4 tmp_xsec;

        tmp_xsec.name=GetParameterName(i);

        // Set the mode of the normalisation parameter
        TVectorD* tempVector = (TVectorD*)(xsec_param_norm_modes->At(i));
        std::vector<int> temp;
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
	  temp.push_back(tempVector[0][j]);
	}
        tmp_xsec.modes=temp;
	temp.clear();

        // Set the target of the normalisation parameter
        tempVector = (TVectorD*)(xsec_param_norm_horncurrents->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
        tmp_xsec.horncurrents=temp;
        temp.clear();

        // Set the target of the normalisation parameter
        tempVector = (TVectorD*)(xsec_param_norm_elem->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
        tmp_xsec.targets=temp;
        temp.clear();

	
        // Set the pdg of the normalisation parameter
        tempVector = (TVectorD*)(xsec_param_norm_nupdg->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
        tmp_xsec.pdgs=temp;
	temp.clear();


        // Set the preoscillation neutrino pdg of the normalisation parameter
        tempVector = (TVectorD*)(xsec_param_norm_preoscnupdg->At(i));
        for (int j = 0; j < tempVector->GetNrows(); ++j) {
          temp.push_back(tempVector[0][j]);
        }
        tmp_xsec.preoscpdgs=temp;
	temp.clear();


	//Next ones are kinematic bounds on where normalisation parameter should apply (at the moment only etrue but hope to add q2
	//We set a bool to see if any bounds exist so we can shortcircuit checking all of them every step
	bool haskinbounds=false;

	////////////////////
	//New generic cuts things 
	// just the same as for the Far affecting params above
	////////////////////
	
	//Only consider the kinematic string and the boundaries if you've actually given it a string to use...
	if( ((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString().Length() > 0 ){
	  haskinbounds = true;
	  //ETA - This can be a vector :) can provide different kinematic variables to cut on
	  tmp_xsec.KinematicVarStr.push_back(std::string(((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString()));

	  //ETA - This can be a vector :) can provide different kinematic variables to cut on
	  std::vector<double> single_selec;

	  //ETA - push back kinematic type with dummy -999 since this needs to be converted into an enum for a kinematic type within
	  //a samplePDFFD daughter class
	  single_selec.push_back(-999);
	  single_selec.push_back((*xsec_kinematic_lb)(i));
	  single_selec.push_back((*xsec_kinematic_ub)(i));

	  tmp_xsec.Selection.push_back(single_selec);

	  std::cout << "Found kinematic bound on " << ((TObjString*)xsec_kinematic_type->At(norm_counter))->GetString() << " from " << single_selec[1] << " to " << single_selec[2] << " for parameter " << i << std::endl;
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
      else if (GetXSecParamID(i, 0) == -2) {//Near affecting func pars
	
        NearfuncParsNames.push_back(GetParameterName(i));
        NearfuncParsIndex.push_back(i);
        nNearFuncParams++;
	
        // If they are greater >= 0 it's a spline parameter
      }//End Near affecting func pars
      else if (GetXSecParamID(i, 0) >= 0) {
        NearsplineParsNames.push_back(GetParameterName(i));
        NearsplineParsIndex.push_back(i);
        nNearSplineParams++;
	
	//Fill the name of the Far spline objects in the spline files
	NearSplineFileParsNames.push_back(std::string(((TObjString*)xsec_param_nd_spline_names->At(i))->GetString()));
      }//End Near affecting spline pars
      else {
	std::cerr << "Found a parameter in covarianceXsec which wasn't -2, -1 or above 0!" << std::endl;
	std::cerr << "This is undefined behaviour currently, and implementation should change" << std::endl;
	std::cerr << "Param " << GetParameterName(i) << " (param " << i << ") = " << GetXSecParamID(i, 0) << std::endl;
	std::cerr << __LINE__ << ":" << __LINE__ << std::endl;
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
// Throw the nominal values according to the cholesky decomposed
void covarianceXsec::throwNominal(bool nomValues, int seed) {
  // ********************************************

  TDecompChol chdcmp(*covMatrix);

  if(!chdcmp.Decompose()) {
    std::cerr << "Cholesky decomposition failed for " << matrixName << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  chel = new TMatrixD(chdcmp.GetU());
  CholeskyDecomp(size, *chel);

  ThrowParms* nom_throws = new ThrowParms((*xsec_param_prior), (*covMatrix));
  nom_throws->SetSeed(seed);
  nominal.clear();
  nominal.resize(size);

  // default is nomValues = true
  if(!nomValues) {
    // Need to throw again if we throw the nominal outside the upper or lower bound
    bool throw_again = true;
    while(throw_again == true) {
      throw_again = false;
      nom_throws->ThrowSet(nominal);
      for (int i = 0; i < int(nominal.size()); i++) {
        // If the parameter is fixed, set to the prior
        if (fParSigma[i] < 0) {
          nominal[i] = (*xsec_param_prior)(i);
          continue;
        }
        // If the throw is outside the lower bound
        if (nominal[i] < (*xsec_param_lb)(i)) throw_again = true;
        // If the throw is outside the upper bound
        if (nominal[i] > (*xsec_param_ub)(i)) throw_again = true;
      }
    }
  } else {
    // Set nominals to prior
    for (int i = 0; i < int(nominal.size()); i++) {
      nominal[i] = (*xsec_param_prior)(i);
      //nominal[i] = (*xsec_param_nom)(i);
    }
  }
  delete nom_throws;
}

// ********************************************
void covarianceXsec::initParams(double fScale) {
  // ********************************************

  for (int i = 0; i < size; ++i) {
    char pname[10];
    fParNames[i] = new Char_t[9];
    sprintf(pname, "xsec_%i", i);
    strcpy(fParNames[i], pname);

    // This just assign the initial parameters to the prior
    fParInit[i] = (*xsec_param_prior)[i];
    //fParInit[i] = (*xsec_param_nom)[i];

    // Any param with nom == 1 should be > 0
    if ((*xsec_param_nom)(i) == 1) {
      // If the fParInit is negative we should try to throw it above this
      // We don't really do this ever...
      while (fParInit[i] <= 0) {
        fParInit[i] = random_number[0]->Gaus((*xsec_param_prior)[i], 0.1*fScale*TMath::Sqrt( (*covMatrix)(i,i) ));
      }
    }

    // Set covarianceBase parameters (Curr = current, Prop = proposed, Sigma = step)
    fParCurr[i] = fParInit[i];
    fParProp[i] = fParCurr[i];
    fParSigma[i] = fScale;
  }

  //DB Set Individual Step scale for PCA parameters to the lastpcadpar fIndivStepScale because the step scale for those parameters is set by 'eigen_values[i]' but needs an overall step scale
  //   However, individual step scale for non-PCA parameters needs to be set correctly
  if (pca) {
    for (int i=firstpcadpar;i<=lastpcadpar;i++) {
      xsec_stepscale_vec[i] = xsec_stepscale_vec[lastpcadpar-1];
    }
  }

  setIndivStepScale(xsec_stepscale_vec);
  genPropKernels();
  randomize();
  throwNominal();
  CorrelateSteps();
}

// ********************************************
// Get the likelihood for moving the cross-section parameters to these values
double covarianceXsec::getLikelihood() {
  // ********************************************

  double xsecLogL = 0.0;

//KS: This brings speed up of the order 2 per thread, since xsec matix can only grow this might be even more profitable in the future
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:xsecLogL)
#endif
  for (int i = 0; i < nPars; i++) {

    // Make sure we're in a good region for parameter i
    if (fParProp[i] > xsec_param_ub_a[i] || fParProp[i] < xsec_param_lb_a[i]) {
      xsecLogL += __LARGE_LOGL__;
    }

    // Apply the prior
    for (int j = 0; j <= i; ++j) {
      // If we want to evaluate likelihood (Gaussian prior essentially, not flat)
      if (fParEvalLikelihood[i] && fParEvalLikelihood[j]) {
        //KS: Since matrix is symetric we can calcaute non daigonal elements only once and multiply by 2, can bring up to factor speed decrease.   
        int scale = 1;
        if(i != j) scale = 2;
        xsecLogL += scale * 0.5*( nominal[i]-fParProp[i])*(nominal[j]-fParProp[j])*(*invCovMatrix)(i,j);
      }
    } // end j for loop
  } // end i for loop


  return xsecLogL;
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
    fParEvalLikelihood[i] = eL;
  }
}

// ********************************************
void covarianceXsec::toggleFixParameter(int i) {
  // ********************************************

  if (i > size) {
    std::cerr << "Can't toggleFixParameter for parameter " << i << " because size of covariance =" << size << std::endl;
    std::cerr << "Fix this in your config file please!" << std::endl;
    exit(-1);
  } else {
    fParSigma[i] *= -1.0;
    std::cout << "Setting " << GetParameterName(i) << " (parameter " << i << ") to fixed at " << fParCurr[i] << std::endl;
  }
}

// ********************************************
void covarianceXsec::setXsecParNames() {
  // ********************************************
  // This feels sort of silly but fine...
  // Shouldn't really be called because a lot of post-processing depends on having name xsec_i
  // Have made covarianceXsec->GetParameterName(i) which returns the cross-section parameter name which is read from the ROOT input file (e.g. xsec_covariance_2015v0.root)
  for (int i = 0; i < size; ++i) {
    fParNames[i] = const_cast<char*>(xsec_param_names[i].c_str());
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
    std::cout << std::left << std::setprecision(3) << std::setw(5) << i << std::setw(2) << "|" << std::setw(25) << GetParameterName(i) << std::setw(2) << "|" << std::setw(10) << (*xsec_param_nom)(i) << std::setw(2) << "|" << std::setw(10) << (*xsec_param_prior)(i) << std::setw(2) << "|" << "+/- " << std::setw(11) << sqrt((*covMatrix)(i,i)) << std::setw(2) << "|" << std::setw(10) << (*xsec_param_lb)(i) << std::setw(2) << "|" << std::setw(10) << (*xsec_param_ub)(i) << "|" << std::setw(15) << fIndivStepScale[i] << std::endl;
  }

  std::cout << std::endl;

  // Output the normalisation parameters as a sanity check!
  std::cout << "Normalisation parameters:" << nNearNormParams << std::endl;
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

  // Start Far priting


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
// Sets the proposed Flux parameters to the nominal values
void covarianceXsec::setFluxOnlyParameters() {
// ********************************************
    for (int i = 0; i < size; i++) 
    {
      if(isFlux[i]) fParProp[i] = nominal[i];
    }
}
// ********************************************
// Sets the proposed Flux parameters to the nominal values
void covarianceXsec::setXsecOnlyParameters() {
// ********************************************
    for (int i = 0; i < size; i++) 
    {
      if(!isFlux[i]) fParProp[i] = nominal[i];
    }
}
