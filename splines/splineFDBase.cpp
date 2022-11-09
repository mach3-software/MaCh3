#include <iostream>
#include <stdlib.h>
#include "splineFDBase.h"
#include "TF1.h"

//#define DEBUG_WEIGHTS


// ----- Constructor (first: original erec version, then 2d version) ---- //

splineFDBase::splineFDBase(const char *name, int ntype, int nevents, int DetID, covarianceXsec* xsec_cov) // constructor for erec spline binning
  : splineBase(name, ntype)
{

  if(xsec_cov == NULL){
    std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__ << " xsec_cov is null!!" << std::endl;
  }
  //DB Need DetID to be before assigned before SetupSplineInfoArray
  nutype = ntype;
  BinningOpt = 0;
  SampleDetID = DetID;

  FindUniqueModes();
  xsec_cov->GetSplineParsNamesFromDetID(25);
  std::cout << "Creating Etrue-Erec splines" << std::endl;
#if USE_SPLINE_FD == USE_TSpline3_FD
  std::cout << "Using normal TSpline3 splines" << std::endl;
#elif USE_SPLINE_FD == USE_TSpline3_red
  std::cout << "Using reduced TSpline3 splines!" << std::endl;
#endif

  std::cout << "About to call SetupSplineInfoArray" << std::endl;
  // ETA - we've passed a xsec_cov to splineFDBase then let's use it!
  // Now we want to try to be smarter and use info in the xsec covariance to 
  // automate the spline Setup
  SetupSplineInfoArray(xsec_cov);

  std::cout << "Number of spline params Set to " << number_parms << std::endl;
  std::cout << "Number of events is " << nevents << std::endl;

  TString name_again = TString(name);
  splinefile = new TFile(name_again, "READ");

}


/// ~~~ check if OK
splineFDBase::splineFDBase(const char *name, int ntype, int nevents, double opt_binning, int DetID, covarianceXsec* xsec_cov) // constructor for 2d spline binning
  : splineBase(name, ntype)
{
  if (opt_binning==0.0) std::cout << std::endl << "Error: 2D constructor for splineFDBase called with BinningOpt = 0.0" << std::endl << std::endl;

  BinningOpt = (int)opt_binning;
  if (BinningOpt==2) {
    std::cout << "Creating Etrue-Erec-theta splines" << std::endl;
  }
  else if (BinningOpt==4) {
    std::cout << "Creating Etrue-Cosz-LepMom atmospheric splines" << std::endl;
  }

  SampleDetID = DetID;
  FindUniqueModes();

#if USE_SPLINE_FD == USE_TSpline3_FD
  std::cout << "Using normal TSpline3 splines" << std::endl;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
  std::cout << "Using reduced TSpline3 splines!" << std::endl;
#endif

  // ETA - we've passed a xsec_cov to splineFDBase then let's use it!
  // Now we want to try to be smarter and use info in the xsec covariance to 
  // automate the spline Setup
  SetupSplineInfoArray(xsec_cov);

  std::cout << "Number of spline params Set to " << number_parms << std::endl;
  std::cout << "Setup splinefile to be " << name << std::endl; 
  splinefile = new TFile(name, "READ");
  nutype = ntype;

}

// ---- Destructor ---- //
splineFDBase::~splineFDBase()
{

}

// ---- SetSplineBinning (first: original erec version, then 2d version) ---- //
void splineFDBase::SetSplineBinning() // erec version
{
  // Get binning from first histogram saved in spline file automatically 
  // (x axis = etrue, y axis = erec)
  if(!splinefile){std::cout << "Couldn't find spline file...." << std::endl;}

  TH2F *hist0 = (TH2F*)splinefile->Get("dev_tmp_0_0");
  std::cout << "Looking in splinefile " << splinefile << "for dev_tmp_0_0" << std::endl;
  if (!hist0){
    std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " could not find dev_tmp_0_0 in spline file. Spline binning will not be Set!" << std::endl;
    throw;
  }
  else{
	std::cout << "Found dev_tmp_0_0 in 1D version" << std::endl;
	std::cout << "Setting spline spline binning!!" << std::endl;
  }

  const int netrue = hist0->GetXaxis()->GetNbins(); 
  const double *etruerange = hist0->GetXaxis()->GetXbins()->GetArray();
  enu_spline = new TAxis(netrue, etruerange);

  const int nvar1 = hist0->GetYaxis()->GetNbins();
  const double *var1_range = hist0->GetYaxis()->GetXbins()->GetArray();
  var1_spline = new TAxis(nvar1, var1_range);

}


void splineFDBase::SetSplineBinning(int opt_binning) // enu-var1-var2 version
{
  TH3D *hist0 = (TH3D*)splinefile->Get("dev_tmp_0_0");
  if (!hist0){
    std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " could not find dev_tmp_0_0 in spline file. Spline binning will not be Set!" << std::endl;
	throw;
  }
  else{
	std::cout << "Found dev_tmp_0_0 in 2D version " << std::endl;
	std::cout << "Setting spline spline binning!!" << std::endl;
  }

  const int netrue = hist0->GetXaxis()->GetNbins(); 
  const double *etruerange = hist0->GetXaxis()->GetXbins()->GetArray();
  enu_spline = new TAxis(netrue, etruerange);
 
  const int nvar1 = hist0->GetYaxis()->GetNbins();
  const double *var1_range = hist0->GetYaxis()->GetXbins()->GetArray();

  const int nvar2 =  hist0->GetZaxis()->GetNbins();
  const double *var2_range =  hist0->GetZaxis()->GetXbins()->GetArray();

  if((netrue <= 1 || nvar1 <= 1 || nvar2 <= 1)){
	std::cerr << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " You're Setting up a 2D spline but one of the axes has only one bin or less..." << std::endl;
	std::cerr << "netrue is " << netrue << std::endl;
	std::cerr << "nvar1 is " << nvar1 << std::endl;
	std::cerr << "nvar2 is " << nvar2 << std::endl;
	std::cerr << "I think you've Set up the wrong spline! Maybe you've used 1D splines by mistake?! " << std::endl;
    throw;	
  }

  var1_spline = new TAxis(nvar1, var1_range);
  var2_spline = new TAxis(nvar2, var2_range);

  return;
}

//ETA spline weight dev
// Calc weight for each event using eventSplines array
// saves having to find the splines which apply to each
// event as we've pre-caclulated this!
void splineFDBase::calcWeights(){

  //ETA - switch based on binning opt in this case we just care if it's 1D or 2D
  switch(BinningOpt){
    //1D i.e. Etrue + one other variable
  case 0:
    //Loop over number of params
    for(unsigned int param_i = 0 ; param_i < dev_1D_w.size() ; ++param_i){
      //Loop over number of modes
      for(unsigned int mode_i = 0 ; mode_i < dev_1D_w[0].size() ; ++mode_i){
	//Loop over number of etrue bins
	for(unsigned int etrue_i = 0 ; etrue_i < dev_1D_w[0][0].size() ; ++etrue_i){
	  //Loop over number of var1 bins
	  for(unsigned int var1_i = 0 ; var1_i < dev_1D_w[0][0][0].size() ; ++var1_i){
	    if(dev_1D_vec[param_i][mode_i][etrue_i][var1_i] == NULL){continue;}
	    dev_1D_w[param_i][mode_i][etrue_i][var1_i] = FastSplineEval(dev_1D_vec[param_i][mode_i][etrue_i][var1_i], param_i);
	    if(dev_1D_w[param_i][mode_i][etrue_i][var1_i]<0){dev_1D_w[param_i][mode_i][etrue_i][var1_i]=0;}
	  }
	}
      }
    }
    break;
    
    //2D i.e. Etrue + two other variables
  default: 
    //Loop over number of params
    for(unsigned int param_i = 0 ; param_i < dev_2D_w.size() ; ++param_i){
      //Loop over number of modes
      for(unsigned int mode_i = 0 ; mode_i < dev_2D_w[0].size() ; ++mode_i){
	//Loop over number of etrue bins
	for(unsigned int etrue_i = 0 ; etrue_i < dev_2D_w[0][0].size() ; ++etrue_i){
	  //Loop over number of var1 bins
	  for(unsigned int var1_i = 0 ; var1_i < dev_2D_w[0][0][0].size() ; ++var1_i){
	    //Loop over number of var2 bins
	    for(unsigned int var2_i = 0 ; var2_i < dev_2D_w[0][0][0][0].size() ; ++var2_i){
	      if(dev_2D_vec[param_i][mode_i][etrue_i][var1_i][var2_i] == NULL){continue;}
	      dev_2D_w[param_i][mode_i][etrue_i][var1_i][var2_i] = FastSplineEval(dev_2D_vec[param_i][mode_i][etrue_i][var1_i][var2_i], param_i);
	      if(dev_2D_w[param_i][mode_i][etrue_i][var1_i][var2_i]<0){dev_2D_w[param_i][mode_i][etrue_i][var1_i][var2_i]=0;}
	    }
	  }
	}
      }
    }
    break;
  }

  return;
}

// ---- GetSplineBins (first: original erec version, then var1-var2 version) ---- //

void splineFDBase::GetSplineBins(int &nutype, bool &sig, double &enu, double &var1, unsigned int &enu_bin, unsigned int &var1_bin) // get bins for etrue-erec splines
{

  enu_bin = enu_spline->FindBin(enu)-1;
  var1_bin = var1_spline->FindBin(var1)-1;

  return;
}

void splineFDBase::GetSplineBins(int &nutype, bool &sig, double &enu, double &var1, double &var2, unsigned int &enu_bin, unsigned int &bin1, unsigned int &bin2) // get bins for etrue-var1-var2 splines
{

  enu_bin = enu_spline->FindBin(enu)-1;
  bin1 = var1_spline->FindBin(var1)-1;
  bin2 = var2_spline->FindBin(var2)-1;

  return;
}

std::vector< std::vector<int> > splineFDBase::getEventSplines(int &event_i, int eventmode, unsigned int &enu_bin, unsigned int &var1_bin){
 
  std::vector< std::vector<int> > returnVec;
  int mode = MaCh3Mode_SplineMode_Map[eventmode];

  //DB Now use detid to determine number of spline systematics, names and corresponding modes
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);

  int ebin = enu_bin;
  if(ebin<0 || ebin>=enu_spline->GetNbins()) 
  {
    return returnVec;  
  }
  
  int vbin = var1_bin;
  if(vbin<0 || vbin>=var1_spline->GetNbins())
  {
    return returnVec;
  }

  for(unsigned param_i = 0 ; param_i < (unsigned)numSplineParams ; param_i++){
    std::vector<int> spline_modes = SplineModeVecs[param_i];
    for(int mode_i = 0 ; (unsigned)mode_i < spline_modes.size() ; mode_i++){
      if(mode == spline_modes[mode_i]){
	std::vector<int> vecBin(4);

	vecBin[0] = param_i;
	vecBin[1] = mode;
	vecBin[2] = ebin;
	vecBin[3] = vbin;
	returnVec.push_back(vecBin);
      } 
    }
  }

  return returnVec;
}

//Now the 2D version of the function above
std::vector< std::vector<int> > splineFDBase::getEventSplines(int &event_i, int eventmode, unsigned int &enu_bin, unsigned int &var1_bin, unsigned int &var2_bin){

  //ETA - only support binning option 2 at the moment, should be easy to add in other though - 18/09/2020: DB We have atmospheric xsec splines now :P
  if((BinningOpt != 2)&&(BinningOpt != 4)){std::cerr << "[ERROR:splineFDBase::getEventSplines() - only binning option 2 and 4 are currently implemented here" << std::endl;}

  //DB Now use detid to determine number of spline systematics, names and corresponding modes
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);

  std::vector< std::vector<int> > returnVec;
  int mode = MaCh3Mode_SplineMode_Map[eventmode];

  int ebin = enu_bin;  
  int enuBinLimit = enu_spline->GetNbins();
  int v1BinLimit = -1;
  int v2BinLimit = -1;

  v1BinLimit = var1_spline->GetNbins();
  v2BinLimit = var2_spline->GetNbins();

  if (v1BinLimit==-1 || v2BinLimit==-1) {
    std::cerr << "[ERROR:splineFDBase::getEventSplines() - binning limits not Set correctly" << std::endl;
    throw;
  }

  if(ebin<0 || ebin>=enuBinLimit) {
    return returnVec;  
  }

  int vbin = var1_bin;
  if(vbin<0 || vbin>=v1BinLimit) {
    return returnVec;
  }

  int v2bin = var2_bin;
  if(v2bin<0 || v2bin>=v2BinLimit) {
    return returnVec; 
  }

  for(unsigned param_i = 0 ; param_i < (unsigned)numSplineParams ; param_i++){
    std::vector<int> spline_modes = SplineModeVecs[param_i];
    for(int mode_i = 0 ; (unsigned)mode_i < spline_modes.size() ; mode_i++){
      if(mode == spline_modes[mode_i]) {
	std::vector<int> vecBin(5);
   
	vecBin[0] = param_i;
	vecBin[1] = mode;
	vecBin[2] = ebin;
	vecBin[3] = vbin;
	vecBin[4] = v2bin;
	returnVec.push_back(vecBin);
      }
    } 
  }
  
  return returnVec;
}

//Needed for FastSplineEval
void splineFDBase::SetupSplineInfoArray(covarianceXsec * xsec){

  if(!xsec){std::cout << "[ERROR]:: " << __FILE__ << ":" << __LINE__ << " xsec cov is NULL" << std::endl;}
  else{
    std::cout << " xsec cov is not null " << std::endl;
  }
  covxsec = xsec;

  std::cout << "Now in SetupSplineInfoArray!!" << std::endl;
  std::cout << "xsec is " << xsec << std::endl;

  //DB Now use detid to determine number of spline systematics, names and corresponding modes 
  std::vector<std::string> splinenames = covxsec->GetSplineParsNamesFromDetID(SampleDetID);
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);

  std::cout << "Found " << splinenames.size() << " spline names" << std::endl;
  std::cout << "Found " << numSplineParams << std::endl;

  // Total XSec params might not be the same as far detector params  
  // So we have to select by hand which ones we want    
  int xsecnParams = numSplineParams;
  number_parms = xsecnParams;
  //std::cout << "xsecnParams: " << xsecnParams << std::endl;
  std::vector<std::string> xsecParnames = SplineFileParsNames;
  std::vector<int> xsecParsIndex = SplineParsIndex;

  //ETA - ADD NUMBER OF SPLINES
  unsigned int n_far_splines = numSplineParams;

  //for(unsigned int i = 0; i < xsecParnames.size(); i++){
  //std::cout << " Xsec parname: \t" <<  xsecParnames[i] << ", index: " << xsecParsIndex[i] << std::endl;
  //}

  // ETA - replacing old hard coded array of spline parameter names
  // we can just get these from the xsec covariance.
  std::string FarXsecPars[n_far_splines];
  for(unsigned int splinePar_i = 0 ; splinePar_i < n_far_splines ; splinePar_i++){
	FarXsecPars[splinePar_i] = xsecParnames[splinePar_i];
  }
  unsigned int nPars = (sizeof(FarXsecPars)/sizeof(*FarXsecPars));

  // Check if an unexpected parameter is loaded
  // ETA - only fill splineParsNames if we haven't done this already
  // this assumes you're not trying to change xsec covariance but 
  // I don't think you'd ever want to?
  if(splineParsNames.size() == 0){
	for(unsigned int i = 0; i < xsecParnames.size(); i++){
	  std::string * p = std::find(FarXsecPars, FarXsecPars+nPars, xsecParnames[i]);
	  if(p != FarXsecPars+nPars){
	    splineParsNames.push_back(xsecParnames[i]);
	    splineParsIndex.push_back(xsecParsIndex[i]);
	  }
	}
  }
  //else{std::cout << "splineParsNames and splineParsIndex have already been filled" << std::endl;}

  // Also check that each of the above expected parameters are actually picked up
  // and unique
  for(unsigned int i = 0; i < nPars; i++){
    std::vector<std::string>::iterator it = std::find(splineParsNames.begin(), splineParsNames.end(), FarXsecPars[i]);
    if(it == splineParsNames.end()){
      std::cerr << "Error: Far detector xsec par " << FarXsecPars[i] << " not loaded!" << std::endl;
      throw;
    }
  }
    
  nSplineParams = splineParsNames.size();

  SplineInfoArray = new FastSplineInfo[nSplineParams];
  for (int i = 0; i < nSplineParams; ++i) {
    SplineInfoArray[i].nPts = -999;
    SplineInfoArray[i].xPts = NULL;
    SplineInfoArray[i].CurrSegment = -999;
	//ETA Set flat member of FastSplineInfo
	SplineInfoArray[i].flat = 0;
  }
  
  return;
}


void splineFDBase::SetSplineInfoArrays(){
  // Here we try to find a non-empty (or larger than 3 knot) Spline object to record the knots
  // ETA - now we Set any flat splines or splines we can't find to NULL
  // so here we now need to loop through the spines until we find one that isn't NULL to Set 
  // the SplineInfoArray. If all the splines are flat for a parameter we Set the flat member
  // of FastSplineInfoArray to be 1 (i.e true). This is then checked later when finding the
  // segment so we don't get a seg fault. NULL splines return 1 in FastSplineEval

  //DB Now use detid to determine number of spline systematics, names and corresponding modes
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);


  //ETA - new way of Setting up splines
  bool is_null = true;
  switch (BinningOpt)
    {
	  //for 1D spline (i.e. splines binned in etrue-var1)
	  case 0:
	  case 1:
	for(unsigned spline_i = 0 ; spline_i < (unsigned)numSplineParams ; spline_i++){
	  is_null = true;
	  //std::cout << "Filled SplineInfoArray at " << spline_i << " for spline " << SplineFileParsNames[spline_i] << std::endl;
	  for(unsigned mode_i = 0 ; mode_i < (unsigned)nUniqueModes; mode_i++){
		for(unsigned enu_i = 0 ; enu_i < (unsigned)enu_spline->GetNbins() ; enu_i++){
		  for(unsigned var1_i = 0 ; var1_i < (unsigned)var1_spline->GetNbins() ; var1_i++){
			if(SetSplineInfoArray(dev_1D_vec[spline_i][mode_i][enu_i][var1_i], spline_i)){		  
			  is_null = false; break;
			}
		  }
		  if(!is_null){break;}
		}
		if(!is_null){break;}
	  }	  
	  if(is_null){SplineInfoArray[spline_i].flat = 1;}
	}
      break;
    //for 2D spline (i.e. splines binned in etrue-var1-var2)
    case 2:
      //std::cout << "BinningOpt: " << BinningOpt << std::endl;                               
	  for(unsigned spline_i = 0 ; spline_i < (unsigned)numSplineParams ; spline_i++){
		is_null = true;
		//std::cout << "Filled SplineInfoArray at " << spline_i << " for spline " << SplineFileParsNames[spline_i] << std::endl;
		for(unsigned mode_i = 0 ; mode_i < (unsigned)nUniqueModes; mode_i++){
		  for(unsigned enu_i = 0 ; enu_i < (unsigned)enu_spline->GetNbins() ; enu_i++){
			for(unsigned var1_i = 0 ; var1_i < (unsigned)var1_spline->GetNbins() ; var1_i++){
			  for(unsigned var2_i = 0 ; var2_i < (unsigned)var2_spline->GetNbins() ; var2_i++){
				if(SetSplineInfoArray(dev_2D_vec[spline_i][mode_i][enu_i][var1_i][var2_i], spline_i)){		  
				  is_null = false; break;
				}
			  }
			  if(!is_null){break;}
			}
			if(!is_null){break;}
		  }	  
		  if(!is_null){break;}
		}
		if(is_null){SplineInfoArray[spline_i].flat = 1;}
	  } 
	  break;
    default:
      std::cout << "Something went wrong in splineFDBase::SetSplineInfoArrays()! BinningOpt: " << BinningOpt << std::endl;
    }

}


// *************************
// Only need to do the binary search once per parameter, not once per bin!
void splineFDBase::FindSplineSegment() {

  // Loop over the splines
  for (int i = 0; i < nSplineParams; i++) {

	//If the spline is always flat then move onto next param
	//You don't need this info for flat splines as they're
	//Set to NULL which gets caught in FastSplineEval
	if(SplineInfoArray[i].flat == 1){continue;}

    const int nPoints = SplineInfoArray[i].nPts;
    const double* xArray = SplineInfoArray[i].xPts;

	// && here is just being too cautious, flat should always be 0
	if ( (nPoints == -999 || xArray == NULL) && SplineInfoArray[i].flat == 0) {
      std::cerr << "ERROR" << std::endl;
      std::cerr << "SplineInfoArray[" << i << "] isn't Set yet" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    // Get the variation for this reconfigure for the ith parameter 
    const int GlobalIndex = splineParsIndex[i];
    double xvar_tmp;
    if(covxsec!=NULL){
      xvar_tmp=covxsec->calcReWeight(GlobalIndex);
    }
    else{
      std::cout << "Null Xsec - Quitting" << std::endl;
      throw;
    }
    const double xvar=xvar_tmp;

    // The segment we're interested in (klow in ROOT code)   
    int segment = 0;
    int kHigh = nPoints-1;

    // If the variation is below the lowest saved spline point
    if (xvar <= xArray[0]) {
      segment = 0;
      // If the variation is above the highest saved spline point
    } else if (xvar >= xArray[nPoints-1]) {
      // Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
      segment = kHigh;
      // If the variation is between the maximum and minimum, perform a binary search     
    } else {
      // The top point we've got                 
      int kHalf = 0;
      // While there is still a difference in the points (we haven't yet found the segment)
      // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment  
      while (kHigh - segment > 1) {
        // Increment the half-step
        kHalf = (segment + kHigh)/2;
        // If our variation is above the kHalf, Set the segment to kHalf   
        if (xvar > xArray[kHalf]) {
          segment = kHalf;
          // Else move kHigh down 
        } else {
          kHigh = kHalf;
        }
      } // End the while: we've now done our binary search
    } // End the else: we've now found our point
    if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;
    // Save the segment for the ith parameter
    SplineInfoArray[i].CurrSegment = segment;
#ifdef DEBUG_WEIGHTS
    if (SplineInfoArray[i].xPts[segment] > xvar && segment != 0) {
      std::cerr << "Found a segment which is _ABOVE_ the variation!" << std::endl;
      std::cerr << "IT SHOULD ALWAYS BE BELOW! (except when segment 0)" << std::endl;
      std::cerr << splineParsNames[i] << std::endl;

      std::cerr << "Found segment   = " << segment << std::endl;
      std::cerr << "Doing variation = " << xvar << std::endl;
      std::cerr << "x in spline     = " << SplineInfoArray[i].xPts[segment] << std::endl;
      for (int j = 0; j < SplineInfoArray[j].nPts; ++j) {
	std::cerr << "    " << j << " = " << SplineInfoArray[i].xPts[j] << std::endl;
      }
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
#endif
  }

}

// This does what it says on the tin, if there are any underlying interaction modes which are being grouped together for the splines then we don't need to keep them both
std::vector< std::vector<int> > splineFDBase::StripDuplicatedModes(std::vector< std::vector<int> > InputVector) {
  
  int InputVectorSize = InputVector.size();
  std::vector< std::vector<int> > ReturnVec(InputVectorSize);

  for (int iVec=0;iVec<InputVectorSize;iVec++) {
    std::vector<int> TmpVec;

    for (unsigned int iMode=0;iMode<InputVector[iVec].size();iMode++) {
      int Mode = InputVector[iVec][iMode];
      bool IncludeMode = true;

      for (unsigned int iDuplicatedMode=0;iDuplicatedMode<DuplicatedFDModes.size();iDuplicatedMode++) {
	if (Mode == DuplicatedFDModes[iDuplicatedMode]) {
	  IncludeMode = false;
	}
      }

      if (IncludeMode) {
	TmpVec.push_back(Mode);
      }
    }
    
    ReturnVec[iVec] = TmpVec;
  }

  return ReturnVec;
}

// ***************************************************************************
// Fast spline evaluation
// Inspired by TSpline3::Eval, which we can speed up considerably
// Main reason is that we know that for one parameter (e.g. MAQE) we will have the same number of points, x min, xmax, etc for all MAQE splines, so we can signficantly reduce number of operations
// The curious can find very similar GPU code in splines/gpuSplineUtils.cu and CPU code in spline/SplineMonolith.cpp::Eval
// I've included it here for more transparency: this kind of eval should be possible for binned splines too
template <class T>
double splineFDBase::FastSplineEval(T* spline, const int SplineNumber) {
  // ***************************************************************************  

  // Check if the  coveriance is NULL 
  if (covxsec == NULL) return 1.0;
  // The segment has already been found in FindSplineSegment()
  int segment = SplineInfoArray[SplineNumber].CurrSegment;

  // These are what we can extract from the TSpline3  
  double x = -999.99;
  double y = -999.99;
  double b = -999.99;
  double c = -999.99;
  double d = -999.99;

  // Now write the coefficients  
  spline->GetCoeff(segment, x, y, b, c, d);
  // Get the variation for this reconfigure for the ith parameter   
  int GlobalIndex = splineParsIndex[SplineNumber];
  double xvar;
  xvar=covxsec->calcReWeight(GlobalIndex);

  // The Delta(x)   
  double dx = xvar - x;

  // The spline weight to return           
  double weight = y+dx*(b+dx*(c+d*dx));


  // Check that eval on the TSpline3 is the same as our weight    
#ifdef DEBUG_WEIGHTS

  // Difference between eval and weight   
  double diff =fabs(spline->Eval(xvar) - weight);

  if (diff > 1.E-7) {

    std::cerr << "TSpline3->Eval() != custom eval: Something is wrong with FastSplineEval!" << std::endl;
    std::cerr << "Difference in spline evaluation > 1.E-5!" << std::endl;

    std::cerr << "Eval      = " << spline->Eval(xvar) << std::endl;
    std::cerr << "Cust      = " << weight << std::endl;
    std::cerr << "diff      = " << diff << std::endl;
    std::cerr << "param     = " << splineParsNames[SplineNumber] << std::endl;
    if(covxsec){
      std::cerr << "variation = " << covxsec->calcReWeight(GlobalIndex) << std::endl;
      std::cerr << "paramVal  = " << covxsec->GetParProp(GlobalIndex) << std::endl;
    }

    // Check we've found the right segment    
    int klow = spline->FindX(xvar);
    if (klow >= spline->GetNp()-1 && spline->GetNp() > 1) klow = spline->GetNp()-2;

    std::cerr << "segment   = " << segment << std::endl;
    std::cerr << "spl segm  = " << klow << std::endl;
    std::cerr << "nPoints   = " << SplineInfoArray[SplineNumber].nPts << std::endl;
    std::cerr << "nPoints sp= " << spline->GetNp() << std::endl;


    std::cerr << "Printing x information:" << std::endl;
    for (int i = 0; i < SplineInfoArray[SplineNumber].nPts; ++i) {
      std::cerr << "   " << i << " = " << SplineInfoArray[SplineNumber].xPts[i] << std::endl;
    }
    std::cerr << "From spline: " << std::endl;
    for (int i = 0; i < spline->GetNp(); ++i) {
      double xtmp, ytmp;
      spline->GetKnot(i, xtmp, ytmp);
      std::cerr << "   " << i << " = " << xtmp << std::endl;
    }

    if (klow != segment) {
      std::cerr << "Have found wrong segment in FindSplineSegment!" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
#endif // DEBUG_WEIGHTS
  return weight;
}
//ETA testing this out
template<class T>
bool splineFDBase::SetSplineInfoArray(T* spline, int isyst){
  // Fill the SplineInfoArray entries with information on each splinified parameter                                                                                                                             
  if (spline == NULL){return false;}
  if (SplineInfoArray[isyst].xPts == NULL) {
    // Fill the number of points 
	SplineInfoArray[isyst].nPts = spline->GetNp();
    if(SplineInfoArray[isyst].nPts == 3){
      return false;
    }else{
      // Fill the x points         
      SplineInfoArray[isyst].xPts = new double[SplineInfoArray[isyst].nPts];
      for (int k = 0; k < SplineInfoArray[isyst].nPts; ++k) {
        double xtemp = -999.99;
        double ytemp = -999.99;
        spline->GetKnot(k, xtemp, ytemp);
        SplineInfoArray[isyst].xPts[k] = xtemp;
      }
    }
  }

  return true;
}
