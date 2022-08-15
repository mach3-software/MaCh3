#include <iostream>
#include <stdlib.h>
#include "splineFDBase.h"
#include "TF1.h"

//#define DEBUG_WEIGHTS


// ----- Constructor (first: original erec version, then 2d version) ---- //

splineFDBase::splineFDBase(const char *name, int ntype, int nevents, int DetID, covarianceXsec* xsec_cov) // constructor for erec spline binning
  : splineBase(name, ntype)
{
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

  splinefile = new TFile(name, "READ");
  //ETA - this is a bit weird and I think should be Setup in a different way

  std::cout << "About to call SetupSplines() " << std::endl;

  SetupSplines();

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
  
  splinefile = new TFile(name, "READ");
  nutype = ntype;
  SetupSplines(BinningOpt); //~~~

}

// ---- Destructor ---- //

splineFDBase::~splineFDBase()
{

}


// ---- SetupSplines (first: original erec version, then 2d version) ---- //
void splineFDBase::SetupSplines()
{
  // ETA - need to think about how to do this configurably. If we store all the dev_blah_sp in one vector then can loop through giving the address of each?
  // Also need to pass in the name of the spline in the splinefile from the xsec xml to the covarianceXsec class
  std::vector<syst*> systs;

#if USE_SPLINE_FD == USE_TSpline3_red_FD
  std::cout << "###########################" << std::endl;
  std::cout << "USING TSPLINE3 RED !!!!!" << std::endl;
  std::cout << "###########################" << std::endl;
#endif


  // Set spline binning
  SetSplineBinning();

  //vector to keep track of which splines are flat. We use this later on to make
  //sure we've loaded everything correctly 
  std::vector<std::vector<std::vector<std::vector<bool> > > > flat_vec;

  //ETA - testing new Setup
  //DB Now use detid to determine number of spline systematics, names and corresponding modes
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::cout << "Length of SplineFileParsNames is " << SplineFileParsNames.size() << std::endl;
  std::cout << "SampleDetID is " << SampleDetID << std::endl;
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);

  std::cout << "Expecting " << numSplineParams << " splines " << std::endl;

  for(int isyst=0; isyst<numSplineParams ; isyst++){  // loop over systematics 
	//std::cout << "On isyst " << isyst << std::endl;
#if USE_SPLINE_FD == USE_TSpline3_FD
	std::vector<std::vector<std::vector<TSpline3*> > > tmp_tmp_imode;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	std::vector<std::vector<std::vector<TSpline3_red*> > > tmp_tmp_imode;
#endif
	std::vector<std::vector<std::vector<bool> > > tmp_flat_mode; 
	//ETA adding in this to store weights for all splines
	std::vector<std::vector<std::vector<double> > > tmp_w_mode; 
        std::cout << "Num of modes: " << nUniqueModes << std::endl;
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
#if USE_SPLINE_FD == USE_TSpline3_FD
	  std::vector<std::vector<TSpline3*> > tmp_mbin;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	  std::vector<std::vector<TSpline3_red*> > tmp_mbin;
#endif
	  std::vector<std::vector<bool> > tmp_flat_enu;
	  std::vector<std::vector<double> > tmp_w_enu;
	  for(int ienu = 0; ienu < enu_spline->GetNbins(); ienu++){ // loop over true nu energy
#if USE_SPLINE_FD == USE_TSpline3_FD
		std::vector<TSpline3*> tmp_enu;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
		std::vector<TSpline3_red*> tmp_enu;
#endif
		std::vector<bool> tmp_flat_var1;
		std::vector<double> tmp_w_erec;
		for(int ierec = 0; ierec < var1_spline->GetNbins(); ierec++){ // loop over 1st variable
#if USE_SPLINE_FD == USE_TSpline3_FD
		  TSpline3 *tmp_erec=NULL;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
		  TSpline3_red *tmp_erec=NULL;
#endif
		  tmp_enu.push_back(tmp_erec);
		  tmp_flat_var1.push_back(false);//assume everything isn't flat intially
		  tmp_w_erec.push_back(1.0); // All weights can just be Set to 1
		} // end ierec loop
		tmp_mbin.push_back(tmp_enu);
		tmp_flat_enu.push_back(tmp_flat_var1);
		tmp_w_enu.push_back(tmp_w_erec);
	  } // end ienu loop               
	  tmp_tmp_imode.push_back(tmp_mbin);
	  tmp_flat_mode.push_back(tmp_flat_enu);
	  tmp_w_mode.push_back(tmp_w_enu);
	}// end of mode loop
	dev_1D_vec.push_back(tmp_tmp_imode);
	flat_vec.push_back(tmp_flat_mode);
	dev_1D_w.push_back(tmp_w_mode);
  }//end of syst loop

  for(int isyst=0 ; isyst < numSplineParams ; isyst++){
    syst* temp = new syst(SplineFileParsNames[isyst], &(dev_1D_vec.at(isyst)));
    systs.push_back(temp);
  }
 
  // Dummy spline: flat     
  TGraph *dummy_gr = new TGraph();
  dummy_gr->SetPoint(0,-99999999999,1);
  dummy_gr->SetPoint(1,0,1);
  dummy_gr->SetPoint(2,99999999999,1);

  /////////////////
  // Now load the splines from the spline file
  ////////////////

  TIter next(splinefile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TSpline3")) continue;

	char* splinename=(char*)key->GetName();
	//std::cout << "Spline is " << splinename << std::endl;
	char* syst;
	char* mode;
	int etruebin;
	int erecbin;

	char* tok= strtok(splinename,"_");//dev
	tok = strtok (NULL, "_");//syst
	syst = tok;

	int systnum=-1;
	for(unsigned isyst=0; isyst<systs.size(); isyst++){  // loop over systematics
	  if(strcmp(syst,systs.at(isyst)->name.c_str())==0){
		systnum=isyst;
		break;
	  }
	}

	//If the syst doesn't match any of the spline names then skip it
	//e.g. LowQ2suppression splines that we don't use in MaCh3
	if(systnum==-1){
	  continue;
	}	

	int modenum=-1;
	mode = strtok (NULL, "_");//mode
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
	  if(strcmp(mode,(UniqueModeFarSplineNames[imode]).c_str())==0){
		modenum=imode;
		break;
	  }
	}

	if(modenum==-1){
	  std::cout << "COULDN'T MATCH " << syst << std::endl;
	  std::cout << "No matching mode found for this spline... this shouldn't happen " << std::endl;
      throw;
	}

	tok = strtok (NULL, "_");//sp
	etruebin = atoi(strtok (NULL, "_"));//x
	erecbin = atoi(strtok (NULL, "_"));//y

	TSpline3 *h = (TSpline3*)key->ReadObj();
#if USE_SPLINE_FD == USE_TSpline3_FD
	TSpline3 *spl=(TSpline3*)h->Clone();
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	TSpline3_red *spl = new TSpline3_red(h);
#endif
	delete h;

	//loop over all the spline knots and check their value
	//if the value is 1 then Set the flat bool to false
	int n_knots = spl->GetNp();
	bool flat = true;
	for(int knot_i = 0 ; knot_i < n_knots ; knot_i++){
	  double x =-999;
	  double y = -999;
	  spl->GetKnot(knot_i, x, y);
	  if(x == -999 || y == -999){
		std::cerr << "Something has gone wrong... knot position is at -999" << std::endl;
		throw;
	  }
	  double eval = spl->Eval(x);
	  if(eval < 0.99999 || eval > 1.00001){flat = false; break;}
	}

	//If the spline is flat Set it to NULL and update the flat vector
	if(flat){
	  systs.at(systnum)->spline->at(modenum).at(etruebin).at(erecbin)=NULL; 
	  flat_vec.at(systnum).at(modenum).at(etruebin).at(erecbin) = flat;
	}
	else{
	  systs.at(systnum)->spline->at(modenum).at(etruebin).at(erecbin)=spl;
	}		
  }

  for(unsigned isyst=0; isyst<systs.size(); isyst++){  // loop over systematics
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
	  for(int ienu = 0; ienu < enu_spline->GetNbins(); ienu++){ // loop over true nu energy
		for(int ierec = 0; ierec < var1_spline->GetNbins(); ierec++){ // loop over 1st variable
		  bool flat_spline = flat_vec.at(isyst).at(imode).at(ienu).at(ierec); // check is the spline is flat
		  // if the spline is not flat and the spline is NULL then we have a problem!
		  if((systs.at(isyst)->spline->at(imode).at(ienu).at(ierec)==NULL) && (!flat_spline)){
			char sname[50];
			sprintf(sname,"dev_%s_%s_sp_%d_%d",systs.at(isyst)->name.c_str(),UniqueModeFarSplineNames[imode].c_str(),ienu,ierec);
			//Special case for params which apply to all modes i.e. I've Set mode = 12 in xsec cov 
			std::vector<int> modes = SplineModeVecs[isyst]; 
			for(unsigned spline_mode_i = 0 ; spline_mode_i < modes.size() ; spline_mode_i++){
			  if(modes[spline_mode_i] == imode){
				std::cerr << "[ERROR:] splineFDBase::SetupSplines() - cannot FIND Erec SPLINE " << sname << std::endl;
				std::cerr << "[ERROR:] check that the spline name given in the xsec covariance matches that in the spline file" << std::endl;
			  }
			}//End of mode that splines apply to
		  }//End of if
		} // end ierec loop
	  } // end ienu loop               
	}//end of imode loop
  }//end of syst loop

  splinefile->Close();                                       

  //We're now done with these structs so lets delete them
  for(unsigned int syst_i = 0 ; syst_i < systs.size() ; syst_i++){
    delete systs[syst_i];
  }

  return;
}

void splineFDBase::SetupSplines(int opt_binning) // 2d version
{  
  BinningOpt = opt_binning;
  SetSplineBinning(BinningOpt); 

  int Nbins_1st_var=__BAD_SPLINE__, Nbins_2nd_var=__BAD_SPLINE__;

  Nbins_1st_var = var1_spline->GetNbins();
  Nbins_2nd_var = var2_spline->GetNbins();

  if (Nbins_1st_var == __BAD_SPLINE__ || Nbins_2nd_var == __BAD_SPLINE__) {
	std::cout << "Error: Nbins_1st_var or Nbins_2nd_var = __BAD_SPLINE__" << std::endl;
	std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
	exit(-1);
  }

  //DB Now use detid to determine number of spline systematics, names and corresponding modes
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);

  // ETA - need to think about how to do this configurably. If we store all the dev_blah_sp in one vector then can loop through giving the address of each?
  // Also need to pass in the name of the spline in the splinefile from the xsec xml to the covarianceXsec class
  std::vector<syst2D*> systs;

#if USE_SPLINE_FD == USE_TSpline3_red_FD
  std::cout << "###########################" << std::endl;
  std::cout << "USING TSPLINE3 RED !!!!!" << std::endl;
  std::cout << "###########################" << std::endl;
#endif

  std::vector<std::vector<std::vector<std::vector<std::vector<bool> > > > > flat_vec;

  for(int isyst=0; isyst< numSplineParams ; isyst++){  // loop over systematics
#if USE_SPLINE_FD == USE_TSpline3_FD
	std::vector<std::vector<std::vector<std::vector<TSpline3*> > > > tmp_tmp_imode;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	std::vector<std::vector<std::vector<std::vector<TSpline3_red*> > > > tmp_tmp_imode;
#endif
	std::vector<std::vector<std::vector<std::vector<bool> > > > tmp_flat_mode;
	//ETA adding this to store the weights for each spline eval
	std::vector<std::vector<std::vector<std::vector<double> > > > tmp_w_mode;
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
#if USE_SPLINE_FD == USE_TSpline3_FD
	  std::vector< std::vector<std::vector<TSpline3*> > > tmp_mbin;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	  std::vector< std::vector<std::vector<TSpline3_red*> > > tmp_mbin;
#endif
	  std::vector<std::vector<std::vector<bool> > > tmp_flat_enu;
	  std::vector<std::vector<std::vector<double> > > tmp_w_enu;
	  for(int ienu = 0; ienu < enu_spline->GetNbins(); ienu++){ // loop over true nu energy
#if USE_SPLINE_FD == USE_TSpline3_FD
		std::vector<std::vector<TSpline3*> > tmp_enu;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
		std::vector<std::vector<TSpline3_red*> > tmp_enu;
#endif
		std::vector<std::vector<bool> > tmp_flat_var1;
		std::vector<std::vector<double> > tmp_w_var1;
		for(int i1 = 0; i1 < Nbins_1st_var; i1++){ // loop over 1st variable
#if USE_SPLINE_FD == USE_TSpline3_FD
		  std::vector<TSpline3*> tmp_var1;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
		  std::vector<TSpline3_red*> tmp_var1;
#endif
		  std::vector<bool> tmp_flat_var2;
		  std::vector<double> tmp_w_var2;
		  for (int i2 = 0; i2 < Nbins_2nd_var; i2++){ // loop over 2nd variable
#if USE_SPLINE_FD == USE_TSpline3_FD
			TSpline3 *tmp_var2=NULL;
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
			TSpline3_red *tmp_var2=NULL;
#endif
			tmp_var1.push_back(tmp_var2);
			tmp_flat_var2.push_back(false);//assume everything is flat initally
			tmp_w_var2.push_back(1.0); // All weights can just be Set to 1
		  } // end i2 loop
		  tmp_enu.push_back(tmp_var1);
		  tmp_flat_var1.push_back(tmp_flat_var2);
		  tmp_w_var1.push_back(tmp_w_var2);
		} // end i1 loop
		tmp_mbin.push_back(tmp_enu);
		tmp_flat_enu.push_back(tmp_flat_var1);	
		tmp_w_enu.push_back(tmp_w_var1);
	  } // end ienu loop               
	  tmp_tmp_imode.push_back(tmp_mbin);
	  tmp_flat_mode.push_back(tmp_flat_enu);
	  tmp_w_mode.push_back(tmp_w_enu);
	}// end of mode loop
	flat_vec.push_back(tmp_flat_mode);
	dev_2D_vec.push_back(tmp_tmp_imode);
	dev_2D_w.push_back(tmp_w_mode);
  }// end of syst loop

  for(int isyst=0 ; isyst < numSplineParams ; isyst++){
	syst2D* temp = new syst2D(SplineFileParsNames[isyst], &(dev_2D_vec.at(isyst)));
	systs.push_back(temp);
  }

  // Dummy splines: flat               
  TGraph *dummy_gr = new TGraph();
  dummy_gr->SetPoint(0,-99999999999,1);
  dummy_gr->SetPoint(1,0,1);
  dummy_gr->SetPoint(2,99999999999,1);

  /////////////////
  // Now load the splines from the spline file
  ////////////////
    
  //get all spline objects from file
  TIter next(splinefile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TSpline3")) continue;

    //std::cout<< "Spline is " << key->GetName()<<std::endl; 

    char* splinename=(char*)key->GetName();
    char* syst;
    char* tok= strtok(splinename,"_");//dev
    tok = strtok (NULL, "_");//syst
	syst=tok;
    

    int systnum=-1;
    for(unsigned isyst=0; isyst<systs.size(); isyst++){  // loop over systematics
      if(strcmp(syst,systs.at(isyst)->name.c_str())==0){
        systnum=isyst;
	break;
      }
    }
    
    //If the syst doesn't match any of the spline names then skip it
    //e.g. LowQ2suppression splines that we don't use in MaCh3
    if(systnum==-1){
      continue;
    }
    
    int modenum=-1;
    char* mode = strtok (NULL, "_");//mode
    for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
      if(strcmp(mode,(UniqueModeFarSplineNames[imode]).c_str())==0){
	modenum=imode;
	break;
      }
    }
    
    if(modenum==-1){
      std::cout << "COULDN'T MATCH " << syst << std::endl;
      std::cout << "No matching mode found for this spline... this shouldn't happen " << std::endl;
      throw;
    }
    
    tok = strtok (NULL, "_");//sp
    int etruebin = atoi(strtok (NULL, "_"));//x
    int var1bin = atoi(strtok (NULL, "_"));//y
    int var2bin = atoi(strtok (NULL, "_"));//z
    
    TSpline3 *h = (TSpline3*)key->ReadObj();
#if USE_SPLINE_FD == USE_TSpline3_FD
    TSpline3 *spl=(TSpline3*)h->Clone();
#elif USE_SPLINE_FD == USE_TSpline3_red_FD
	TSpline3_red *spl = new TSpline3_red(h);
#endif
	delete h;
	//std::cout << "address is " << spl << std::endl;
	
	//loop over all the spline knots and check their value
	//if the value is 1 then Set the flat bool to false
	int n_knots = spl->GetNp();
	bool flat = true;
	for(int knot_i = 0 ; knot_i < n_knots ; knot_i++){
	  double x =-999;
	  double y = -999;
	  spl->GetKnot(knot_i, x, y);
	  if(x == -999 || y == -999){
		std::cerr << "Something has gone wrong... knot position is at -999" << std::endl;
		throw;
	  }
	  double eval = spl->Eval(x);
	  if(eval < 0.99999 || eval > 1.00001){flat = false; break;}
	}

	//If the spline is flat Set it to NULL and update the flat vector
	if(flat){
	  systs.at(systnum)->spline->at(modenum).at(etruebin).at(var1bin).at(var2bin)=NULL;
	  flat_vec.at(systnum).at(modenum).at(etruebin).at(var1bin).at(var2bin) = true;
	}
	else{
	  systs.at(systnum)->spline->at(modenum).at(etruebin).at(var1bin).at(var2bin)=spl;
	} 

  } 

  for(unsigned isyst=0; isyst<systs.size(); isyst++){  // loop over systematics
	for(int imode = 0; imode<nUniqueModes; imode++){ // loop over modes
	  for(int ienu = 0; ienu < enu_spline->GetNbins(); ienu++){ // loop over true nu energy
		for(int i1 = 0; i1 < Nbins_1st_var; i1++){ // loop over 1st variable
		  for (int i2 = 0; i2 < Nbins_2nd_var; i2++){ // loop over 2nd variable
			bool flat_spline = flat_vec.at(isyst).at(imode).at(ienu).at(i1).at(i2);//check to see if the spline is flat
		  // if the spline is not flat and the spline is NULL then we have a problem!
			if((systs.at(isyst)->spline->at(imode).at(ienu).at(i1).at(i2)==NULL) && (!flat_spline)){
			  char sname[50];
			  sprintf(sname,"dev_%s_%s_sp_%d_%d_%d",systs.at(isyst)->name.c_str(),UniqueModeFarSplineNames[imode].c_str(),ienu,i1,i2);
			  //Special case for params which apply to all modes i.e. I've Set mode = 12 in xsec cov 
			  std::vector<int> modes = SplineModeVecs[isyst];
			  for(unsigned spline_mode_i = 0 ; spline_mode_i < modes.size() ; spline_mode_i++){
				if(modes[spline_mode_i] == imode){
				  std::cerr << "[ERROR:] splineFDBase::SetupSplines() - cannot FIND Erec SPLINE " << sname << std::endl;
				  std::cerr << "[ERROR:] check that the spline name given in the xsec covariance matches that in the spline file" << std::endl;
				}
			  }//End of mode that splines apply to
			}
		  } // end i2 loop
		} // end i1 loop
	  } // end ienu loop               
	}
  }

  splinefile->Close();      

  //We're now done with these structs so lets delete them
  for(unsigned int syst_i = 0 ; syst_i < systs.size() ; syst_i++){
    delete systs[syst_i];
  }

  return;
}


// ---- SetSplineBinning (first: original erec version, then 2d version) ---- //
void splineFDBase::SetSplineBinning() // erec version
{
  // Get binning from first histogram saved in spline file automatically 
  // (x axis = etrue, y axis = erec)
  TH2D *hist0 = (TH2D*)splinefile->Get("dev_tmp_0_0");
  if (!hist0){
    std::cout << "Error: could not find dev_tmp_0_0 in spline file. Spline binning will not be Set!" << std::endl;
    throw;
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
  if (!hist0)
    std::cout << "Error: could not find dev_tmp_0_0 in spline file. Spline binning will not be Set!" << std::endl;

  const int netrue = hist0->GetXaxis()->GetNbins(); 
  const double *etruerange = hist0->GetXaxis()->GetXbins()->GetArray();
  enu_spline = new TAxis(netrue, etruerange);
 
  const int nvar1 = hist0->GetYaxis()->GetNbins();
  const double *var1_range = hist0->GetYaxis()->GetXbins()->GetArray();

  const int nvar2 =  hist0->GetZaxis()->GetNbins();
  const double *var2_range =  hist0->GetZaxis()->GetXbins()->GetArray();

  if((netrue <= 1 || nvar1 <= 1 || nvar2 <= 1)){
	std::cerr << "[ERROR] - You're Setting up a 2D spline but one of the axes has only one bin or less..." << std::endl;
	std::cerr << "I think you've Set up the wrong spline! Maybe you've used 1D splines by mistake?! " << std::endl;
    throw;	
  }

  var1_spline = new TAxis(nvar1, var1_range);
  var2_spline = new TAxis(nvar2, var2_range);

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
}

void splineFDBase::GetSplineBins(int &nutype, bool &sig, double &enu, double &var1, double &var2, unsigned int &enu_bin, unsigned int &bin1, unsigned int &bin2) // get bins for etrue-var1-var2 splines
{
  enu_bin = enu_spline->FindBin(enu)-1;
  bin1 = var1_spline->FindBin(var1)-1;
  bin2 = var2_spline->FindBin(var2)-1;

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
  covxsec = xsec;

  //DB Now use detid to determine number of spline systematics, names and corresponding modes
  std::vector<std::string> splinenames = covxsec->GetSplineParsNamesFromDetID(SampleDetID);
  int numSplineParams = covxsec->GetNumSplineParamsFromDetID(SampleDetID);
  std::vector<std::string> SplineFileParsNames = covxsec->GetSplineFileParsNamesFromDetID(SampleDetID);
  std::vector< std::vector<int> > SplineModeVecs = StripDuplicatedModes(covxsec->GetSplineModeVecFromDetID(SampleDetID));
  std::vector<int> SplineParsIndex = covxsec->GetSplineParsIndexFromDetID(SampleDetID);

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

//TODO (ETA) - need to pass number of interaction modes and unique spline modes to spline object
//the mode inofrmation etc. will be defined for each experiment

void splineFDBase::FindUniqueModes() {

/*  
  nUniqueModes = 0;

  for (int iMode=0;iMode<kMaCh3_nModes;iMode++) {
    if (MaCh3Mode_to_SplineMode(iMode)==iMode) {
      nUniqueModes += 1;
    }
  }

  int Counter = 0;
  UniqueModeFarSplineNames.resize(nUniqueModes);

  for (int iMode=0;iMode<kMaCh3_nModes;iMode++) {
    if (MaCh3Mode_to_SplineMode(iMode)==iMode) {
      UniqueModeFarSplineNames[Counter] = MaCh3mode_ToString((MaCh3_Mode)iMode);
    } else {
      DuplicatedFDModes.push_back(iMode);
    }
    Counter += 1;
  }

  MaCh3Mode_SplineMode_Map.resize(kMaCh3_nModes);
  for (int iMode=0;iMode<kMaCh3_nModes;iMode++) {
    MaCh3Mode_SplineMode_Map[iMode] = MaCh3Mode_to_SplineMode(iMode);
  }
   
*/
 for (int i = 0; i < 1;i++) {
   MaCh3Mode_SplineMode_Map.push_back(i); }
   nUniqueModes = 1;
   UniqueModeFarSplineNames.push_back("ccqe"); 
   std::cout << "Temp implementation" << std::endl;

}


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
