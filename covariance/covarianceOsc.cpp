#include "covarianceOsc.h"

covarianceOsc::covarianceOsc(const char* name, const char *file, TH2D *hist_dcpth13NH, TH2D *hist_dcpth13IH, TH2D *hist_23)
: covarianceBase(name, file) {

  if (hist_dcpth13NH) {
    h_dcpth13NH = hist_dcpth13NH;
    std::cout << "Using delta cp and th13 correlated prior from histogram" << std::endl;
    std::cout << "NOTE: THIS WILL IGNORE ANY CORRELATIONS BETWEEN THESE PARAMETERS AND OTHER PARAMETERS" << std::endl;
  } else {
    h_dcpth13NH = NULL;
  }

  if (hist_dcpth13IH) {
    h_dcpth13IH = hist_dcpth13IH;
    std::cout << "Using delta cp and th13 IH correlated prior from histogram" << std::endl;
    std::cout << "NOTE: THIS WILL IGNORE ANY CORRELATIONS BETWEEN THESE PARAMETERS AND OTHER PARAMETERS" << std::endl;
  } else {
    h_dcpth13IH = NULL;
  }

  if (hist_23) {
    h_23 = hist_23;
    std::cout << "Using th23 and dm23 correlated prior from histogram" << std::endl;
    std::cout << "NOTE: THIS WILL IGNORE ANY CORRELATIONS BETWEEN th_23 PARAMETERS AND OTHER PARAMETERS" << std::endl;
  } else {
    h_23 = NULL;
  }

  //Read in osc pars from xml file
  TFile *infile = new TFile(file, "READ");
  osc_prior = (TVectorD*)infile->Get("osc_nom");
  TVectorD* osc_stepscale = (TVectorD*)infile->Get("osc_stepscale");
  TVectorD* osc_sigma = (TVectorD*)infile->Get("osc_sigma");
  TVectorD* osc_flat_prior = (TVectorD*)infile->Get("osc_flat_prior");

  TObjArray* objarr_name = (TObjArray*)(infile->Get("osc_param_names"));

  TVectorD* osc_baseline = (TVectorD*)infile->Get("osc_baseline");
  TVectorD* osc_density = (TVectorD*)infile->Get("osc_density");
  double fScale = 1.0;

  //KS: Save all neccesary information from covariance
  for(int io = 0; io <size; io++)
  {

	_fNames[io] = std::string(((TObjString*)objarr_name->At(io))->GetString());
    
    _fPreFitValue[io]  = (*osc_prior)(io);
    _fCurrVal[io] = _fPropVal[io] = _fPreFitValue[io];
    _fError[io] = (*osc_sigma)(io);
 
    _fIndivStepScale[io] = fScale * (*osc_stepscale)(io);
    
	std::cout << "Setting step scale to " << _fIndivStepScale[io] << std::endl; 

    //KS: Set flat prior
    if( (bool)((*osc_flat_prior)(io)) ) setEvalLikelihood(io,false);
  }
    
  L = (*osc_baseline)(0);
  density = (*osc_density)(0);

  flipdelM=false;
  reactorPrior = false;
  flipBeta=false;

  fixdm23NH = -999;
  fixdm23IH = -999;
  fixth23NH = -999;
  fixth23IH = -999;

  randomize();

  oscpars1 = new double[10];
  
  //KS:those are constant no need to overwrite them each time
  oscpars1[6] = 2;
  oscpars1[7] = L;
  oscpars1[8] = density;
  
  Print();
  CheckOrderOfParams();
    
  infile->Close();
  delete infile;
  
  std::cout << "created oscillation parameter handler" << std::endl;
}


covarianceOsc::~covarianceOsc()
{
}

double covarianceOsc::GetLikelihood() {
  double logL = 0.0;
  if(size==6) {
    #ifdef MULTITHREAD
    #pragma omp parallel for reduction(+:logL)
    #endif
    for(int i = 0; i < size; i++) 
    {
      for(int j = 0; j <= i; j++)
      {
        // If parameter 23 histogram exists, use that instead of matrix
        if (h_23 && (i==1 || i==4 || j==1 || j==4))
          continue;
        // If th13-dcp histogram exists, use that instead of matrix
        if ((h_dcpth13NH || h_dcpth13IH) && (i==2 || j==2)) // dcp not in loop
          continue;

        if(!_fFlatPrior[i] && !_fFlatPrior[j])
        {
          /*std::cout << i << " " << j << " _fPreFitValue i " << _fPreFitValue[i] << " _fPreFitValue j " << _fPreFitValue[j] << std::endl;
            std::cout << "parcurr " << _fPropVal[i] << " invmatrix " << InvertCovMatrix[i][j] << std::endl;
            std::cout << "diff " << _fPropVal[i] - _fPreFitValue[i] << std::endl;
            std::cout << "likelihood " << 0.5*(_fPropVal[i] - _fPreFitValue[i])*(_fPropVal[j] - _fPreFitValue[j])*InvertCovMatrix[i][j]; << std::endl;*/

          //check
          if (h_23 && (i==1 || i==4 || j==1 || j==4))
            std::cout << "Error: using matrix when parameter-23 histogram exists" << std::endl;
          if ((h_dcpth13NH || h_dcpth13IH) && (i==2 || j==2))
            std::cout << "Error: using matrix when theta13-dcp histogram exists" << std::endl;

          int scale = 1;
          if(i != j) scale = 2;
          logL += scale * 0.5*(_fPropVal[i] - _fPreFitValue[i])*(_fPropVal[j] - _fPreFitValue[j])*InvertCovMatrix[i][j];
        }
      }
    }

    // If parameter 23 histogram exists, use that instead of the matrix
    if (h_23)
    {
      if (!_fFlatPrior[1] && !_fFlatPrior[4])
        logL+=TMath::Log(h_23->Interpolate(_fPropVal[1],_fPropVal[4]))*-1.0;
      else if (!_fFlatPrior[1] || !_fFlatPrior[4])
        std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta23 and dm23 for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
    }

    // If 23-parameters are fixed but you're flipping MH, evaluate the MH prior
    if (_fError[1]<0.0 && _fError[4]<0.0 && flipdelM)
    {
      if (_fPropVal[4]==fixdm23NH)
        logL+=TMath::Log(0.684)*-1.0;
      else if (_fPropVal[4]==fixdm23IH)
        logL+=TMath::Log(0.316)*-1.0;
      else
        std::cerr << "ERROR: dm23 = " << _fPropVal[4] << ", fixdm23NH = " << fixdm23NH << ", fixdm23IH = " << fixdm23IH << std::endl;
    }

    // If theta13-dcp histogram exists, use that instead of the matrix
    if (h_dcpth13IH && _fPropVal[4]<0) // if IH
    {
      if (!_fFlatPrior[2] && !_fFlatPrior[5])
        logL += TMath::Log(h_dcpth13IH->Interpolate(_fPropVal[2],_fPropVal[5]))*-1.0;
      else if (!_fFlatPrior[2] || !_fFlatPrior[5])
        std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta13 and dcp for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
    }      
    else if (h_dcpth13NH) // if NH or IH histogram doesn't exist (default to NH)
    {
      if (!_fFlatPrior[2] && !_fFlatPrior[5])
        logL += TMath::Log(h_dcpth13NH->Interpolate(_fPropVal[2],_fPropVal[5]))*-1.0;
      else if (!_fFlatPrior[2] || !_fFlatPrior[5])
        std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta13 and dcp for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
    }
    else if (!_fFlatPrior[5]) // Evaluate likelihood for dcp (if histogram doesn't exist, th13 prior will already have been evaluated from matrix)
      logL+=1/(2.0*TMath::Pi()); // flat prior for delta
  }

  else if(size==7)
  {
    for(int i = 0; i < size; ++i)
    {
      for(int j = 0; j <= i; ++j)
	  {
        // If parameter 23 histogram exists, use that instead of matrix
        if (h_23 && (i==1 || i==4 || j==1 || j==4))
          continue;
        // If th13-dcp histogram exists, use that instead of matrix
        if ((h_dcpth13NH || h_dcpth13IH) && (i==2 || j==2)) // dcp not in loop
          continue;

        if(!_fFlatPrior[i] && !_fFlatPrior[j])
        {
          /*std::cout << i << " " << j << " _fPreFitValue i " << _fPreFitValue[i] << " _fPreFitValue j " << _fPreFitValue[j] << std::endl;
            std::cout << "parcurr " << fParProp[i] << " invmatrix " << InvertCovMatrix[i][j] << std::endl;
            std::cout << "diff " << fParProp[i] - _fPreFitValue[i] << std::endl;
            std::cout << "likelihood " << 0.5*(_fPropVal[i] - _fPreFitValue[i])*(_fPropVal[j] - _fPreFitValue[j])*InvertCovMatrix[i][j]; << std::endl;*/

          //check
          if (h_23 && (i==1 || i==4 || j==1 || j==4))
            std::cout << "Error: using matrix when parameter-23 histogram exists" << std::endl;
          if ((h_dcpth13NH || h_dcpth13IH) && (i==2 || j==2))
            std::cout << "Error: using matrix when theta13-dcp histogram exists" << std::endl;
          int scale = 1;
          if(i != j) scale = 2;
          logL += scale * 0.5*(_fPropVal[i] - _fPreFitValue[i])*(_fPropVal[j] - _fPreFitValue[j])*InvertCovMatrix[i][j];

        }
      }
    }

    // If parameter 23 histogram exists, use that instead of the matrix
    if (h_23)
    {
      if (!_fFlatPrior[1] && !_fFlatPrior[4])
        logL+=TMath::Log(h_23->Interpolate(_fPropVal[1],_fPropVal[4]))*-1.0;
      else if (!_fFlatPrior[1] || !_fFlatPrior[4])
        std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta23 and dm23 for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
    }

    // If 23-parameters are fixed but you're flipping MH, evaluate the MH prior
	if (_fError[1]<0.0 && _fError[4]<0.0 && flipdelM)
	{
	  if (_fPropVal[4]==fixdm23NH)
		logL+=TMath::Log(0.684)*-1.0;
	  else if (_fPropVal[4]==fixdm23IH)
		logL+=TMath::Log(0.316)*-1.0;
	  else
		std::cerr << "ERROR: dm23 = " << _fPropVal[4] << ", fixdm23NH = " << fixdm23NH << ", fixdm23IH = " << fixdm23IH << std::endl;
	}

	// If theta13-dcp histogram exists, use that instead of the matrix
	// IH
	if (h_dcpth13IH && _fPropVal[4]<0) {
	  if (!_fFlatPrior[2] && !_fFlatPrior[5]) {
		logL += TMath::Log(h_dcpth13IH->Interpolate(_fPropVal[2],_fPropVal[5]))*-1.0;
	  } else if (!_fFlatPrior[2] || !_fFlatPrior[5]) {
		std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta13 and dcp for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
	  } 
	  // if NH or only one histogram given (default to NH)
	} else if (h_dcpth13NH) { 
	  if (!_fFlatPrior[2] && !_fFlatPrior[5]) {
		logL += TMath::Log(h_dcpth13NH->Interpolate(_fPropVal[2],_fPropVal[5]))*-1.0;
	  } else if (!_fFlatPrior[2] || !_fFlatPrior[5]) {
		std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta13 and dcp for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
	  }
	} else if (!_fFlatPrior[5]) {// Evaluate likelihood for dcp (if histogram doesn't exist, th13 prior will already have been evaluated from matrix)
	  logL+=1/(2.0*TMath::Pi()); // flat prior for delta
	}
    
    // Evaluate likelihood for beta
    // Commented out = flat prior 
    //if (fParEvalLikelihood[6])
    // logL+=1/(3.0); 
  }
  //std::cout << "oscpars: " << _fPropVal[0] << "  " << _fPropVal[1] << "  " << _fPropVal[2] << "  " << _fPropVal[3] << "  " << _fPropVal[4] << "  " << _fPropVal[5] <<std::endl;
  //std::cout << "oscllh pre-RC: " << logL << std::endl;

  // reactor prior
  // ETA - we can remove this
  if (reactorPrior)
  {
    // Reactor prior from 2013 PDG: sin^2(2theta13) = 0.095 +/- 0.01
    // Reactor prior from 2014 PDG: sin^2(2theta13) = 0.093 +/- 0.008
    // Reactor prior from 2015 PDG: sin^2(2theta13) = 0.085 +/- 0.005
    // Reactor prior from 2016ca PDG: sin^2(2theta13) = 0.0857 +/- 0.0046
    // Next convert between single angle (what we have) and double angle (what the PDG gives)
    // double dblang = 4*_fPropVal[2]*(1-_fPropVal[2]);
    // double tmp = (dblang-0.0857) / 0.0046;
    
    // Reactor prior from 2018 PDG: sin^2(theta13) = 0.0212 +/- 0.0008
    // Reactor prior from 2019 PDG: sin^2(theta13) = 0.0218 +/- 0.0007
    // Reactor prior from 2021 PDG: sin^2(theta13) = 0.0220 +/- 0.0007
    // This time we don't have to convert between single<->double angle, PDG gives RC in single angle.

    // Now calculate penalty           
    double tmp = (_fPropVal[2]-0.0220) / 0.0007;
    //double tmp = (dblang-0.095) / 0.01;

    // this line for T2K joint fit result, NOT REACTOR CONSTRAINT!
    //double tmp = (_fPropVal[2] - 0.04571) / 0.01125; 

    // Finally: add to likelihood
    logL += 0.5 * tmp * tmp;
    //std::cout << "oscllh post-RC: " << logL << std::endl;
  }
  return logL;
}

double *covarianceOsc::getPropPars()
{
  for(int i = 0; i < 6; i++)
    oscpars1[i] = _fPropVal[i];

  //Those are constant we initalised them already
  //oscpars1[6] = 2;
  //oscpars1[7] = L;
  //oscpars1[8] = density;
  if(size==7)
    oscpars1[9]=_fPropVal[6];
  else
    oscpars1[9]=1;

  return oscpars1;
}

void covarianceOsc::proposeStep() {
  if(use_adaptive && total_steps<upper_adapt) updateAdaptiveCovariance();

  randomize();
  CorrelateSteps();

  //ETA
  //this won't work if abs(_fPropVal) > 2pi so we should consider
  //plonking a while here
  if(_fPropVal[5] > TMath::Pi()) {
	//std::cout << "PROPOSE STEP:_fPropVal[5] was " << _fPropVal[5] << std::endl;
    _fPropVal[5] = (-2.*TMath::Pi() + _fPropVal[5]);
	//std::cout << "Wrapping delta-CP: _fPropVal[5] is now " << _fPropVal[5] << std::endl;
  } else if (_fPropVal[5] < -TMath::Pi()) {
	//std::cout << "PROPOSE STEP:_fPropVal[5] was " << _fPropVal[5] << std::endl;
    _fPropVal[5] = (2.*TMath::Pi() + _fPropVal[5]);
	//std::cout << "Wrapping delta-CP: _fPropVal[5] is now " << _fPropVal[5] << std::endl;
  }
  
  // Okay now we've done the standard steps, we can add in our nice flips
  // hierarchy flip first
  if(random_number[0]->Uniform()<0.5 && flipdelM){
    _fPropVal[4]*=-1;
  }
  // now octant flip
  if(random_number[0]->Uniform()<0.5){
    // flip octant around point of maximal disappearance (0.5112)
    // this ensures we move to a parameter value which has the same oscillation probability
    _fPropVal[1] = 0.5112 - (_fPropVal[1] - 0.5112);
  }
  // And the beta flip (not sure if beta will still work...)
  if(random_number[0]->Uniform()<0.5 && flipBeta){
    _fPropVal[6]=1-_fCurrVal[6];
  }


  // HI This bit lived outside the loop anywhere, let's leave it!
  // if flipdelM and parameter 4 is fixed, flip between fixed parameters for two hierarchies (Note: this will flip parameters 4 and 1 - dm23 *and* theta23).
  if (_fError[4] < 0.0 && flipdelM)
  {
    // First: check that the NH and IH fixed values were set
    if (fixdm23NH == -999 || fixdm23IH == -999 || fixth23NH == -999 || fixth23IH == -999)
      std::cerr << "ERROR: cannot flip between fixed values of normal and inverted heirarchy if values are not assigned. You have provided: fixdm23NH = " << fixdm23NH << ", fixdm23IH = " << fixdm23IH << ", fixth23NH = " << fixth23NH << ", fixth23IH = " << fixth23IH << std::endl;

    if (_fError[1] > 0.0) std::cout << "ERROR: cannot use this code to flip heirarchies unless BOTH dm23 and th23 are fixed" << std::endl;

    // Check _fCurrVal corresponds to one of the values given
    if (!(_fCurrVal[4] == fixdm23NH || _fCurrVal[4] == fixdm23IH))
      std::cout << "ERROR: _fCurrVal[4] = " << _fCurrVal[4] << " is not equal to fixdm23NH (" << fixdm23NH << ") or fixdm23IH (" << fixdm23IH << "). Not changing this parameter." << std::endl;

    if (!(_fCurrVal[1] == fixth23IH || _fCurrVal[1] == fixth23NH))
      std::cout << "ERROR: _fCurrVal[1] = " << _fCurrVal[1] << " is not equal to fixth23NH (" << fixth23NH << ") or fixth23IH (" << fixth23IH << "). Not changing this parameter." << std::endl;

    // Flip parameters together
    double a = random_number[0]->Uniform();
    if (a<0.5) 
    {
      if (_fCurrVal[4] == fixdm23NH)
      {
        _fPropVal[4] = fixdm23IH;
        _fPropVal[1] = fixth23IH;
      }
      else if (_fCurrVal[4] == fixdm23IH) 
      {
        _fPropVal[4] = fixdm23NH;
        _fPropVal[1] = fixth23NH;
      }
    }
    else 
    {
      _fPropVal[4] = _fCurrVal[4];
      _fPropVal[1] = _fCurrVal[1];
    }
    //std::cout << a << "\t" << _fPropVal[4] << "\t" << _fPropVal[1] << std::endl;
  }
}

std::vector<double> covarianceOsc::defaultPars(bool doubled)
{
  std::vector<double> oscpars;
  if(doubled)
  {
    oscpars.push_back(0.857);
    oscpars.push_back(1.0);
    oscpars.push_back(0.098);
    oscpars.push_back(7.5E-5);
    oscpars.push_back(2.5E-3);
    oscpars.push_back(0);
  }
  else
  {
    oscpars.push_back(0.311);
    oscpars.push_back(0.5);
    oscpars.push_back(0.0251);
    oscpars.push_back(7.5e-05);
    oscpars.push_back(0.0024);
    oscpars.push_back(0);
  }
  return oscpars;
}

void covarianceOsc::setExtraBranches(TTree &tree)
{
  // set branches to save current and proposed osc pars for dm_32 and th_23
  // Is this any different for any of the other parameters?
  // Or can this SetExtraBranches just be in the base class?
  for (int i = 0; i<size; ++i)
  {
    if (!(i==1 || i==4)) continue;

    char bit_c[1024] = "c_";
    strcat(bit_c, _fNames[i].c_str());
    strcat(bit_c, "/D");
    char name_c[1024] = "c_";
    strcat(name_c, _fNames[i].c_str());
    tree.Branch(name_c, (double*)&_fCurrVal[i], bit_c);

    char bit_p[1024] = "p_";
    strcat(bit_p, _fNames[i].c_str());
    strcat(bit_p, "/D");
    char name_p[1024] = "p_";
    strcat(name_p, _fNames[i].c_str());
    tree.Branch(name_p, (double*)&_fPropVal[i], bit_p);
  }
}

// Overload setFlipDeltaM23 to provide values for fixed NH and IH dm23 and th23
void covarianceOsc::setFlipDeltaM23(double dm23NH, double dm23IH, double th23NH, double th23IH)
{
  fixdm23NH = dm23NH;
  fixdm23IH = dm23IH;
  fixth23NH = th23NH;
  fixth23IH = th23IH;

  flipdelM = true;
}

//KS: Print all usefull informations after initialization
void covarianceOsc::Print() {
  std::cout << "Number of pars: " << size << std::endl;
  std::cout << "current " << matrixName << " parameters:" << std::endl;
  std::cout << std::left << std::setw(5) << "#" << std::setw(2) << "|" << std::setw(25) << "Name" << std::setw(2) << "|" << std::setw(10) << "Nom." << std::setw(2) << "|" << std::setw(15) << "IndivStepScale" << std::setw(2) << "|" <<std::setw(15) << "_fError"  << std::endl;
  for(int i = 0; i < size; i++) {
    std::cout << std::fixed << std::setprecision(5) << std::left << std::setw(5) << i << std::setw(2) << "|" << std::setw(25) << _fNames[i].c_str() << std::setw(2) << "|" << std::setw(10) << _fPreFitValue[i]<< std::setw(2) << "|" << std::setw(15) << _fIndivStepScale[i] << std::setw(2) << "|" << std::setw(15)<< _fError[i]<< std::endl;
  }
  
  std::cout<<"Baseline: "<<L<<std::endl;
  std::cout<<"Earth Density: "<<density<<std::endl;
}


//KS: Currently prob3++/probgp requiers particular order so we need to check this is the case
void covarianceOsc::CheckOrderOfParams() 
{
    std::vector<int> wrongParam;
    bool wrongMatrix = false;
    if(strcmp( _fNames[0].c_str(),"sin2th_12") != 0 ){wrongParam.push_back(0); wrongMatrix = true;};
    if(strcmp( _fNames[1].c_str(),"sin2th_23") != 0 ){wrongParam.push_back(1); wrongMatrix = true;};
    if(strcmp( _fNames[2].c_str(),"sin2th_13") != 0 ){wrongParam.push_back(2); wrongMatrix = true;};
    if(strcmp( _fNames[3].c_str(),"delm2_12")  != 0 ){wrongParam.push_back(3); wrongMatrix = true;};
    if(strcmp( _fNames[4].c_str(),"delm2_23")  != 0 ){wrongParam.push_back(4); wrongMatrix = true;};
    if(strcmp( _fNames[5].c_str(),"delta_cp")  != 0 ){wrongParam.push_back(5); wrongMatrix = true;};
    if(size == 7 && strcmp( _fNames[6].c_str(),"beta")  != 0 ){wrongParam.push_back(6); wrongMatrix = true;};
    
    if(wrongMatrix)
    {
        for(unsigned int i =0; i < wrongParam.size(); i++ )
        {
            std::cerr << "Osc Patameter "<< _fNames[i].c_str() <<" isn't in good order"<<std::endl;  
        }
        std::cerr << "Currently prob3++/probgp requiers particular order"<< std::endl;
        std::cerr << "Please modify XML and make new matrix with good order"<< std::endl;
        std::cerr << "Find me here "<<__FILE__ << ":" << __LINE__ << std::endl;
        throw;  
    }

}
