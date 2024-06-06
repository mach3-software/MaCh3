#include "covarianceOsc.h"

// *************************************
covarianceOsc::covarianceOsc(const char* name, const char *file)
: covarianceBase(name, file) {
// *************************************

  //Read in osc pars from xml file
  TFile *infile = new TFile(file, "READ");
  TVectorD* osc_prior = (TVectorD*)infile->Get("osc_nom");
  TVectorD* osc_stepscale = (TVectorD*)infile->Get("osc_stepscale");
  TVectorD* osc_sigma = (TVectorD*)infile->Get("osc_sigma");
  TVectorD* osc_flat_prior = (TVectorD*)infile->Get("osc_flat_prior");

  TObjArray* objarr_name = (TObjArray*)(infile->Get("osc_param_names"));

  TVectorD* osc_baseline = (TVectorD*)infile->Get("osc_baseline");
  TVectorD* osc_density = (TVectorD*)infile->Get("osc_density");
  double fScale = 1.0;

  //KS: Save all necessary information from covariance
  for(int io = 0; io < _fNumPar; io++)
  {
    _fNames[io] = std::string(((TObjString*)objarr_name->At(io))->GetString());
    
    _fPreFitValue[io]  = (*osc_prior)(io);
    _fCurrVal[io] = _fPropVal[io] = _fPreFitValue[io];
    _fError[io] = (*osc_sigma)(io);
 
    _fIndivStepScale[io] = fScale * (*osc_stepscale)(io);
    
    //KS: Set flat prior
    if( (bool)((*osc_flat_prior)(io)) ) setEvalLikelihood(io, false);
  }

  kDeltaCP = -999;
  kDeltaM23 = -999;
  kSinTheta23 = -999;
  kBeta = -999;
  PerformBetaStudy = false;
  for(int io = 0; io < _fNumPar; io++)
  {
    if(_fNames[io] == "delta_cp")  kDeltaCP = io;
    if(_fNames[io] == "delm2_23")  kDeltaM23 = io;
    if(_fNames[io] == "sin2th_23") kSinTheta23 = io;
    if(_fNames[io] == "beta")
    {
      PerformBetaStudy = true;
      kBeta = io;
    }
  }

  L = (*osc_baseline)(0);
  density = (*osc_density)(0);

  flipdelM = false;
  reactorPrior = false;
  flipBeta = false;

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
  
  MACH3LOG_INFO("created oscillation parameter handler");
}

// *************************************
covarianceOsc::~covarianceOsc() {
// *************************************

  delete[] oscpars1;
}

// *************************************
int covarianceOsc::CheckBounds() {
// *************************************
  int NOutside = 0;

  // ensure osc params don't go unphysical
  if (_fPropVal[0] > 1.0 || _fPropVal[0] < 0 ||
      _fPropVal[kSinTheta23] > 1.0 || _fPropVal[kSinTheta23] < 0 ||
      _fPropVal[2] > 1.0 || _fPropVal[2] < 0
      //|| _fPropVal[kDeltaM23] < 0.0 || _fPropVal[kDeltaM23] > 20E-3 // don't let dm32 go to IH
      //|| fabs(_fPropVal[kDeltaM23]) > 0.004 || fabs(_fPropVal[kDeltaM23]) < 0.001
      //|| _fPropVal[kDeltaCP] < -1*TMath::Pi() || _fPropVal[kDeltaCP] > TMath::Pi()
     ){ NOutside++;}

  if(PerformBetaStudy) {
    if(_fPropVal[kBeta] < 0) { // Don't let beta be less than 0 (no upper limit)
      NOutside++;
    }
  }

  return NOutside;
}
// *************************************
double covarianceOsc::GetLikelihood() {
// *************************************

  double logL = 0.0;
  logL += covarianceBase::GetLikelihood();

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

// *************************************
double *covarianceOsc::getPropPars() {
// *************************************

  for(int i = 0; i < 6; i++)
    oscpars1[i] = _fPropVal[i];

  //Those are constant we initialised them already
  //oscpars1[6] = 2;
  //oscpars1[7] = L;
  //oscpars1[8] = density;
  if(PerformBetaStudy)
    oscpars1[9] = _fPropVal[kBeta];
  else
    oscpars1[9] = 1;

  return oscpars1;
}
// *************************************
void covarianceOsc::proposeStep() {
// *************************************

  covarianceBase::proposeStep();

  //ETA
  //this won't work if abs(_fPropVal) > 2pi so we should consider
  //plonking a while here
  if(_fPropVal[kDeltaCP] > TMath::Pi()) {
    _fPropVal[kDeltaCP] = (-2.*TMath::Pi() + _fPropVal[kDeltaCP]);
  } else if (_fPropVal[kDeltaCP] < -TMath::Pi()) {
    _fPropVal[kDeltaCP] = (2.*TMath::Pi() + _fPropVal[kDeltaCP]);
  }
  
  // Okay now we've done the standard steps, we can add in our nice flips
  // hierarchy flip first
  if(random_number[0]->Uniform() < 0.5 && flipdelM){
    _fPropVal[kDeltaM23] *= -1;
  }
  // now octant flip
  if(random_number[0]->Uniform() < 0.5){
    // flip octant around point of maximal disappearance (0.5112)
    // this ensures we move to a parameter value which has the same oscillation probability
    _fPropVal[kSinTheta23] = 0.5112 - (_fPropVal[kSinTheta23] - 0.5112);
  }
  // And the beta flip (not sure if beta will still work...)
  if(flipBeta)
  {
    if(random_number[0]->Uniform() < 0.5) _fPropVal[kBeta] = 1-_fCurrVal[kBeta];
  }
}

// *************************************
void covarianceOsc::setExtraBranches(TTree &tree) {
// *************************************

  // set branches to save current and proposed osc pars for dm_32 and th_23
  // Is this any different for any of the other parameters?
  // Or can this SetExtraBranches just be in the base class?
  for (int i = 0; i < _fNumPar; ++i)
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
// *************************************
//KS: Print all usefull informations after initialization
void covarianceOsc::Print() {
// *************************************

  MACH3LOG_INFO("Number of pars: {}", _fNumPar);
  MACH3LOG_INFO("Current: {} parameters:", matrixName);
  std::cout << std::left << std::setw(5) << "#" << std::setw(2) << "|" << std::setw(25) << "Name" << std::setw(2) << "|" << std::setw(10) << "Nom." << std::setw(2) << "|" << std::setw(15) << "IndivStepScale" << std::setw(2) << "|" <<std::setw(15) << "_fError"  << std::endl;
  for(int i = 0; i < _fNumPar; i++) {
    std::cout << std::fixed << std::setprecision(5) << std::left << std::setw(5) << i << std::setw(2) << "|" << std::setw(25) << _fNames[i].c_str() << std::setw(2) << "|" << std::setw(10) << _fPreFitValue[i]<< std::setw(2) << "|" << std::setw(15) << _fIndivStepScale[i] << std::setw(2) << "|" << std::setw(15)<< _fError[i]<< std::endl;
  }
  
  MACH3LOG_INFO("Baseline: {}", L);
  MACH3LOG_INFO("Earth Density: {}", density);
}

// *************************************
//KS: Currently prob3++/probgpu requires particular order so we need to check this is the case
void covarianceOsc::CheckOrderOfParams()  {
// *************************************

  std::vector<int> wrongParam;
  bool wrongMatrix = false;
  if(strcmp( _fNames[0].c_str(),"sin2th_12") != 0 ){wrongParam.push_back(0); wrongMatrix = true;};
  if(strcmp( _fNames[1].c_str(),"sin2th_23") != 0 ){wrongParam.push_back(1); wrongMatrix = true;};
  if(strcmp( _fNames[2].c_str(),"sin2th_13") != 0 ){wrongParam.push_back(2); wrongMatrix = true;};
  if(strcmp( _fNames[3].c_str(),"delm2_12")  != 0 ){wrongParam.push_back(3); wrongMatrix = true;};
  if(strcmp( _fNames[4].c_str(),"delm2_23")  != 0 ){wrongParam.push_back(4); wrongMatrix = true;};
  if(strcmp( _fNames[5].c_str(),"delta_cp")  != 0 ){wrongParam.push_back(5); wrongMatrix = true;};
  if(PerformBetaStudy && strcmp( _fNames[6].c_str(),"beta")  != 0 ){wrongParam.push_back(6); wrongMatrix = true;};

  if(wrongMatrix)
  {
    for(unsigned int i = 0; i < wrongParam.size(); i++ )
    {
      MACH3LOG_ERROR("Osc Parameter  {} isn't in good order", _fNames[i].c_str());
    }
    MACH3LOG_ERROR("Currently prob3++/probgp requires particular order");
    MACH3LOG_ERROR("Please modify XML and make new matrix with good order");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
}
