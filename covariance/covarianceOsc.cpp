#include "covarianceOsc.h"

// *************************************
covarianceOsc::covarianceOsc(const std::vector<std::string>& YAMLFile, const char *name, double threshold, int FirstPCA, int LastPCA)
: covarianceBase(YAMLFile, name, threshold, FirstPCA, LastPCA){
// *************************************

  kDeltaCP = -999;
  kDeltaM23 = -999;
  kSinTheta23 = -999;
  kSinTheta23 = -999;
  kSinTheta23 = -999;
  for(int io = 0; io < _fNumPar; io++)
  {
    _fNames[io] = _fFancyNames[io];

    if(_fNames[io] == "delta_cp")  kDeltaCP = io;
    if(_fNames[io] == "delm2_23")  kDeltaM23 = io;
    if(_fNames[io] == "sin2th_23") kSinTheta23 = io;
    if(_fNames[io] == "baseline") kBaseline = io;
    if(_fNames[io] == "density") kDensity = io;

    // Set covarianceBase parameters (Curr = current, Prop = proposed, Sigma = step)
    _fCurrVal[io] = _fPreFitValue[io];
    _fPropVal[io] = _fCurrVal[io];
  }

  /// WARNING HARDCODED
  MACH3LOG_CRITICAL("Fixing baseline and density, it is hardcoded sry");
  toggleFixParameter("baseline");
  toggleFixParameter("density");

  /// @todo KS: Technically if we would like to use PCA we have ot initialise parts here...
  flipdelM = false;

  randomize();

  Print();

  MACH3LOG_INFO("Created oscillation parameter handler");
}

// *************************************
covarianceOsc::~covarianceOsc() {
// *************************************

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
  return NOutside;
}

// *************************************
void covarianceOsc::proposeStep() {
// *************************************

  covarianceBase::proposeStep();

  // HW :: This method is a tad hacky but modular arithmetic gives me a headache.
  //        It should now automatically set dcp to be with [-pi, pi]
  if(_fPropVal[kDeltaCP] > TMath::Pi()) {
    _fPropVal[kDeltaCP] = -1*TMath::Pi() + std::fmod(_fPropVal[kDeltaCP], TMath::Pi());
  } else if (_fPropVal[kDeltaCP] < -TMath::Pi()) {
    _fPropVal[kDeltaCP] = TMath::Pi() + std::fmod(_fPropVal[kDeltaCP], TMath::Pi());
  }
  
  // Okay now we've done the standard steps, we can add in our nice flips
  // hierarchy flip first
  if(random_number[0]->Uniform() < 0.5 && flipdelM){
    _fPropVal[kDeltaM23] *= -1;
  }
  // now octant flip
  if(random_number[0]->Uniform() < 0.5) {
    // flip octant around point of maximal disappearance (0.5112)
    // this ensures we move to a parameter value which has the same oscillation probability
    _fPropVal[kSinTheta23] = 0.5112 - (_fPropVal[kSinTheta23] - 0.5112);
  }
}

// *************************************
//KS: Print all useful information's after initialization
void covarianceOsc::Print() {
// *************************************

  MACH3LOG_INFO("Number of pars: {}", _fNumPar);
  MACH3LOG_INFO("Current: {} parameters:", matrixName);

  MACH3LOG_INFO("{:<5} | {:<25} | {:<10} | {:<15} | {:<15} | {:<10}",
                "#", "Name", "Nom.", "IndivStepScale", "_fError", "FlatPrior");

  for(int i = 0; i < _fNumPar; i++) {
    MACH3LOG_INFO("{:<5} | {:<25} | {:<10.4f} | {:<15.2f} | {:<15.4f} | {:<10}",
                  i, _fNames[i].c_str(), _fPreFitValue[i], _fIndivStepScale[i], _fError[i], _fFlatPrior[i]);
  }
}
