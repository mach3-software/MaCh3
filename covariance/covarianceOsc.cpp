#include "covariance/covarianceOsc.h"
#include "unordered_map"

// *************************************
covarianceOsc::covarianceOsc(const std::vector<std::string>& YAMLFile, std::string name, double threshold, int FirstPCA, int LastPCA)
: covarianceBase(YAMLFile, name, threshold, FirstPCA, LastPCA){
// *************************************
  kDeltaCP = -999;
  kDeltaM23 = -999;
  kDeltaM21 = -999;
  kSinTheta23 = -999;
  kSinTheta12 = -999;
  kSinTheta13 = -999;

  for(int io = 0; io < _fNumPar; io++)
  {
    _fNames[io] = _fFancyNames[io];
    if(_fNames[io] == "delta_cp")  kDeltaCP = io;
    if(_fNames[io] == "delm2_23")  kDeltaM23 = io;
    if(_fNames[io] == "sin2th_23") kSinTheta23 = io;
    if(_fNames[io] == "delm2_12")  kDeltaM21 = io;
    if(_fNames[io] == "sin2th_12") kSinTheta12 = io;
    if(_fNames[io] == "sin2th_13") kSinTheta13 = io;

    // Set covarianceBase parameters (Curr = current, Prop = proposed, Sigma = step)
    _fCurrVal[io] = _fPreFitValue[io];
    _fPropVal[io] = _fCurrVal[io];
  }
  // KS: Check if params have been found. Otherwise things can go wrong
  auto CheckInitialisation = [](const std::string& paramName, int paramValue) {
    if (paramValue == -999) {
      MACH3LOG_WARN("Parameter {} has not been initialized properly.", paramName);
      MACH3LOG_WARN("Things are unpredictable");
    }
  };
  CheckInitialisation("delta_cp", kDeltaCP);
  CheckInitialisation("delm2_23", kDeltaM23);
  CheckInitialisation("sin2th_23", kSinTheta23);
  CheckInitialisation("delm2_21", kDeltaM21);
  CheckInitialisation("sin2th_12", kSinTheta12);
  CheckInitialisation("sin2th_13", kSinTheta13);

  /// @todo KS: Technically if we would like to use PCA we have to initialise parts here...
  flipdelM = true;

  randomize();
  Print();
}

// *************************************
covarianceOsc::~covarianceOsc() {
// *************************************

}

// *************************************
void covarianceOsc::proposeStep() {
// *************************************
  covarianceBase::proposeStep();

  // HW It should now automatically set dcp to be with [-pi, pi]
  CircularPrior(kDeltaCP, -TMath::Pi(), TMath::Pi());

  // Okay now we've done the standard steps, we can add in our nice flips

  // hierarchy flip first
  if(random_number[0]->Uniform() < 0.5 && flipdelM){
    double sin2th12 = TMath::Sin(2 * TMath::ASin(TMath::Sqrt(_fPropVal[kSinTheta12])));
    double tanth23 = TMath::Tan(TMath::ASin(TMath::Sqrt(_fPropVal[kSinTheta23])));
    double sinth13 = TMath::Sqrt(_fPropVal[kSinTheta13]);
    double cosdcp = TMath::Cos(_fPropVal[kDeltaCP]);
    double delm2_21 = _fPropVal[kDeltaM21];
    double correction = -delm2_21 * (1 - ((1-_fPropVal[kSinTheta12]) - cosdcp * sinth13 * tanth23 * sin2th12)); // why did I have to add this minus sign???

    _fPropVal[kDeltaM23] = -_fPropVal[kDeltaM23] + 2*correction;

    //_fPropVal[kDeltaM23] *= -1;
    // _fPropVal[kDeltaM23] = -0.0000325 - (_fPropVal[kDeltaM23] + 0.0000325); // flip around the mean value of NO vs IO
    
  }
  // now octant flip
  if(random_number[0]->Uniform() < 0.5) {
    // flip octant around point of maximal disappearance (0.5112)
    // this ensures we move to a parameter value which has the same oscillation probability
    _fPropVal[kSinTheta23] = 0.5112 - (_fPropVal[kSinTheta23] - 0.5112);
  }

  if(random_number[0]->Uniform() < 0.5 && flipWindow){
    // flip delta_cp around mean value
    _fPropVal[kDeltaCP] = multicanonicalSeparateMean - (_fPropVal[kDeltaCP] - multicanonicalSeparateMean);
  }


}

// *************************************
//HW: This method is a tad hacky but modular arithmetic gives me a headache.
void covarianceOsc::CircularPrior(const int index, const double LowBound, const double UpBound) {
// *************************************
  if(_fPropVal[index] > UpBound) {
    _fPropVal[index] = LowBound + std::fmod(_fPropVal[index] - UpBound, UpBound - LowBound);
  } else if (_fPropVal[index] < LowBound) {
    _fPropVal[index] = UpBound - std::fmod(LowBound - _fPropVal[index], UpBound - LowBound);
  }
}

void covarianceOsc::initialiseAdaption(const YAML::Node& adapt_manager){
    covarianceBase::initialiseAdaption(adapt_manager);
    std::unordered_map<int, double> flip_indices;
    flip_indices[kDeltaM23] = 0.0;
    flip_indices[kSinTheta23] = 0.5112;
    AdaptiveHandler.AddFlipParameters(flip_indices);
}

// *************************************
//KS: Print all useful information's after initialization
void covarianceOsc::Print() {
// *************************************
  MACH3LOG_INFO("Number of pars: {}", _fNumPar);
  MACH3LOG_INFO("Current: {} parameters:", matrixName);

  MACH3LOG_INFO("=================================================================================================================================");
  MACH3LOG_INFO("{:<5} {:2} {:<25} {:2} {:<10} {:2} {:<15} {:2} {:<15} {:2} {:<10} {:2} {:<10}",
                "#", "|", "Name", "|", "Prior", "|", "IndivStepScale", "|", "Error", "|", "FlatPrior", "|", "SampleNames");
  MACH3LOG_INFO("---------------------------------------------------------------------------------------------------------------------------------");
  for (int i = 0; i < _fNumPar; i++) {
    std::string SampleNameString = "";
    for (const auto& SampleName : _fSampleNames[i]) {
      if (!SampleNameString.empty()) {
        SampleNameString += ", ";
      }
      SampleNameString += SampleName;
    }

    MACH3LOG_INFO("{:<5} {:2} {:<25} {:2} {:<10.4f} {:2} {:<15.2f} {:2} {:<15.4f} {:2} {:<10} {:2} {:<10}",
                  i, "|", _fNames[i].c_str(), "|", _fPreFitValue[i], "|", _fIndivStepScale[i], "|", _fError[i], "|", _fFlatPrior[i], "|", SampleNameString);
  }
  MACH3LOG_INFO("=================================================================================================================================");
}

// ********************************************
// DB Grab the Normalisation parameters for the relevant sample name
std::vector<const double*> covarianceOsc::GetOscParsFromSampleName(const std::string& SampleName) {
// ********************************************
  std::vector<const double*> returnVec;
  for (int i = 0; i < _fNumPar; ++i) {
    if (AppliesToSample(i, SampleName)) {
      returnVec.push_back(retPointer(i));
    }
  }
  return returnVec;
}
