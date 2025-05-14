#include "OscillationHandler.h"

_MaCh3_Safe_Include_Start_ //{
#include "Oscillator/OscillatorFactory.h"
#include "Constants/OscillatorConstants.h"
_MaCh3_Safe_Include_End_ //}

// ************************************************
OscillationHandler::OscillationHandler(const std::string& NuOscillatorConfigFile, bool BinningPerOscChannel_,
                                       std::vector<const double*> OscParams_, const int SubChannels) {
// ************************************************
  EqualBinningPerOscChannel = BinningPerOscChannel_;
  OscParams = OscParams_;
  auto OscillFactory = std::make_unique<OscillatorFactory>();

  //DB's explanation of EqualBinningPerOscChannel:
  //In the situation where we are applying binning oscillation probabilities to a SampleHandler object, it maybe the case that there is identical binning per oscillation channel
  //In which case, and remembering that each NuOscillator::Oscillator object calculate the oscillation probabilities for all channels, we just have to create one Oscillator object and use the results from that
  //This means that we can get up to a factor of 12 reduction in the calculation time of the oscillation probabilities, because we don't need to repeat the operation per oscillation channel

  if (EqualBinningPerOscChannel) {
    NuOscProbCalcers.resize(1);
    LoggerPrint("NuOscillator",
                [](const std::string& message) { MACH3LOG_INFO("{}", message); },
                [this, &OscillFactory, &NuOscillatorConfigFile]() {
                  this->NuOscProbCalcers[0] = std::unique_ptr<OscillatorBase>(OscillFactory->CreateOscillator(NuOscillatorConfigFile));
                });

    if (!NuOscProbCalcers[0]->EvalPointsSetInConstructor()) {
      MACH3LOG_ERROR("Attempted to use equal binning per oscillation channel, but not binning has been set in the NuOscillator::Oscillator object");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    NuOscProbCalcers[0]->Setup();
  } else {
    NuOscProbCalcers.resize(SubChannels);
    for (int iSample = 0; iSample < SubChannels; iSample++) {
      MACH3LOG_INFO("Setting up NuOscillator::Oscillator object in OscillationChannel: {}/{}", iSample, SubChannels);

      LoggerPrint("NuOscillator",
                  [](const std::string& message) { MACH3LOG_INFO("{}", message); },
                  [this, iSample, &OscillFactory, &NuOscillatorConfigFile]() {
                    this->NuOscProbCalcers[iSample] = std::unique_ptr<OscillatorBase>(
                      OscillFactory->CreateOscillator(NuOscillatorConfigFile));
                  });
    }
  }
}

// ************************************************
OscillationHandler::~OscillationHandler() {
// ************************************************

}

// ************************************************
void OscillationHandler::Evaluate() {
// ************************************************
  std::vector<M3::float_t> OscVec(OscParams.size());
  for (size_t iPar = 0; iPar < OscParams.size(); ++iPar) {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wuseless-cast"
    OscVec[iPar] = M3::float_t(*OscParams[iPar]);
    #pragma GCC diagnostic pop
  }

  if (EqualBinningPerOscChannel) {
    NuOscProbCalcers[0]->CalculateProbabilities(OscVec);
  } else {
    for (size_t iChannel = 0; iChannel < NuOscProbCalcers.size(); iChannel++) {
      NuOscProbCalcers[iChannel]->CalculateProbabilities(OscVec);
    }
  }
}


// ************************************************
const M3::float_t* OscillationHandler::GetNuOscillatorPointers(int Channel, int InitFlav, int FinalFlav, FLOAT_T TrueEnu, FLOAT_T TrueCosZenith) {
// ************************************************
  int Index = 0;
  if (!EqualBinningPerOscChannel) {
    Index = Channel;
  }

  if(TrueCosZenith != -999) {
    return NuOscProbCalcers[Index]->ReturnWeightPointer(InitFlav ,FinalFlav, TrueEnu, TrueCosZenith);
  } else {
    return NuOscProbCalcers[Index]->ReturnWeightPointer(InitFlav ,FinalFlav, TrueEnu);
  }
}


// ************************************************
void OscillationHandler::SetOscillatorBinning(const int Channel, const std::vector<M3::float_t>& EnergyArray, const std::vector<M3::float_t>& CosineZArray) {
// ************************************************
  if (!NuOscProbCalcers[Channel]->EvalPointsSetInConstructor()) {
    NuOscProbCalcers[Channel]->SetEnergyArrayInCalcer(EnergyArray);
    if(CosineZArray.size() != 0) NuOscProbCalcers[Channel]->SetCosineZArrayInCalcer(CosineZArray);
  }
  NuOscProbCalcers[Channel]->Setup();
}
