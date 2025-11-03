#include "OscillationHandler.h"

_MaCh3_Safe_Include_Start_ //{
#include "Oscillator/OscillatorFactory.h"
#include "Constants/OscillatorConstants.h"
_MaCh3_Safe_Include_End_ //}

// ************************************************
OscillationHandler::OscillationHandler(const std::string& NuOscillatorConfigFile, bool BinningPerOscChannel_,
                                       std::vector<const double*> OscParams_, const int SubChannels) {
// ************************************************

  MaCh3_ProfileScope;
  
  EqualBinningPerOscChannel = BinningPerOscChannel_;
  OscParams = OscParams_;
  // Add first sample
  NuOscProbCalcers.resize(1);

  auto OscillFactory = std::make_unique<OscillatorFactory>();
  //DB's explanation of EqualBinningPerOscChannel:
  //In the situation where we are applying binning oscillation probabilities to a SampleHandler object, it maybe the case that there is identical binning per oscillation channel
  //In which case, and remembering that each NuOscillator::Oscillator object calculate the oscillation probabilities for all channels, we just have to create one Oscillator object and use the results from that
  //This means that we can get up to a factor of 12 reduction in the calculation time of the oscillation probabilities, because we don't need to repeat the operation per oscillation channel

  if (EqualBinningPerOscChannel) {
    NuOscProbCalcers[0].resize(1);
    LoggerPrint("NuOscillator",
                [](const std::string& message) { MACH3LOG_INFO("{}", message); },
                [this, &OscillFactory, &NuOscillatorConfigFile]() {
                  this->NuOscProbCalcers[0][0] = std::unique_ptr<OscillatorBase>(OscillFactory->CreateOscillator(NuOscillatorConfigFile));
                });

    if (!NuOscProbCalcers[0][0]->EvalPointsSetInConstructor()) {
      MACH3LOG_ERROR("Attempted to use equal binning per oscillation channel, but not binning has been set in the NuOscillator::Oscillator object");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    NuOscProbCalcers[0][0]->Setup();
  } else {
    NuOscProbCalcers[0].resize(SubChannels);
    for (int iChannel = 0; iChannel < SubChannels; iChannel++) {
      MACH3LOG_INFO("Setting up NuOscillator::Oscillator object in OscillationChannel: {}/{}", iChannel, SubChannels);

      LoggerPrint("NuOscillator",
                  [](const std::string& message) { MACH3LOG_INFO("{}", message); },
                  [this, iChannel, &OscillFactory, &NuOscillatorConfigFile]() {
                    this->NuOscProbCalcers[0][iChannel] = std::unique_ptr<OscillatorBase>(
                      OscillFactory->CreateOscillator(NuOscillatorConfigFile));
                  });
    }
  }
}

// ************************************************
OscillationHandler::~OscillationHandler() {
// ************************************************

  MaCh3_ProfileScope;

}

// ************************************************
void OscillationHandler::AddSample(const std::string& NuOscillatorConfigFile, const int SubChannels) {
// ************************************************
  if(EqualBinningPerOscChannel){
    MACH3LOG_ERROR("Trying to add sample while EqualBinningPerOscChannel is enabled. This will not work...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  auto OscillFactory = std::make_unique<OscillatorFactory>();

  std::vector<std::unique_ptr<OscillatorBase>> OscProbCalcersTemp(SubChannels);

  for (int iChannel = 0; iChannel < SubChannels; iChannel++) {
    MACH3LOG_INFO("Setting up NuOscillator::Oscillator object in OscillationChannel: {}/{}", iChannel, SubChannels);

    LoggerPrint("NuOscillator",
                [](const std::string& message) { MACH3LOG_INFO("{}", message); },
                [&OscProbCalcersTemp, iChannel, &OscillFactory, &NuOscillatorConfigFile]() {
                  OscProbCalcersTemp[iChannel] = std::unique_ptr<OscillatorBase>(
                    OscillFactory->CreateOscillator(NuOscillatorConfigFile));
                });
  }
  NuOscProbCalcers.push_back(std::move(OscProbCalcersTemp));
}

// ************************************************
void OscillationHandler::Evaluate() {
// ************************************************

  MaCh3_ProfileScope;
  
  std::vector<M3::float_t> OscVec(OscParams.size());
  for (size_t iPar = 0; iPar < OscParams.size(); ++iPar) {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wuseless-cast"
    OscVec[iPar] = static_cast<M3::float_t>(*OscParams[iPar]);
    #pragma GCC diagnostic pop
  }

  if (EqualBinningPerOscChannel) {
    NuOscProbCalcers[0][0]->CalculateProbabilities(OscVec);
  } else {
    for (size_t iSample = 0; iSample < NuOscProbCalcers.size(); iSample++) {
      for (size_t iChannel = 0; iChannel < NuOscProbCalcers[iSample].size(); iChannel++) {
        NuOscProbCalcers[iSample][iChannel]->CalculateProbabilities(OscVec);
      }
    }
  }
}


// ************************************************
const M3::float_t* OscillationHandler::GetNuOscillatorPointers(const int Sample, const int Channel, const int InitFlav,
                                                               const int FinalFlav, const FLOAT_T TrueEnu, const FLOAT_T TrueCosZenith) {
// ************************************************

  MaCh3_ProfileScope;
  
  int IndexSample = 0;
  int IndexChannel = 0;
  if (!EqualBinningPerOscChannel) {
    IndexSample = Sample;
    IndexChannel = Channel;
  }

  if(TrueCosZenith != -999) {
    return NuOscProbCalcers[IndexSample][IndexChannel]->ReturnWeightPointer(InitFlav ,FinalFlav, TrueEnu, TrueCosZenith);
  } else {
    return NuOscProbCalcers[IndexSample][IndexChannel]->ReturnWeightPointer(InitFlav ,FinalFlav, TrueEnu);
  }
}


// ************************************************
void OscillationHandler::SetOscillatorBinning(const int Sample, const int Channel, const std::vector<M3::float_t>& EnergyArray,
                                              const std::vector<M3::float_t>& CosineZArray) {
// ************************************************

  MaCh3_ProfileScope;
  
  if (!NuOscProbCalcers[Sample][Channel]->EvalPointsSetInConstructor()) {
    NuOscProbCalcers[Sample][Channel]->SetEnergyArrayInCalcer(EnergyArray);
    if(CosineZArray.size() != 0) NuOscProbCalcers[Sample][Channel]->SetCosineZArrayInCalcer(CosineZArray);
  }
  NuOscProbCalcers[Sample][Channel]->Setup();
}
