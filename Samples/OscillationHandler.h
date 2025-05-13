#pragma once

//MaCh3 includes
#include "Samples/Structs.h"

//forward declare so we don't bleed NuOscillator headers
class OscillatorBase;

/// @brief Interface between NuOscillator and MaCh3, meant to compute oscillation weights for events/bin
/// @author Dan Barrow
class OscillationHandler
{
 public:
   OscillationHandler();
  /// @brief destructor
  virtual ~OscillationHandler();



  /// @brief flag used to define whether all oscillation channels have a probability calculated using the same binning
  bool EqualBinningPerOscChannel = false;
  /// If using shared NuOsc
  bool SharedNuOsc = false;

  //===============================================================================
  /// DB Variables required for oscillation
  std::vector<OscillatorBase*> NuOscProbCalcers;
  std::string NuOscillatorConfigFile;
  //===============================================================================
};
