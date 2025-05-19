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
  /// @brief Constructor
  /// @param ConfigFile name/path to NuOscillator Config
  /// @param EqualBinningPerChannel whether to use same binning per each oscillation channel
  /// @param OscParams_ Pointers to values of oscillation parameters
  /// @param SubChannels Number of oscillation channels
  /// @warning if @EqualBinningPerChannel is true then argument SubChannels make no difference
  OscillationHandler(const std::string& ConfigFile, bool EqualBinningPerChannel,
                     std::vector<const double*> OscParams_, const int SubChannels);
  /// @brief Destructor
  virtual ~OscillationHandler();
  /// @brief check if same binning is used for multiple oscillation channels
  bool isEqualBinningPerOscChannel() {return EqualBinningPerOscChannel;}
  /// @brief DB Evaluate oscillation weights for each defined event/bin
  void Evaluate();
  /// @brief Get pointer to oscillation weight
  const M3::float_t* GetNuOscillatorPointers(int Channel, int InitFlav, int FinalFlav, FLOAT_T TrueEnu, FLOAT_T TrueCosZenith = -999);

  /// @brief Setup binning, arrays correspond to events and their energy bins
  void SetOscillatorBinning(const int Channel, const std::vector<M3::float_t>& EnergyArray, const std::vector<M3::float_t>& CosineZArray);

  /// @brief return size of oscillation parameter pointer vector
  unsigned int GetOscParamsSize() const {return static_cast<unsigned int>(OscParams.size());};
 private:
  /// flag used to define whether all oscillation channels have a probability calculated using the same binning
  bool EqualBinningPerOscChannel;

  /// DB Variables required for oscillation
  std::vector<std::unique_ptr<OscillatorBase>> NuOscProbCalcers;

  /// pointer to osc params, since not all params affect every sample, we perform some operations before hand for speed
  std::vector<const double*> OscParams;
};
