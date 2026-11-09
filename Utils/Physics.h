#pragma once

// MaCh3 includes
#include "Manager/Manager.h"
#include "Manager/Random.h"

/// @brief DB Get the Cherenkov momentum threshold in MeV
/// @param PDG PDG code of the particle for which the Cherenkov threshold is requested.
double returnCherenkovThresholdMomentum(const int PDG);

/// @brief Recalculate Q^2 after Eb shift. Takes in shifted lepton momentum, lepton angle, and true neutrino energy
/// @param PLep Shifted outgoing lepton momentum (MeV/c).
/// @param PUpd Updated lepton angle or related kinematic quantity used in the recalculation.
/// @param EnuTrue True neutrino energy (MeV).
/// @param InitialQ2 Optional initial Q^2 value (used as seed or fallback, default = 0).
double CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2 = 0.0);

/// @brief Recalculate Enu after Eb shift. Takes in shifted lepton momentum, lepton angle, and binding energy change, and if nu/anu
/// @param PLep Shifted outgoing lepton momentum (MeV/c).
/// @param cosTheta Cosine of the lepton scattering angle.
/// @param EB Binding energy shift applied to the interaction (MeV).
/// @param neutrino True if the interaction is neutrino, false if antineutrino.
double CalculateEnu(double PLep, double cosTheta, double EB, bool neutrino);
