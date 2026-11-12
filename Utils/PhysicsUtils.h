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

namespace M3 {
namespace Utils {
// *****************************
/// @brief Return mass for given PDG
/// @note Get the mass of a particle from the PDG In GeV, not MeV!
/// @cite pdg2024 (particle masses)
/// @cite ame2020 (nuclear masses)
inline constexpr double GetMassFromPDG(const int PDG) {
// *****************************
  switch (abs(PDG)) {
    // Leptons
    case 11: return 0.00051099895; // e
    case 13: return 0.1056583755;  // mu
    case 15: return 1.77693;       // tau
    // Neutrinos
    case 12:
    case 14:
    case 16:
      return 0.;
      // Photon
    case 22: return 0.;
    // Mesons
    case 211: return 0.13957039; // pi_+/-
    case 111: return 0.1349768;  // pi_0
    case 221: return 0.547862;   // eta
    case 331: return 0.95778;    // eta'
    case 311:                    // K_0
    case 130:                    // K_0_L
    case 310:                    // K_0_S
      return 0.497611;
    case 321: return 0.493677;   // K_+/-
    case 113: return 0.77526;    // rho_0
    case 213: return 0.77511;    // rho_+/-
    case 223: return 0.78266;    // omega
    case 411: return 1.86966;    // D_+/-
    case 421: return 1.86484;    // D_0
    case 431: return 1.96835;    // D_s+/-
    // Baryons
    case 2112: return 0.939565; // n
    case 2212: return 0.938272; // p
    case 3122: return 1.115683; // Lambda
    case 3222: return 1.118937; // Sig_+
    case 3112: return 1.197449; // Sig_-
    case 3212: return 1.192642; // Sig_0
    case 3312: return 1.32171;  // Xi_+/-
    case 3322: return 1.31486;  // Xi_0
    case 3334: return 1.67245;  // Omega_+/-
    case 4122: return 2.28646;  // Lambda_c+
    case 4212: return 2.45265;  // Sigma_c+
    case 4222: return 2.45397;  // Sigma_c++
    // Nuclei
    case 1000050110: return 10.255103; // Boron-11
    case 1000060120: return 11.177929; // Carbon-12
    case 1000070140: return 13.043781; // Nitrogen-14
    case 1000080160: return 14.899169; // Oxygen-16
    case 1000090190: return 17.696901; // Fluorine-19
    case 1000110230: return 21.414835; // Sodium-23
    case 1000130270: return 25.133144; // Aluminum-27
    case 1000140280: return 26.060342; // Silicon-28
    case 1000190390: return 36.294463; // Potassium-39
    case 1000180400: return 37.224724; // Argon-40
    case 1000220480: return 44.663224; // Titanium-48
    case 1000300640: return 59.549619; // Zinc-64
    default:
      MACH3LOG_ERROR("Haven't got a saved mass for PDG: {}", PDG);
      MACH3LOG_ERROR("Please implement me!");
      throw MaCh3Exception(__FILE__, __LINE__);
  } // End switch
  return 0;
}
} // end Utils
} // end M3
