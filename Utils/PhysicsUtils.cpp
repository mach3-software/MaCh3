#include "Utils/PhysicsUtils.h"

// ****************
//DB Get the Cherenkov momentum threshold in MeV
double returnCherenkovThresholdMomentum(const int PDG) {
// ****************
  constexpr double refractiveIndex = 1.334; //DB From https://github.com/fiTQun/fiTQun/blob/646cf9c8ba3d4f7400bcbbde029d5ca15513a3bf/fiTQun_shared.cc#L757
  double mass = M3::Utils::GetMassFromPDG(PDG)*1e3;
  double momentumThreshold = mass/sqrt(refractiveIndex*refractiveIndex-1.);
  return momentumThreshold;
}

// **************************************************************************
// Recalculate Q^2 after Eb shift. Takes in shifted lepton momentum, lepton angle, and true neutrino energy
double CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2){
// ***************************************************************************
  constexpr double MLep = 0.10565837;
  constexpr double MLep2 = MLep * MLep;

  // Calculate muon energy
  double ELep = std::sqrt((MLep2)+(PLep*PLep));

  double CosTh = (2*EnuTrue*ELep - MLep2 - InitialQ2)/(2*EnuTrue*PLep);

  ELep = std::sqrt((MLep2)+(PUpd*PUpd));

  // Calculate the new Q2
  double Q2Upd = -(MLep2) + 2.0*EnuTrue*(ELep - PUpd*CosTh);

  return Q2Upd - InitialQ2;
}

// **************************************************************************
// Recalculate Enu after Eb shift. Takes in shifted lepton momentum, lepton angle, and binding energy change, and if nu/anu
double CalculateEnu(double PLep, double costh, double Eb, bool neutrino){
// ***************************************************************************
  double mNeff = 0.93956536 - Eb / 1000.;
  double mNoth = 0.93827203;

  if (!neutrino) {
    mNeff = 0.93827203 - Eb / 1000.;
    mNoth = 0.93956536;
  }
  /// @todo WARNING this is hardcoded
  constexpr double mLep = 0.10565837;
  constexpr double mLep2 = mLep * mLep;

  const double eLep = sqrt(PLep * PLep + mLep2);
  double Enu = (2 * mNeff * eLep - mLep2 + mNoth * mNoth - mNeff * mNeff) /(2 * (mNeff - eLep + PLep * costh));

  return Enu;
}
