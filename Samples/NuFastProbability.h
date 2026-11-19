#pragma once

#include "Samples/OscillationProbability.h"

/// @file NuFastProbability.h
/// @brief NuFAST's long-baseline oscillation probabilities in constant matter (Denton and Parke,
/// Probability_Matter_LBL in NuFast.cpp, which NuOscillator's NuFASTLinear calls), as a template on
/// the number type, so dual numbers carry the derivatives with respect to the oscillation parameters
/// through it.
/// @details A line by line port, same order of operations: with T = double it gives NuOscillator's
/// NuFASTLinear weights, which OscillationGradient checks. All nine channels come out of one call at
/// one energy, from an approximate matter eigenvalue improved by a fixed number of Newton steps and a
/// few dozen products: far cheaper than CUDAProb3's transition matrix, in dual numbers especially.
/// There are no branches on the parameters (only on N_Newton), so the derivative is that of the
/// computed probability everywhere.
/// @author Henry Wallace

namespace M3 {
namespace Osc {

/// NuFast.cpp's constants
constexpr double kNuFastLover4E = 1e-9 / 1.97327e-7 * 1e3 / 4;
constexpr double kNuFastYerhoE2a = 1.52588e-4;

/// @brief All nine probabilities, P[from][to] (0 e, 1 mu, 2 tau), as NuFast.cpp's probs_returned
/// @param E Energy in GeV, negative for antineutrinos (NuFAST's convention)
/// @param N_Newton Newton iterations on lambda3; must be zero or more (the exact, acos branch isn't ported)
template <typename T>
M3_OSC_HD void NuFastProbabilities(const T& s12sq, const T& s13sq, const T& s23sq, const T& delta, const T& Dmsq21,
                                   const T& Dmsq31, const T& L, const double E, const T& rho, const T& Ye,
                                   const int N_Newton, T P[3][3]) {
  T c13sq, sind, cosd, Jrr, Jmatter, Dmsqee, Amatter;
  T Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq;
  T A, B, C;
  T See, Tee, Smm, Tmm;
  T xmat, lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32;
  T Xp2, Xp3, PiDlambdaInv, tmp;
  T Lover4E, D21, D32;
  T sinD21, sinD31, sinD32;
  T sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin;
  T Pme_CPC, Pme_CPV, Pmm, Pee;

  c13sq = 1. - s13sq;
  Ue2sq = c13sq * s12sq;
  Ue3sq = s13sq;
  Um3sq = c13sq * s23sq;
  Ut2sq = s13sq * s12sq * s23sq;
  Um2sq = (1. - s12sq) * (1. - s23sq);

  Jrr = Sqrt(Um2sq * Ut2sq);
  sind = Sin(delta);
  cosd = Cos(delta);

  Um2sq = Um2sq + Ut2sq - 2. * Jrr * cosd;
  Jmatter = 8. * Jrr * c13sq * sind;
  Amatter = Ye * rho * E * kNuFastYerhoE2a;
  Dmsqee = Dmsq31 - s12sq * Dmsq21;

  A = Dmsq21 + Dmsq31;
  See = A - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq;
  Tmm = Dmsq21 * Dmsq31;
  Tee = Tmm * (1. - Ue3sq - Ue2sq);
  C = Amatter * Tee;
  A = A + Amatter;

  xmat = Amatter / Dmsqee;
  tmp = 1. - xmat;
  lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1. + Sqrt(tmp * tmp + 4. * s13sq * xmat));

  B = Tmm + Amatter * See;
  for (int i = 0; i < N_Newton; i++)
    lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2. * (lambda3 - A) + lambda3) + B);

  tmp = A - lambda3;
  Dlambda21 = Sqrt(tmp * tmp - 4. * C / lambda3);
  lambda2 = 0.5 * (A - lambda3 + Dlambda21);
  Dlambda32 = lambda3 - lambda2;
  Dlambda31 = Dlambda32 + Dlambda21;

  PiDlambdaInv = 1. / (Dlambda31 * Dlambda32 * Dlambda21);
  Xp3 = PiDlambdaInv * Dlambda21;
  Xp2 = -PiDlambdaInv * Dlambda31;

  Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
  Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

  Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq;
  Tmm = Tmm * (1. - Um3sq - Um2sq) + Amatter * (See + Smm - A);

  Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
  Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

  Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv;

  Ue1sq = 1. - Ue3sq - Ue2sq;
  Um1sq = 1. - Um3sq - Um2sq;
  Ut3sq = 1. - Um3sq - Ue3sq;
  Ut2sq = 1. - Um2sq - Ue2sq;
  Ut1sq = 1. - Um1sq - Ue1sq;

  Lover4E = kNuFastLover4E * L / E;
  D21 = Dlambda21 * Lover4E;
  D32 = Dlambda32 * Lover4E;
  sinD21 = Sin(D21);
  sinD31 = Sin(D32 + D21);
  sinD32 = Sin(D32);
  triple_sin = sinD21 * sinD31 * sinD32;
  sinsqD21_2 = 2. * sinD21 * sinD21;
  sinsqD31_2 = 2. * sinD31 * sinD31;
  sinsqD32_2 = 2. * sinD32 * sinD32;

  Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
          + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
          + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
  Pme_CPV = -Jmatter * triple_sin;
  Pmm = 1. - 2. * (Um2sq * Um1sq * sinsqD21_2 + Um3sq * Um1sq * sinsqD31_2 + Um3sq * Um2sq * sinsqD32_2);
  Pee = 1. - 2. * (Ue2sq * Ue1sq * sinsqD21_2 + Ue3sq * Ue1sq * sinsqD31_2 + Ue3sq * Ue2sq * sinsqD32_2);

  P[0][0] = Pee;
  P[0][1] = Pme_CPC - Pme_CPV;
  P[0][2] = 1. - Pee - P[0][1];
  P[1][0] = Pme_CPC + Pme_CPV;
  P[1][1] = Pmm;
  P[1][2] = 1. - P[1][0] - Pmm;
  P[2][0] = 1. - Pee - P[1][0];
  P[2][1] = 1. - P[0][1] - Pmm;
  P[2][2] = 1. - P[0][2] - P[1][2];
}

}  // namespace Osc
}  // namespace M3
