#pragma once

#include <cmath>

/// @file OscillationProbability.h
/// @brief CUDAProb3Linear's oscillation probability (three flavours, constant density, straight
/// baseline), for one neutrino at a time, on the host and on the device, with exact derivatives.
/// @details A line by line port of NuOscillator's OscProbCalcerCUDAProb3Linear::CalculateProbabilities
/// and cudaprob3linear::physics::get_transition_matrix (setMNSMatrix, setNeutrinoMasses, setMixMatrix,
/// prepare_getMfast, getMfast, get_product, getA), same order of operations, in double precision.
///
/// Everything is a template on the number type T:
///  - T = double gives the probability, which SampleHandlerBase checks against NuOscillator's own
///  - T = OscDual gives the probability and its derivative with respect to each of the eight
///    oscillation parameters (OscDual1 just one of them), carried through the very same arithmetic (forward mode automatic
///    differentiation). That is the exact derivative of the computed probability, so a leapfrog
///    using it stays a proper gradient flow, which finite differences wouldn't be.
///
/// The parameters, in CUDAProb3Linear's order: sin2_th12, sin2_th23, sin2_th13, dm2_12, dm2_23,
/// delta_cp, path_length, matter_density.
/// @author Henry Wallace

#ifdef __CUDACC__
#define M3_OSC_HD __host__ __device__
#else
#define M3_OSC_HD
#endif

namespace M3 {
namespace Osc {

/// Number of oscillation parameters CUDAProb3Linear takes
constexpr int kNPars = 8;

// ================================================================================================
// Dual numbers: a value and its derivative with respect to each oscillation parameter
// ================================================================================================

/// @brief A value and its N partial derivatives, with the arithmetic that carries them along. Each
/// derivative is carried on its own (no operation mixes two of them), so derivative k of Dual<1> seeded
/// with tangent k is bit for bit derivative k of Dual<kNPars>: the device runs one tangent per thread
template <int N>
struct Dual {
  double v;
  double d[N];
  M3_OSC_HD Dual() : v(0.) { for (int k = 0; k < N; ++k) d[k] = 0.; }
  M3_OSC_HD Dual(const double Value) : v(Value) { for (int k = 0; k < N; ++k) d[k] = 0.; }
  /// The variable of parameter k: value x, derivative 1 with respect to itself
  M3_OSC_HD static Dual Variable(const double x, const int k) { Dual r(x); r.d[k] = 1.; return r; }
  M3_OSC_HD Dual& operator+=(const Dual& b) { v += b.v; for (int k = 0; k < N; ++k) d[k] += b.d[k]; return *this; }
  M3_OSC_HD Dual& operator-=(const Dual& b) { v -= b.v; for (int k = 0; k < N; ++k) d[k] -= b.d[k]; return *this; }
  M3_OSC_HD Dual& operator*=(const Dual& b) {
    for (int k = 0; k < N; ++k) d[k] = d[k] * b.v + v * b.d[k];
    v *= b.v;
    return *this;
  }
  M3_OSC_HD Dual& operator/=(const Dual& b) {
    const double q = v / b.v;
    for (int k = 0; k < N; ++k) d[k] = (d[k] - q * b.d[k]) / b.v;
    v = q;
    return *this;
  }
};
/// All eight derivatives at once (the host)
using OscDual = Dual<kNPars>;
/// One derivative (the device, one thread per event and parameter)
using OscDual1 = Dual<1>;

template <int N> M3_OSC_HD inline Dual<N> operator-(const Dual<N>& a) { Dual<N> r(a); r.v = -r.v; for (int k = 0; k < N; ++k) r.d[k] = -r.d[k]; return r; }
template <int N> M3_OSC_HD inline Dual<N> operator+(Dual<N> a, const Dual<N>& b) { a += b; return a; }
template <int N> M3_OSC_HD inline Dual<N> operator-(Dual<N> a, const Dual<N>& b) { a -= b; return a; }
template <int N> M3_OSC_HD inline Dual<N> operator*(Dual<N> a, const Dual<N>& b) { a *= b; return a; }
template <int N> M3_OSC_HD inline Dual<N> operator/(Dual<N> a, const Dual<N>& b) { a /= b; return a; }
template <int N> M3_OSC_HD inline Dual<N> operator+(Dual<N> a, const double b) { a.v += b; return a; }
template <int N> M3_OSC_HD inline Dual<N> operator+(const double a, Dual<N> b) { b.v += a; return b; }
template <int N> M3_OSC_HD inline Dual<N> operator-(Dual<N> a, const double b) { a.v -= b; return a; }
template <int N> M3_OSC_HD inline Dual<N> operator-(const double a, const Dual<N>& b) { return Dual<N>(a) - b; }
template <int N> M3_OSC_HD inline Dual<N> operator*(Dual<N> a, const double b) { a.v *= b; for (int k = 0; k < N; ++k) a.d[k] *= b; return a; }
template <int N> M3_OSC_HD inline Dual<N> operator*(const double a, Dual<N> b) { return b * a; }
template <int N> M3_OSC_HD inline Dual<N> operator/(Dual<N> a, const double b) { a.v /= b; for (int k = 0; k < N; ++k) a.d[k] /= b; return a; }
template <int N> M3_OSC_HD inline Dual<N> operator/(const double a, const Dual<N>& b) { return Dual<N>(a) / b; }

/// Chain rule: f(a) with f'(a) = fp
template <int N>
M3_OSC_HD inline Dual<N> Chain(const Dual<N>& a, const double f, const double fp) {
  Dual<N> r(f);
  for (int k = 0; k < N; ++k) r.d[k] = fp * a.d[k];
  return r;
}

// The functions the probability needs, for both number types. The value parts are exactly the double versions
M3_OSC_HD inline double Value(const double x) { return x; }
template <int N> M3_OSC_HD inline double Value(const Dual<N>& x) { return x.v; }
M3_OSC_HD inline double Sqrt(const double x) { return sqrt(x); }
template <int N> M3_OSC_HD inline Dual<N> Sqrt(const Dual<N>& x) { const double s = sqrt(x.v); return Chain(x, s, 0.5 / s); }
M3_OSC_HD inline double Sin(const double x) { return sin(x); }
template <int N> M3_OSC_HD inline Dual<N> Sin(const Dual<N>& x) { return Chain(x, sin(x.v), cos(x.v)); }
M3_OSC_HD inline double Cos(const double x) { return cos(x); }
template <int N> M3_OSC_HD inline Dual<N> Cos(const Dual<N>& x) { return Chain(x, cos(x.v), -sin(x.v)); }
M3_OSC_HD inline double Asin(const double x) { return asin(x); }
template <int N> M3_OSC_HD inline Dual<N> Asin(const Dual<N>& x) { return Chain(x, asin(x.v), 1. / sqrt(1. - x.v * x.v)); }
M3_OSC_HD inline double Acos(const double x) { return acos(x); }
template <int N> M3_OSC_HD inline Dual<N> Acos(const Dual<N>& x) { return Chain(x, acos(x.v), -1. / sqrt(1. - x.v * x.v)); }
M3_OSC_HD inline double Fabs(const double x) { return fabs(x); }
template <int N> M3_OSC_HD inline Dual<N> Fabs(const Dual<N>& x) { return (x.v < 0.) ? -x : x; }
M3_OSC_HD inline void SinCos(const double x, double* s, double* c) {
#ifdef __CUDA_ARCH__
  sincos(x, s, c);
#else
  *s = sin(x);
  *c = cos(x);
#endif
}
template <int N>
M3_OSC_HD inline void SinCos(const Dual<N>& x, Dual<N>* s, Dual<N>* c) {
  double sv, cv;
  SinCos(x.v, &sv, &cv);
  *s = Chain(x, sv, cv);
  *c = Chain(x, cv, -sv);
}

// ================================================================================================
// The constants of one calculation, from the oscillation parameters (host or device)
// ================================================================================================

/// @brief PMNS matrix, mass differences and CUDAProb3's precomputed products, for neutrinos ([0]) and
/// antineutrinos ([1], delta_cp flipped). Order is which matter eigenvalue goes with which vacuum one,
/// which is discrete, so it has no derivative
template <typename T>
struct Constants {
  T URe[2][9];
  T UIm[2][9];
  T AX[2][324];
  T DM[9];
  int Order[3];
  /// Matter density times CUDAProb3's density_convert
  T Rho;
  /// Baseline in km
  T PathLength;
};

/// Derivative k of every constant, as one-tangent duals (bit for bit what Dual<1> seeded with k gives)
inline Constants<OscDual1> SliceTangent(const Constants<OscDual>& C, const int k) {
  Constants<OscDual1> S;
  auto Cut = [k](const OscDual& x) { OscDual1 r(x.v); r.d[0] = x.d[k]; return r; };
  for (int t = 0; t < 2; ++t) {
    for (int i = 0; i < 9; ++i) { S.URe[t][i] = Cut(C.URe[t][i]); S.UIm[t][i] = Cut(C.UIm[t][i]); }
    for (int i = 0; i < 324; ++i) S.AX[t][i] = Cut(C.AX[t][i]);
  }
  for (int i = 0; i < 9; ++i) S.DM[i] = Cut(C.DM[i]);
  for (int i = 0; i < 3; ++i) S.Order[i] = C.Order[i];
  S.Rho = Cut(C.Rho);
  S.PathLength = Cut(C.PathLength);
  return S;
}

/// CUDAProb3's 2*sqrt(2)*G_F and (1/2)/(hbar c) in its units, and pi as M_PI
constexpr double kTwoRootTwoGf = 1.52588e-4;
constexpr double kLoEFactor = 2.534;
constexpr double kPi = 3.14159265358979323846;

/// @brief The constants for parameters P[0..7], as NuOscillator and CUDAProb3Linear work them out.
/// Returns false if a mixing parameter is negative (NuOscillator throws there)
template <typename T>
M3_OSC_HD bool MakeConstants(const T* P, Constants<T>& C) {
  for (int i = 0; i < 3; ++i) if (Value(P[i]) < 0) return false;
  const T theta12 = Asin(Sqrt(P[0]));
  const T theta23 = Asin(Sqrt(P[1]));
  const T theta13 = Asin(Sqrt(P[2]));
  const T dm12sq = P[3];
  const T dm23sq = P[4];

  // setNeutrinoMasses
  T mVac[3] = {T(0.0), dm12sq, dm12sq + dm23sq};
  const double delta = 5.0e-9;
  if (Value(dm12sq) == 0.0) mVac[0] -= delta;
  if (Value(dm23sq) == 0.0) mVac[2] += delta;
  C.DM[0 * 3 + 0] = T(0.0); C.DM[1 * 3 + 1] = T(0.0); C.DM[2 * 3 + 2] = T(0.0);
  C.DM[0 * 3 + 1] = mVac[0] - mVac[1]; C.DM[1 * 3 + 0] = -C.DM[0 * 3 + 1];
  C.DM[0 * 3 + 2] = mVac[0] - mVac[2]; C.DM[2 * 3 + 0] = -C.DM[0 * 3 + 2];
  C.DM[1 * 3 + 2] = mVac[1] - mVac[2]; C.DM[2 * 3 + 1] = -C.DM[1 * 3 + 2];

  // setMNSMatrix for each neutrino type, then setMixMatrix's precomputed factors
  for (int t = 0; t < 2; ++t) {
    T dCP = P[5];
    if (t == 1) dCP = dCP * -1.;
    const T s12 = Sin(theta12), s13 = Sin(theta13), s23 = Sin(theta23);
    const T c12 = Cos(theta12), c13 = Cos(theta13), c23 = Cos(theta23);
    const T sd = Sin(dCP), cd = Cos(dCP);
    T* re = C.URe[t];
    T* im = C.UIm[t];
    re[0] = c12 * c13;                         im[0] = T(0.0);
    re[1] = s12 * c13;                         im[1] = T(0.0);
    re[2] = s13 * cd;                          im[2] = -s13 * sd;
    re[3] = -s12 * c23 - c12 * s23 * s13 * cd; im[3] = -c12 * s23 * s13 * sd;
    re[4] = c12 * c23 - s12 * s23 * s13 * cd;  im[4] = -s12 * s23 * s13 * sd;
    re[5] = s23 * c13;                         im[5] = T(0.0);
    re[6] = s12 * s23 - c12 * c23 * s13 * cd;  im[6] = -c12 * c23 * s13 * sd;
    re[7] = -c12 * s23 - s12 * c23 * s13 * cd; im[7] = -s12 * c23 * s13 * sd;
    re[8] = c23 * c13;                         im[8] = T(0.0);
    for (int n = 0; n < 3; n++) {
      for (int m = 0; m < 3; m++) {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            T* f = C.AX[t] + n * 3 * 3 * 3 * 4 + m * 3 * 3 * 4 + i * 3 * 4 + j * 4;
            f[0] = re[n * 3 + i] * re[m * 3 + j] + im[n * 3 + i] * im[m * 3 + j];
            f[1] = re[n * 3 + i] * im[m * 3 + j] - im[n * 3 + i] * re[m * 3 + j];
            f[2] = im[n * 3 + i] * im[m * 3 + j] + re[n * 3 + i] * re[m * 3 + j];
            f[3] = im[n * 3 + i] * re[m * 3 + j] - re[n * 3 + i] * im[m * 3 + j];
          }
        }
      }
    }
  }

  // prepare_getMfast: the ordering, from the values alone
  const double DM01 = Value(C.DM[0 * 3 + 1]), DM02 = Value(C.DM[0 * 3 + 2]);
  const double alphaV = DM01 + DM02;
  const double betaV = DM01 * DM02;
  const double gammaV = 0.0;
  double tmpV = alphaV * alphaV - 3.0 * betaV;
  double argV = (2.0 * alphaV * alphaV * alphaV - 9.0 * alphaV * betaV + 27.0 * gammaV) / (2.0 * sqrt(tmpV * tmpV * tmpV));
  if (fabs(argV) > 1.0) argV = argV / fabs(argV);
  const double theta0V = acos(argV) / 3.0;
  const double theta1V = theta0V - (2.0 * kPi / 3.0);
  const double theta2V = theta0V + (2.0 * kPi / 3.0);
  double mMatV[3];
  mMatV[0] = mMatV[1] = mMatV[2] = -(2.0 / 3.0) * sqrt(tmpV);
  mMatV[0] *= cos(theta0V); mMatV[1] *= cos(theta1V); mMatV[2] *= cos(theta2V);
  tmpV = Value(C.DM[0]) - alphaV / 3.0;
  mMatV[0] += tmpV; mMatV[1] += tmpV; mMatV[2] += tmpV;
  for (int i = 0; i < 3; i++) {
    tmpV = fabs(Value(C.DM[i * 3 + 0]) - mMatV[0]);
    int k = 0;
    for (int j = 1; j < 3; j++) {
      const double tmp = fabs(Value(C.DM[i * 3 + 0]) - mMatV[j]);
      if (tmp < tmpV) { k = j; tmpV = tmp; }
    }
    C.Order[i] = k;
  }

  // CUDAProb3's density_convert is 0.5
  C.PathLength = P[6];
  C.Rho = P[7] * 0.5;
  return true;
}

// ================================================================================================
// One probability
// ================================================================================================

template <typename T>
struct Complex { T re, im; };

/// getMfast: the matter mass squared differences, dmMatMat[i][j] = m_i - m_j and dmMatVac[i][j] = m_i - M_j
template <typename T>
M3_OSC_HD void MatterMasses(const Constants<T>& C, const double Enu, const T& rho, const int t, const bool Anti,
                            T dmMatMat[3][3], T dmMatVac[3][3]) {
  const T* DM = C.DM;
  const T* URe = C.URe[t];
  const T* UIm = C.UIm[t];
  const T fac = Anti ? kTwoRootTwoGf * Enu * rho : -kTwoRootTwoGf * Enu * rho;
  const T alpha = fac + DM[0 * 3 + 1] + DM[0 * 3 + 2];
  const T beta = DM[0 * 3 + 1] * DM[0 * 3 + 2] +
                 fac * (DM[0 * 3 + 1] * (1.0 - URe[1] * URe[1] - UIm[1] * UIm[1]) +
                        DM[0 * 3 + 2] * (1.0 - URe[2] * URe[2] - UIm[2] * UIm[2]));
  const T gamma = fac * DM[0 * 3 + 1] * DM[0 * 3 + 2] * (URe[0] * URe[0] + UIm[0] * UIm[0]);
  const T tmp = Value(alpha * alpha - 3.0 * beta) < 0 ? T(0) : alpha * alpha - 3.0 * beta;
  const T argtmp = (2.0 * alpha * alpha * alpha - 9.0 * alpha * beta + 27.0 * gamma) / (2.0 * Sqrt(tmp * tmp * tmp));
  const T arg = (Value(Fabs(argtmp)) > 1.0) ? argtmp / Fabs(argtmp) : argtmp;
  const T theta0 = Acos(arg) / 3.0;
  const T theta1 = theta0 - (2.0 * kPi / 3.0);
  const T theta2 = theta0 + (2.0 * kPi / 3.0);
  T mMatU[3];
  mMatU[0] = -(2.0 / 3.0) * Sqrt(tmp);
  mMatU[1] = -(2.0 / 3.0) * Sqrt(tmp);
  mMatU[2] = -(2.0 / 3.0) * Sqrt(tmp);
  mMatU[0] *= Cos(theta0);
  mMatU[1] *= Cos(theta1);
  mMatU[2] *= Cos(theta2);
  const T tmp2 = DM[0] - alpha / 3.0;
  mMatU[0] += tmp2;
  mMatU[1] += tmp2;
  mMatU[2] += tmp2;
  T mMat[3];
  for (int i = 0; i < 3; i++) mMat[i] = mMatU[C.Order[i]];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dmMatMat[i][j] = mMat[i] - mMat[j];
      dmMatVac[i][j] = mMat[i] - DM[j * 3 + 0];
    }
  }
}

/// get_product followed by getA: the transition amplitude matrix A
template <typename T>
M3_OSC_HD void TransitionMatrix(const Constants<T>& C, const T& L, const double E, const T& rho, const int t, const bool Anti,
                                const T dmMatVac[3][3], const T dmMatMat[3][3], Complex<T> A[3][3]) {
  const T* URe = C.URe[t];
  const T* UIm = C.UIm[t];
  const T fac = Anti ? kTwoRootTwoGf * E * rho : -kTwoRootTwoGf * E * rho;

  // get_product: the matrix 2EH - M_j and the products of Eq. (11)
  Complex<T> twoEHmM[3][3][3];
  for (int n = 0; n < 3; n++) {
    for (int m = 0; m < 3; m++) {
      const T re = -fac * (URe[n] * URe[m] + UIm[n] * UIm[m]);
      const T im = -fac * (URe[n] * UIm[m] - UIm[n] * URe[m]);
      for (int j = 0; j < 3; j++) { twoEHmM[n][m][j].re = re; twoEHmM[n][m][j].im = im; }
    }
  }
  for (int j = 0; j < 3; j++) {
    twoEHmM[0][0][j].re -= dmMatVac[j][0];
    twoEHmM[1][1][j].re -= dmMatVac[j][1];
    twoEHmM[2][2][j].re -= dmMatVac[j][2];
  }
  Complex<T> product[3][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) { product[i][j][k].re = T(0); product[i][j][k].im = T(0); }
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        product[i][j][0].re += twoEHmM[i][k][1].re * twoEHmM[k][j][2].re - twoEHmM[i][k][1].im * twoEHmM[k][j][2].im;
        product[i][j][0].im += twoEHmM[i][k][1].re * twoEHmM[k][j][2].im + twoEHmM[i][k][1].im * twoEHmM[k][j][2].re;
        product[i][j][1].re += twoEHmM[i][k][2].re * twoEHmM[k][j][0].re - twoEHmM[i][k][2].im * twoEHmM[k][j][0].im;
        product[i][j][1].im += twoEHmM[i][k][2].re * twoEHmM[k][j][0].im + twoEHmM[i][k][2].im * twoEHmM[k][j][0].re;
        product[i][j][2].re += twoEHmM[i][k][0].re * twoEHmM[k][j][1].re - twoEHmM[i][k][0].im * twoEHmM[k][j][1].im;
        product[i][j][2].im += twoEHmM[i][k][0].re * twoEHmM[k][j][1].im + twoEHmM[i][k][0].im * twoEHmM[k][j][1].re;
      }
      product[i][j][0].re /= (dmMatMat[0][1] * dmMatMat[0][2]);
      product[i][j][0].im /= (dmMatMat[0][1] * dmMatMat[0][2]);
      product[i][j][1].re /= (dmMatMat[1][2] * dmMatMat[1][0]);
      product[i][j][1].im /= (dmMatMat[1][2] * dmMatMat[1][0]);
      product[i][j][2].re /= (dmMatMat[2][0] * dmMatMat[2][1]);
      product[i][j][2].im /= (dmMatMat[2][0] * dmMatMat[2][1]);
    }
  }

  // getA: the sum with the exponential factors, Eq. (11), then Eq. (10)
  Complex<T> X[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) { X[i][j].re = T(0); X[i][j].im = T(0); }
  }
  for (int k = 0; k < 3; k++) {
    // CUDAProb3 adds a phase offset of 0 to the last term, which changes nothing
    const T arg = -kLoEFactor * dmMatVac[k][0] * L / E;
    T c, sn;
    SinCos(arg, &sn, &c);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        X[i][j].re += c * product[i][j][k].re - sn * product[i][j][k].im;
        X[i][j].im += c * product[i][j][k].im + sn * product[i][j][k].re;
      }
    }
  }
  const T* AX = C.AX[t];
  for (int n = 0; n < 3; n++) {
    for (int m = 0; m < 3; m++) {
      A[n][m].re = T(0);
      A[n][m].im = T(0);
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          const T* f = AX + n * 3 * 3 * 3 * 4 + m * 3 * 3 * 4 + i * 3 * 4 + j * 4;
          A[n][m].re += f[0] * X[i][j].re + f[1] * X[i][j].im;
          A[n][m].im += f[2] * X[i][j].im + f[3] * X[i][j].re;
        }
      }
    }
  }
}

/// The (n, m) element of the 2EH - M_j matrices of get_product, from the shared part Base and the
/// diagonal's dmMatVac, as TransitionMatrix fills twoEHmM[n][m][j] (the same operations)
template <typename T>
M3_OSC_HD Complex<T> TwoEHmM(const Complex<T> Base[3][3], const T dmMatVac[3][3], const int n, const int m, const int j) {
  Complex<T> r = Base[n][m];
  if (n == m) r.re -= dmMatVac[j][n];
  return r;
}

/// The one amplitude A[n][m] of TransitionMatrix, with every product, X and A term computed as
/// TransitionMatrix computes it and summed in the same order, so the result is the same to the bit.
/// Nothing but the 3x3 shared matrix, three phases and the accumulator is kept: TransitionMatrix holds
/// the 27 products, X and all nine amplitudes, which in dual numbers spill out of registers
template <typename T>
M3_OSC_HD Complex<T> Amplitude(const Constants<T>& C, const T& L, const double E, const T& rho, const int t, const bool Anti,
                               const T dmMatVac[3][3], const T dmMatMat[3][3], const int n, const int m) {
  const T* URe = C.URe[t];
  const T* UIm = C.UIm[t];
  const T fac = Anti ? kTwoRootTwoGf * E * rho : -kTwoRootTwoGf * E * rho;
  Complex<T> Base[3][3];
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      Base[a][b].re = -fac * (URe[a] * URe[b] + UIm[a] * UIm[b]);
      Base[a][b].im = -fac * (URe[a] * UIm[b] - UIm[a] * URe[b]);
    }
  }
  T c[3], sn[3];
  for (int k = 0; k < 3; k++) {
    const T arg = -kLoEFactor * dmMatVac[k][0] * L / E;
    SinCos(arg, &sn[k], &c[k]);
  }
  const T Denom[3] = {dmMatMat[0][1] * dmMatMat[0][2], dmMatMat[1][2] * dmMatMat[1][0], dmMatMat[2][0] * dmMatMat[2][1]};
  // product[i][j][p] multiplies matrices J1[p] and J2[p]: (1, 2), (2, 0), (0, 1)
  const int J1[3] = {1, 2, 0};
  const int J2[3] = {2, 0, 1};
  const T* AX = C.AX[t];
  Complex<T> Amp;
  Amp.re = T(0);
  Amp.im = T(0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Complex<T> X;
      X.re = T(0);
      X.im = T(0);
      for (int p = 0; p < 3; p++) {
        Complex<T> prod;
        prod.re = T(0);
        prod.im = T(0);
        for (int k = 0; k < 3; k++) {
          const Complex<T> a = TwoEHmM(Base, dmMatVac, i, k, J1[p]);
          const Complex<T> b = TwoEHmM(Base, dmMatVac, k, j, J2[p]);
          prod.re += a.re * b.re - a.im * b.im;
          prod.im += a.re * b.im + a.im * b.re;
        }
        prod.re /= Denom[p];
        prod.im /= Denom[p];
        X.re += c[p] * prod.re - sn[p] * prod.im;
        X.im += c[p] * prod.im + sn[p] * prod.re;
      }
      const T* f = AX + n * 3 * 3 * 3 * 4 + m * 3 * 3 * 4 + i * 3 * 4 + j * 4;
      Amp.re += f[0] * X.re + f[1] * X.im;
      Amp.im += f[2] * X.im + f[3] * X.re;
    }
  }
  return Amp;
}

/// @brief The probability of From -> To (0 = e, 1 = mu, 2 = tau) at energy E (GeV), for a neutrino
/// (Anti false) or antineutrino, clamped to [0, 1] as NuOscillator clamps CUDAProb3Linear's output
/// (a clamped probability has zero derivative). CUDAProb3Linear's "from -> to" is |A[to][from]|^2
/// (calculate_beam fills result[from * 3 + to] with Prob[to][from])
template <typename T>
M3_OSC_HD T ProbabilityFromAmplitude(const Complex<T>& a) {
  const T P = a.re * a.re + a.im * a.im;
  if (!(Value(P) > 0.0)) return T(0.0);
  if (!(Value(P) < 1.0)) return T(1.0);
  return P;
}

template <typename T>
M3_OSC_HD T Probability(const Constants<T>& C, const double E, const bool Anti, const int From, const int To) {
  const int t = Anti ? 1 : 0;
  T dmMatMat[3][3], dmMatVac[3][3];
  MatterMasses(C, E, C.Rho, t, Anti, dmMatMat, dmMatVac);
  return ProbabilityFromAmplitude(Amplitude(C, C.PathLength, E, C.Rho, t, Anti, dmMatVac, dmMatMat, To, From));
}

/// @brief Every channel's amplitude at one energy: A[to][from], the same numbers Probability's
/// Amplitude gives one at a time (TransitionMatrix, same order of operations). Much cheaper when
/// several channels share an energy, as binned oscillation's weights do
template <typename T>
M3_OSC_HD void AllAmplitudes(const Constants<T>& C, const double E, const bool Anti, Complex<T> A[3][3]) {
  const int t = Anti ? 1 : 0;
  T dmMatMat[3][3], dmMatVac[3][3];
  MatterMasses(C, E, C.Rho, t, Anti, dmMatMat, dmMatVac);
  TransitionMatrix(C, C.PathLength, E, C.Rho, t, Anti, dmMatVac, dmMatMat, A);
}

}  // namespace Osc
}  // namespace M3
