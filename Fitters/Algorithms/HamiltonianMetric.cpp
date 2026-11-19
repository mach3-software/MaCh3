#include "Fitters/Algorithms/HamiltonianMetric.h"

#include <algorithm>
#include <cmath>

#include "Parameters/ParameterHandlerUtils.h"

// *******************
void MetricBlock::Build(const TMatrixDSym& Throw, const std::vector<char>& Frozen, const double Shrinkage,
                        const std::vector<double>& PriorWidth, const std::string& Name) {
// *******************
  N = Throw.GetNrows();
  TMatrixDSym Inverse(Throw);

  // Keep the variances, which adaption estimates well, and damp the correlations by however much
  // it hasn't resolved them. See HamiltonianMC::metric_shrinkage
  if (Shrinkage > 0.) {
    const double Keep = 1. - Shrinkage;
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < i; ++j) {
        Inverse(i, j) *= Keep;
        Inverse(j, i) *= Keep;
      }
    }
  }

  // A frozen parameter mustn't move, so decouple it and leave a unit mass on its diagonal to keep
  // the matrix decomposable. The sampler zeroes its velocity as well
  for (int i = 0; i < N; ++i) {
    if (!Frozen[static_cast<size_t>(i)]) continue;
    for (int j = 0; j < N; ++j) {
      Inverse(i, j) = 0.;
      Inverse(j, i) = 0.;
    }
    Inverse(i, i) = 1.;
  }

  //KS: Flat row major copies, ROOT's indexing is too slow for the per kick matrix multiply
  const size_t Size = static_cast<size_t>(N);
  //HW Fully shrunk means there is nothing off the diagonal left, so the solves collapse to O(N)
  Diagonal = (Shrinkage >= 1.);
  InvMass.assign(Size * Size, 0.);
  Chol.assign(Size * Size, 0.);

  if (Diagonal) {
    //HW L is the square roots of the variances, so there is nothing to decompose
    int NFallback = 0;
    for (int i = 0; i < N; ++i) {
      double Var = Inverse(i, i);
      if (!std::isfinite(Var) || Var <= 0.) {
        //HW A direction adaption never saw move has no variance to scale by. Fall back to the prior
        //width, which is what the metric would have been without adaption
        const double Width = PriorWidth[static_cast<size_t>(i)];
        Var = (std::isfinite(Width) && Width > 0.) ? Width * Width : 1.;
        ++NFallback;
      }
      const size_t Diag = static_cast<size_t>(i) * Size + static_cast<size_t>(i);
      InvMass[Diag] = Var;
      Chol[Diag] = std::sqrt(Var);
    }
    if (NFallback > 0) {
      MACH3LOG_WARN("{} had {} parameters with no usable variance for the diagonal metric, "
                    "they keep their prior width", Name, NFallback);
    }
  } else {
    if (!M3::CanDecomposeMatrix(Inverse)) {
      MACH3LOG_WARN("Throw matrix of {} isn't positive definite, nudging the diagonal", Name);
      M3::MakeMatrixPosDef(&Inverse);
    }
    const auto L = M3::GetCholeskyDecomposedMatrix(Inverse, Name);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        const size_t Element = static_cast<size_t>(i) * Size + static_cast<size_t>(j);
        InvMass[Element] = Inverse(i, j);
        Chol[Element] = L[static_cast<size_t>(i)][static_cast<size_t>(j)];
      }
    }
  }
  BuildSkyline();
}

// *******************
void MetricBlock::BuildSkyline() {
// *******************
  const size_t Size = static_cast<size_t>(N);
  CholRowFirst.assign(Size, 0);
  CholColLast.assign(Size, 0);
  InvRowFirst.assign(Size, 0);
  InvRowLast.assign(Size, 0);
  for (int i = 0; i < N; ++i) {
    int First = i;
    for (int j = 0; j < i; ++j) if (Chol[static_cast<size_t>(i) * Size + static_cast<size_t>(j)] != 0.) { First = j; break; }
    CholRowFirst[static_cast<size_t>(i)] = First;
    int Last = i;
    for (int j = N - 1; j > i; --j) if (Chol[static_cast<size_t>(j) * Size + static_cast<size_t>(i)] != 0.) { Last = j; break; }
    CholColLast[static_cast<size_t>(i)] = Last;
    int InvFirst = N, InvLast = -1;
    for (int j = 0; j < N; ++j) if (InvMass[static_cast<size_t>(i) * Size + static_cast<size_t>(j)] != 0.) { InvFirst = std::min(InvFirst, j); InvLast = j; }
    InvRowFirst[static_cast<size_t>(i)] = (InvLast < 0) ? 0 : InvFirst;
    InvRowLast[static_cast<size_t>(i)] = InvLast;
  }
}

// *******************
void MetricBlock::ThrowVelocity(const double* _restrict_ Throws, double* _restrict_ Velocity) const {
// *******************
  const double* _restrict_ L = Chol.data();
  const int NPars = N;
  if (Diagonal) {
    //HW L is diagonal, so each throw is simply scaled by its own standard deviation
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < NPars; ++i) Velocity[i] = L[i * NPars + i] * Throws[i];
    return;
  }
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < NPars; ++i) {
    double Sum = 0.;
    //HW Lower triangular, so row i only touches the first i + 1 throws, from its skyline on
    for (int j = CholRowFirst[static_cast<size_t>(i)]; j <= i; ++j) Sum += L[i * NPars + j] * Throws[j];
    Velocity[i] = Sum;
  }
}

// *******************
void MetricBlock::Reflect(double* _restrict_ Velocity, const int i) const {
// *******************
  const size_t Row = static_cast<size_t>(i) * static_cast<size_t>(N);
  const double Variance = InvMass[Row + static_cast<size_t>(i)];
  if (Diagonal || !(Variance > 0.)) {
    Velocity[i] = -Velocity[i];
    return;
  }
  //HW M^-1 is symmetric, so its column i is row i, and the skyline bounds the non-zero part
  const double Scale = 2. * Velocity[i] / Variance;
  for (int j = InvRowFirst[static_cast<size_t>(i)]; j <= InvRowLast[static_cast<size_t>(i)]; ++j) {
    Velocity[j] -= Scale * InvMass[Row + static_cast<size_t>(j)];
  }
}

// *******************
void MetricBlock::Kick(double* _restrict_ Velocity, const double* _restrict_ Gradient, const double Step) const {
// *******************
  const double* _restrict_ Inverse = InvMass.data();
  const int NPars = N;
  if (Diagonal) {
    //HW The hot one: O(N) rather than O(N^2), twice per leapfrog step
    #ifdef MULTITHREAD
    #pragma omp parallel for simd
    #endif
    for (int i = 0; i < NPars; ++i) Velocity[i] -= 0.5 * Step * Inverse[i * NPars + i] * Gradient[i];
    return;
  }
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < NPars; ++i) {
    double Sum = 0.;
    #ifdef MULTITHREAD
    #pragma omp simd
    #endif
    for (int j = InvRowFirst[static_cast<size_t>(i)]; j <= InvRowLast[static_cast<size_t>(i)]; ++j) {
      Sum += Inverse[i * NPars + j] * Gradient[j];
    }
    Velocity[i] -= 0.5 * Step * Sum;
  }
}

// *******************
double MetricBlock::AddKineticEnergy(const double* _restrict_ Velocity, double* _restrict_ u, double ke) const {
// *******************
  //HW M^-1 = L L^T, so forward substituting L u = w gives 1/2 w^T M w as 1/2 |u|^2
  const double* _restrict_ L = Chol.data();
  const int NPars = N;
  if (Diagonal) {
    //HW No substitution to do, u[i] is just w[i] over its standard deviation
    for (int i = 0; i < NPars; ++i) {
      const double ui = Velocity[i] / L[i * NPars + i];
      ke += 0.5 * ui * ui;
    }
    return ke;
  }
  for (int i = 0; i < NPars; ++i) {
    double Sum = Velocity[i];
    for (int j = CholRowFirst[static_cast<size_t>(i)]; j < i; ++j) Sum -= L[i * NPars + j] * u[j];
    u[i] = Sum / L[i * NPars + i];
    ke += 0.5 * u[i] * u[i];
  }
  return ke;
}

// *******************
double MetricBlock::AddKineticEnergyAndMomentum(const double* _restrict_ Velocity, double* _restrict_ u,
                                                double* _restrict_ Momentum, double ke) const {
// *******************
  //HW M = L^-T L^-1. Forward substituting L u = w gives the kinetic energy as 1/2 |u|^2, and back
  //substituting L^T p = u then gives p = M w, both without forming M
  const double* _restrict_ L = Chol.data();
  const int NPars = N;
  if (Diagonal) {
    //HW M is the reciprocal of M^-1 element by element, so p = w / variance directly
    for (int i = 0; i < NPars; ++i) {
      const double Sd = L[i * NPars + i];
      const double ui = Velocity[i] / Sd;
      ke += 0.5 * ui * ui;
      Momentum[i] = ui / Sd;
    }
    return ke;
  }
  for (int i = 0; i < NPars; ++i) {
    double Sum = Velocity[i];
    for (int j = CholRowFirst[static_cast<size_t>(i)]; j < i; ++j) Sum -= L[i * NPars + j] * u[j];
    u[i] = Sum / L[i * NPars + i];
    ke += 0.5 * u[i] * u[i];
  }
  //HW L^T is upper triangular, so row i only touches the momenta already solved for
  for (int i = NPars - 1; i >= 0; --i) {
    double Sum = u[i];
    for (int j = i + 1; j <= CholColLast[static_cast<size_t>(i)]; ++j) Sum -= L[j * NPars + i] * Momentum[j];
    Momentum[i] = Sum / L[i * NPars + i];
  }
  return ke;
}
