#pragma once

#include <cstddef>
#include <vector>

/// @brief The Gaussian prior \f$\frac{1}{2}\Delta^T M^{-1}\Delta\f$ and its gradient, summed over
/// only the entries of \f$M^{-1}\f$ that can contribute
///
/// @details A parameter with a flat prior drops out of the penalty, and so does every entry of the
/// inverse covariance that is exactly zero. With hundreds of flat prior parameters, or a block
/// diagonal covariance, that leaves a small fraction of the matrix. This keeps those entries in
/// compressed rows, in ascending column order, which is the order the dense loops visit them in.
/// Every sum therefore adds the same numbers in the same order as the dense version, less terms
/// that are exactly zero. On one thread that gives the same result bit for bit. It does not under
/// MULTITHREAD: LogLikelihood() is an OpenMP reduction, whose cross-thread combining order is
/// implementation defined, and a 16 thread run on Rorqual (19 Nov) was measured varying in the last
/// digits between two runs of the same build.
///
/// The object remembers which parameters were flat when it was built, so Matches() can tell a
/// caller that the flags have changed underneath it.
///
/// @author Henry Wallace
class SparsePrior {
 public:
  /// @brief Collect the contributing entries of an inverse covariance
  /// @param InverseCovariance \f$M^{-1}\f$, square, one row per parameter
  /// @param FlatPrior Per parameter, whether it has a flat prior and so no penalty
  void Build(const std::vector<std::vector<double>>& InverseCovariance, const std::vector<bool>& FlatPrior) {
    const size_t NPars = FlatPrior.size();
    BuiltWithFlatPrior = FlatPrior;
    RowOf.assign(NPars, kNoRow);
    RowStart.assign(1, 0);
    LowerCount.clear();
    Column.clear();
    Value.clear();
    for (size_t i = 0; i < NPars; ++i) {
      if (FlatPrior[i]) continue;
      RowOf[i] = static_cast<int>(LowerCount.size());
      int Lower = 0;
      for (size_t j = 0; j < NPars; ++j) {
        if (FlatPrior[j] || InverseCovariance[i][j] == 0.) continue;
        Column.push_back(static_cast<int>(j));
        Value.push_back(InverseCovariance[i][j]);
        if (j <= i) ++Lower;
      }
      LowerCount.push_back(Lower);
      RowStart.push_back(static_cast<int>(Column.size()));
    }
    Built = true;
  }

  /// @brief Whether Build() has run for these flat prior flags
  bool Matches(const std::vector<bool>& FlatPrior) const { return Built && BuiltWithFlatPrior == FlatPrior; }

  /// @brief Forget the entries, e.g. because the covariance changed
  void Invalidate() { Built = false; }

  /// @brief Number of matrix entries kept
  size_t GetNEntries() const { return Value.size(); }

  /// @brief \f$\frac{1}{2}\Delta^T M^{-1}\Delta\f$ with \f$\Delta\f$ = Values - Central
  /// @details Sums the lower triangle and doubles the off diagonal terms, as the dense loop does
  template <typename ValueVector, typename CentralVector>
  double LogLikelihood(const ValueVector& Values, const CentralVector& Central) const {
    const int NPars = static_cast<int>(RowOf.size());
    double LogL = 0.;
    #ifdef MULTITHREAD
    #pragma omp parallel for reduction(+:LogL)
    #endif
    for (int i = 0; i < NPars; ++i) {
      const int Row = RowOf[static_cast<size_t>(i)];
      if (Row == kNoRow) continue;
      const double Diff = Values[static_cast<size_t>(i)] - Central[static_cast<size_t>(i)];
      const int Begin = RowStart[static_cast<size_t>(Row)];
      const int End = Begin + LowerCount[static_cast<size_t>(Row)];
      for (int k = Begin; k < End; ++k) {
        const size_t j = static_cast<size_t>(Column[static_cast<size_t>(k)]);
        const double Scale = (static_cast<int>(j) != i) ? 1. : 0.5;
        LogL += Scale * Diff * (Values[j] - Central[j]) * Value[static_cast<size_t>(k)];
      }
    }
    return LogL;
  }

  /// @brief The gradient of LogLikelihood(), \f$M^{-1}\Delta\f$, zero for flat prior parameters
  /// @param Out Resized to one entry per parameter and overwritten; pass the same buffer each time
  /// to keep a sampler's gradient evaluations from allocating, see ParameterHandlerBase::CalcGradient
  template <typename ValueVector, typename CentralVector>
  void Gradient(const ValueVector& Values, const CentralVector& Central, std::vector<double>& Out) const {
    const int NPars = static_cast<int>(RowOf.size());
    Out.assign(static_cast<size_t>(NPars), 0.);
    double* Grad = Out.data();
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < NPars; ++i) {
      const int Row = RowOf[static_cast<size_t>(i)];
      if (Row == kNoRow) continue;
      double Sum = 0.;
      for (int k = RowStart[static_cast<size_t>(Row)]; k < RowStart[static_cast<size_t>(Row) + 1]; ++k) {
        const size_t j = static_cast<size_t>(Column[static_cast<size_t>(k)]);
        Sum += Value[static_cast<size_t>(k)] * (Values[j] - Central[j]);
      }
      Grad[static_cast<size_t>(i)] = Sum;
    }
  }

  /// @brief Gradient() into a fresh vector, for callers that aren't on a hot path
  template <typename ValueVector, typename CentralVector>
  std::vector<double> Gradient(const ValueVector& Values, const CentralVector& Central) const {
    std::vector<double> Out;
    Gradient(Values, Central, Out);
    return Out;
  }

 private:
  /// RowOf entry of a parameter with a flat prior
  static constexpr int kNoRow = -1;
  /// Whether the lists below are filled
  bool Built = false;
  /// The flat prior flags the lists were built with
  std::vector<bool> BuiltWithFlatPrior;
  /// Compressed row of each parameter, kNoRow for a flat prior
  std::vector<int> RowOf;
  /// Row r's entries are RowStart[r] to RowStart[r+1]-1
  std::vector<int> RowStart;
  /// How many of a row's entries have column <= row: the lower triangle, which comes first
  std::vector<int> LowerCount;
  /// Column of each entry, ascending within a row
  std::vector<int> Column;
  /// Value of each entry
  std::vector<double> Value;
};
