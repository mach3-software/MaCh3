#pragma once

#include <cstddef>
#include <memory>
#include <vector>

/// @file DeviceReweight.h
/// @brief SampleHandlerBase's reweight on the GPU: both the cached, grouped pass a leapfrog makes
/// (the functional, shift, parameters sit still and only the weights move) and the full pass of a
/// random walk sweep (the shifts move, so every event is shifted, cut and binned again).
/// @details Plain C++ on purpose, so SampleHandlerBase.cpp needn't be compiled by nvcc. The
/// implementation, DeviceReweight.cu, only exists in a CUDA build. Everything is in double
/// precision and follows the CPU code's arithmetic, so the two agree to rounding error.
///
/// Every varying weight a signature multiplies, and every derivative, is a value in one device
/// array, the pool:
///   [0, NPars)                        parameter values, in ParameterHandler order
///   [NPars, NPars+NSplines)           spline weights, evaluated on the device each reweight
///   [NPars+NSplines, NPars+2NSplines) d(spline weight)/d(its parameter)
///   then NGathered values the host gathers each reweight (oscillation weights and anything else)
///   then the two constants 1 and 0
/// so a signature is a list of pool indices, which SampleHandlerBase builds once from its pointers.
///
/// A leapfrog reweight (Run): evaluate the splines (the cubic of BinnedSplineHandler::CalcSplineWeights),
/// weigh every group, sum each bin's groups into the MC prediction (and W2), then accumulate
/// d(-lnL)/d(parameter) straight into NPars numbers with the Poisson pre-factor 1 - data/MC that
/// SampleHandlerBase::GetGradLikelihood uses. That avoids the NBins x NPars gradient array the CPU
/// path fills, and the sums are regrouped by signature and then by parameter, so the parameters
/// every event shares don't turn into atomic contention.
///
/// A full pass (RunFull): every live event's shifted variables from per event coefficients (the
/// shifts have to be additive and linear, which SampleHandlerBase checks before using this), its
/// cuts, its bin and its weight. The event bins stay on the device, and the groups for the
/// leapfrogs that follow are built there too (BuildGroupsFromFullPass), so nothing per event
/// crosses to the host.
///
/// Each reweight sends the parameter values (and spline segments, gathered weights, shift values)
/// and brings back NBins, 2 NBins or 2 NBins + NPars doubles, on one stream with one wait.
/// @author Henry Wallace

/// @brief Everything the device needs that only changes when the model does, uploaded once
struct DeviceReweightStatic {
  // ----- Sizes and data
  /// Number of bins over all samples of the handler
  int NBins = 0;
  /// Data in each bin
  std::vector<double> Data;
  /// Number of parameters, which is also the length of the gradient
  int NPars = 0;
  /// Below this MC prediction a bin is clamped and contributes no gradient (M3::_LOW_MC_BOUND_)
  double LowMCBound = 0.;

  // ----- Splines, as BinnedSplineHandler lays them out (all empty and NSplines = 0 without splines)
  /// Number of splines
  int NSplines = 0;
  /// Number of spline parameters
  int NSplineParams = 0;
  /// Most knots any spline parameter has
  int MaxKnots = 0;
  /// Knot positions, MaxKnots per spline parameter
  std::vector<double> KnotX;
  /// The y, b, c, d coefficients of every spline segment, four per segment
  std::vector<double> SegmentCoeffs;
  /// Where each spline's segments start in SegmentCoeffs (counted in segments)
  std::vector<unsigned int> SplineCoeffStart;
  /// Which spline parameter each spline belongs to
  std::vector<short int> SplineParam;
  /// Index of each spline parameter in the parameter handler, so its value can be read from the pool
  std::vector<int> SplineParamIndex;

  // ----- Values the host gathers into the pool each reweight (oscillation weights and the like)
  int NGathered = 0;

  // ----- Signatures: each is a list of pool indices, flattened. Signature s owns entries Start[s] to Start[s+1]-1
  /// Its varying weights
  std::vector<int> SigWeightStart;
  /// See SigWeightStart
  std::vector<int> SigWeightIdx;
  /// Its weights that have a derivative: the parameter, the derivative's pool index and the weight's pool index
  std::vector<int> SigDerivStart;
  /// See SigDerivStart
  std::vector<int> SigDerivPar;
  /// See SigDerivStart
  std::vector<int> SigDerivIdx;
  /// See SigDerivStart
  std::vector<int> SigDerivWeightIdx;
};

/// @brief Groups of events built on the host (after a CPU full pass), sorted by bin
struct DeviceReweightGroups {
  /// Signature of each group
  std::vector<int> Sig;
  /// Sum over the group's events of their constant weights
  std::vector<double> Const;
  /// Sum of the squared constant weights, for W2
  std::vector<double> Const2;
  /// Bin b's groups are BinStart[b] to BinStart[b+1]-1, NBins + 1 entries
  std::vector<int> BinStart;
};

/// @brief What a full event pass needs on the device, for the reweights where the functional
/// (shift) parameters move: the random walk sweeps
/// @details Only for handlers whose shifts are additive and linear in their parameter: every shifted
/// variable that an event's cuts and binning read is then its nominal value plus, for each of its
/// shifts in turn, shift value times a fixed coefficient. SampleHandlerBase works the coefficients out
/// by applying each shift at 1 to zeroed variables, and checks them against its own ApplyShifts at
/// random values before using them. Only live events are included (those not ruled out once and for
/// all by a cut on a variable no shift changes)
struct DeviceSweepStatic {
  // ----- Live events
  int NLive = 0;
  /// Sample of each live event
  std::vector<int> Sample;
  /// Signature of each live event
  std::vector<int> Sig;
  /// Product of each live event's constant weights
  std::vector<double> ConstWeight;

  // ----- Shifted variables: variable 0 is the binning variable, variable 1 + j the j-th shifted cut of the sample
  /// Variables per event (the same for every event, unused ones are left at zero)
  int NVariables = 1;
  /// Nominal value of variable v of live event i at i * NVariables + v
  std::vector<double> Nominal;
  /// Coefficients of variable v of live event i: entries CoefStart[i*NVariables+v] to CoefStart[i*NVariables+v+1]-1,
  /// in the order the CPU applies the shifts
  std::vector<int> CoefStart;
  /// Which shift each coefficient multiplies
  std::vector<int> CoefShift;
  /// The coefficients
  std::vector<double> Coef;
  /// Number of shift parameters
  int NShifts = 0;

  // ----- Cuts on shifted variables, per sample: sample s's cuts are CutStart[s] to CutStart[s+1]-1
  std::vector<int> CutStart;
  /// Which variable each cut reads
  std::vector<int> CutVariable;
  /// Accepted range of each cut, [CutLo, CutHi)
  std::vector<double> CutLo;
  /// See CutLo
  std::vector<double> CutHi;

  // ----- Binning, per sample: uniform 1D, sample s's edges are EdgeStart[s] to EdgeStart[s+1]-1
  std::vector<int> EdgeStart;
  std::vector<double> Edges;
  /// Global bin = bin within the sample * Stride + Offset
  std::vector<int> Stride;
  /// See Stride
  std::vector<int> Offset;

  // ----- The codes the CPU stores for an event outside the binning and for one a cut rejects
  int UnderOverFlowBin = -1;
  int NotSelected = -2;
};

class DeviceReweightImpl;

/// @brief Owns the device memory and runs the reweight, see the file description
class DeviceReweight {
 public:
  /// @brief Allocate and upload everything in Static
  explicit DeviceReweight(const DeviceReweightStatic& Static);
  ~DeviceReweight();
  DeviceReweight(const DeviceReweight&) = delete;
  DeviceReweight& operator=(const DeviceReweight&) = delete;

  /// @brief Replace the groups, and the data (which is cheap, and may have been reset since)
  void SetGroups(const DeviceReweightGroups& FromHost, const std::vector<double>& Data);
  /// @brief Replace the data, NBins values
  void SetData(const std::vector<double>& Data);

  /// @brief One reweight
  /// @param ParValues NPars current parameter values
  /// @param Segments NSplineParams spline segments, as SplineBase::FindSplineSegment finds them
  /// @param Gathered NGathered values the host gathered
  /// @param FillW2 Also fill the sum of squared weights per bin
  /// @param Gradient Also compute d(-lnL)/d(parameter) into GradOut
  /// @param MCOut NBins predicted events per bin
  /// @param W2Out NBins sums of squared weights, written only with FillW2
  /// @param GradOut NPars gradient entries, written only with Gradient
  void Run(const double* ParValues, const short int* Segments, const double* Gathered, bool FillW2, bool Gradient,
           double* MCOut, double* W2Out, double* GradOut);

  /// @brief Total device memory allocated, in bytes
  std::size_t DeviceBytes() const;

  /// @brief Upload what full passes need, see DeviceSweepStatic
  void SetSweep(const DeviceSweepStatic& Sweep);
  /// @brief One full pass: shifts, cuts, binning and weights of every live event
  /// @param ShiftValues NShifts current shift parameter values
  /// @details Arguments otherwise as Run. Each live event's bin stays on the device for GetEventBins
  void RunFull(const double* ParValues, const short int* Segments, const double* Gathered, const double* ShiftValues,
               bool FillW2, bool Gradient, double* MCOut, double* W2Out, double* GradOut);
  /// @brief The bin (or code) of every live event from the last RunFull, NLive entries
  void GetEventBins(int* Out) const;
  /// @brief Build the groups on the device from the last RunFull's event bins, so the cached passes
  /// that follow need nothing from the host. Same groups, in the same order, as SampleHandlerBase's
  /// BuildCachedEventList and BuildEventGroups make from the same bins. Returns how many
  int BuildGroupsFromFullPass();

 private:
  std::unique_ptr<DeviceReweightImpl> Impl;
};
