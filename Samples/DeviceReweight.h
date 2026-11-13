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
///   then NOther values the host gathers each reweight (oscillation weights and anything else)
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
  /// Number of bins over all samples of the handler, and the data in each
  int NBins = 0;
  std::vector<double> Data;
  /// Number of parameters, the length of the gradient
  int NPars = 0;
  /// Splines, as BinnedSplineHandler lays them out. NSplines = 0 if there are none
  int NSplines = 0;
  int NSplineParams = 0;
  int MaxKnots = 0;
  /// x of every knot of every spline parameter, NSplineParams * MaxKnots
  std::vector<double> XCoeff;
  /// y, b, c, d of every spline segment
  std::vector<double> ManyCoeff;
  /// Where each spline's coefficients start in ManyCoeff, in segments, and which spline parameter it is
  std::vector<unsigned int> NKnots;
  std::vector<short int> ParamNo;
  /// Global parameter index of each spline parameter, so its value can be read from the pool
  std::vector<int> SplineGlobalPar;
  /// Number of host gathered values in the pool
  int NOther = 0;
  /// Signatures: varying weight pool indices, flat, SigWStart[s] to SigWStart[s+1]-1
  std::vector<int> SigWStart, SigWIdx;
  /// Signatures: for each weight with a derivative, the parameter, the derivative's pool index and
  /// the weight's pool index, flat, SigDStart[s] to SigDStart[s+1]-1
  std::vector<int> SigDStart, SigDPar, SigDDeriv, SigDWeight;
  /// Below this the prediction is clamped and the bin contributes no gradient
  double LowMCBound = 0.;
};

/// @brief Groups of events, which change whenever the shifts do (once per NUTS trajectory)
struct DeviceReweightGroups {
  std::vector<int> Sig;
  std::vector<double> Const, Const2;
  /// Groups are sorted by bin: bin b's groups are BinStart[b] to BinStart[b+1]-1, NBins + 1 entries
  std::vector<int> BinStart;
};

/// @brief What a full event pass needs on the device, for the reweights where the functional
/// (shift) parameters move: the random walk sweeps
/// @details Only for handlers whose shifts are additive and linear in their parameter, so every
/// shifted variable an event's cuts and binning read is its nominal value plus, for each of its
/// shifts in order, parameter value times a fixed coefficient. SampleHandlerBase works the
/// coefficients out by applying each shift at 1 to zeroed variables, and checks the result against
/// its own ApplyShifts at random values before using them. Only live events are here (those not
/// ruled out once and for all by a cut on an unshifted variable)
struct DeviceSweepStatic {
  /// Live events: their index in the handler, sample, signature and product of constant weights
  int NLive = 0;
  std::vector<int> EventIndex, Sample, Sig;
  std::vector<double> ConstW;
  /// Variables per event: slot 0 is the binning variable, slot 1 + j the j-th shifted cut of its sample
  int NSlots = 1;
  /// NLive * NSlots nominal values
  std::vector<double> Nominal;
  /// Coefficients of slot j of live event i: CoefStart[i*NSlots+j] to CoefStart[i*NSlots+j+1]-1,
  /// in the order the CPU applies the shifts
  std::vector<int> CoefStart, CoefShift;
  std::vector<double> Coef;
  int NShifts = 0;
  /// Per sample: its shifted cuts, CutStart[s] to CutStart[s+1]-1, as the slot each reads and
  /// the accepted range [Lo, Hi)
  std::vector<int> CutStart, CutSlot;
  std::vector<double> CutLo, CutHi;
  /// Per sample: uniform 1D binning, edges EdgeStart[s] to EdgeStart[s+1]-1, global bin = bin * Stride + Offset
  std::vector<int> EdgeStart, Stride, Offset;
  std::vector<double> Edges;
  /// Codes the CPU uses: outside the binning, and rejected by a cut
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
  void SetGroups(const DeviceReweightGroups& Groups, const std::vector<double>& Data);
  /// @brief Replace the data, NBins values
  void SetData(const std::vector<double>& Data);

  /// @brief One reweight
  /// @param ParValues NPars current parameter values
  /// @param Segments NSplineParams spline segments, as SplineBase::FindSplineSegment finds them
  /// @param Other NOther host gathered values
  /// @param FillW2 Also fill the sum of squared weights per bin
  /// @param Gradient Also compute d(-lnL)/d(parameter) into GradOut
  /// @param MCOut NBins predicted events per bin
  /// @param W2Out NBins sums of squared weights, written only with FillW2
  /// @param GradOut NPars gradient entries, written only with Gradient
  void Run(const double* ParValues, const short int* Segments, const double* Other, bool FillW2, bool Gradient,
           double* MCOut, double* W2Out, double* GradOut);

  /// @brief Total device memory allocated, in bytes
  std::size_t DeviceBytes() const;

  /// @brief Upload what full passes need, see DeviceSweepStatic
  void SetSweep(const DeviceSweepStatic& Sweep);
  /// @brief One full pass: shifts, cuts, binning and weights of every live event
  /// @param ShiftValues NShifts current shift parameter values
  /// @details Arguments otherwise as Run. Each live event's bin stays on the device for GetEventBins
  void RunFull(const double* ParValues, const short int* Segments, const double* Other, const double* ShiftValues,
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
