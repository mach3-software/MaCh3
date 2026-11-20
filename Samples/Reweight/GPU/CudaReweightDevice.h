#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "Samples/Reweight/SweepTables.h"

/// @brief The device side of CudaReweight: memory and kernels for the cached pass and the gradient
///
/// @details Plain C++ with only standard types, so the CUDA code (CudaReweightDevice.cu) never has
/// to compile MaCh3's headers, and the rest of MaCh3 never has to see CUDA's. Everything is in
/// double precision and follows CpuReweight's arithmetic, so the two agree to rounding.
///
/// The inputs a pass needs every time (the pool's values, the derivatives) are written into pinned
/// host buffers the caller fills directly, then copied in one go.
///
/// @author Henry Wallace
class CudaReweightDevice {
 public:
  /// @brief Allocate and upload what never changes: the signatures
  /// @param SigStart,SigPool Signature s multiplies pool entries SigPool[SigStart[s]] onwards
  /// @param NPool Number of pool entries
  /// @param NBins Number of bins
  CudaReweightDevice(const std::vector<int>& SigStart, const std::vector<int>& SigPool, std::size_t NPool, int NBins);
  ~CudaReweightDevice();
  CudaReweightDevice(const CudaReweightDevice&) = delete;
  CudaReweightDevice& operator=(const CudaReweightDevice&) = delete;

  /// @brief Replace the groups (after a full pass rebuilt them)
  /// @param BinStart Bin b's groups are BinStart[b] onwards
  /// @param Sig,Const,Const2 Each group's signature and sums of constant weights and their squares
  /// @param SigGroupStart,SigGroup Signature s's groups are SigGroup[SigGroupStart[s]] onwards
  /// @param GroupBin Each group's bin
  void SetGroups(const std::vector<int>& BinStart, const std::vector<int>& Sig, const std::vector<double>& Const,
                 const std::vector<double>& Const2, const std::vector<int>& SigGroupStart, const std::vector<int>& SigGroup,
                 const std::vector<int>& GroupBin);

  /// @brief Upload the gradient's terms, grouped by parameter
  /// @param ParStart Parameter p's terms are ParStart[p] onwards
  /// @param TermSig,TermDeriv,TermPool Each term's signature, derivative and pool entry
  /// @param NDeriv Number of derivatives the terms index
  void SetTerms(const std::vector<int>& ParStart, const std::vector<int>& TermSig, const std::vector<int>& TermDeriv,
                const std::vector<int>& TermPool, std::size_t NDeriv);

  /// @brief Where to write the pool's values before Fill(), NPool of them (pinned host memory)
  double* PoolBuffer();
  /// @brief Where to write the derivatives before Gradient(), NDeriv of them (pinned host memory)
  double* DerivBuffer();

  /// @brief A cached pass: multiply out the signatures and sum each bin's groups
  /// @param MC,W2 NBins entries each, filled
  void Fill(double* MC, double* W2);
  /// @brief The gradient at the last Fill()
  /// @param BinFactor NBins per bin factors
  /// @param Grad Filled with one entry per parameter
  void Gradient(const std::vector<double>& BinFactor, std::vector<double>& Grad);

  /// @brief Upload what a full pass on the device needs, see SweepTables
  /// @param Tables The live events' variables, cuts and binning
  /// @param Sig,Const Each live event's signature and constant weight, as ReweightModel has them
  void SetSweep(const SweepTables& Tables, const std::vector<int>& Sig, const std::vector<double>& Const);
  /// @brief Whether SetSweep was called
  bool HasSweep() const;
  /// @brief A full pass: cut and bin every live event at the given shifts, and sum the weights
  /// @param Delta The shifts' offsets from the table's base point, SweepTables::NShifts of them
  /// @param MC,W2 NBins entries each, filled
  /// @details Reads the pool from PoolBuffer() like Fill(). Each live event's bin stays on the
  /// device (GetEventBins), and the bins are summed in a fixed order, so the result is the same
  /// from run to run
  void FullPass(const std::vector<double>& Delta, double* MC, double* W2);
  /// @brief Each live event's bin from the last FullPass, -1 for one that isn't selected or binned
  void GetEventBins(int* Out);

  /// @brief The device's name
  static std::string DeviceName();
  /// @brief Whether a CUDA device is available
  static bool Available();
  /// @brief Device memory held, in bytes
  std::size_t Bytes() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> D;
};
