#include "Samples/Reweight/GPU/CudaReweight.h"

#include "Samples/Reweight/GPU/CudaReweightDevice.h"
#include "Samples/Reweight/FastEventAccess.h"

// ************************************************
CudaReweight::CudaReweight(ReweightModel Model_, ReweightEventAccess& Access_)
    : CpuReweight(std::move(Model_), Access_),
      GPU(std::make_unique<CudaReweightDevice>(Model.SigStart, Model.SigPool, Model.PoolPointer.size(), Model.NBins)) {
// ************************************************
  DeviceMC.assign(static_cast<std::size_t>(Model.NBins), 0.);
  DeviceW2.assign(static_cast<std::size_t>(Model.NBins), 0.);
  // Full passes on the GPU, if the shifts' effect is known in closed form
  if (const auto* Fast = dynamic_cast<const FastEventAccess*>(&Access)) {
    SweepTables Tables;
    if (Fast->ExportSweepTables(Model.Event, Tables)) {
      GPU->SetSweep(Tables, Model.Sig, Model.Const);
      Sweep = Fast;
    }
  }
}

// ************************************************
CudaReweight::~CudaReweight() = default;
// ************************************************

// ************************************************
bool CudaReweight::Available() {
// ************************************************
  return CudaReweightDevice::Available();
}

// ************************************************
std::string CudaReweight::Summary() const {
// ************************************************
  return CpuReweight::Summary() + (Sweep ? "; full passes, cached passes and gradients on " : "; cached passes and gradients on ") +
         CudaReweightDevice::DeviceName() + (Sweep ? "" : " (full passes on the CPU: they need linear shifts for every live event)");
}

// ************************************************
std::size_t CudaReweight::Bytes() const {
// ************************************************
  return CpuReweight::Bytes() + GPU->Bytes();
}

// ************************************************
void CudaReweight::ReadPool() {
// ************************************************
  double* Pool = GPU->PoolBuffer();
  const int NPool = static_cast<int>(Model.PoolPointer.size());
  const M3::float_t* const* Pointer = Model.PoolPointer.data();
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int k = 0; k < NPool; ++k) Pool[k] = static_cast<double>(*Pointer[k]);
}

// ************************************************
void CudaReweight::SyncLiveBin() {
// ************************************************
  if (!BinsOnDevice) return;
  GPU->GetEventBins(LiveBin.data());
  BinsOnDevice = false;
}

// ************************************************
void CudaReweight::Fill(const std::vector<double>& ShiftKey, const bool FillW2, double* MC, double* W2) {
// ************************************************
  if (!BinsCurrent(ShiftKey)) {
    if (Sweep == nullptr) {
      // The shifts moved: a full pass, on the CPU, through the experiment's own code
      LastFillOnGPU = false;
      CpuReweight::Fill(ShiftKey, FillW2, MC, W2);
      return;
    }
    // The shifts moved: a full pass on the GPU from the linear shift table
    Sweep->ShiftDeltas(Delta);
    ReadPool();
    GPU->FullPass(Delta, DeviceMC.data(), DeviceW2.data());
    BinKey = ShiftKey;
    BinsValid = true;
    GroupsValid = false;
    BinsOnDevice = true;
  } else {
    EnsureGroups();
    if (UploadedGroups != GroupsVersion) {
      std::vector<int> GroupBin, SigGroupStart, SigGroup;
      IndexGroups(GroupBin, SigGroupStart, SigGroup);
      GPU->SetGroups(Groups.BinStart, Groups.Sig, Groups.Const, Groups.Const2, SigGroupStart, SigGroup, GroupBin);
      UploadedGroups = GroupsVersion;
    }
    ReadPool();
    GPU->Fill(DeviceMC.data(), DeviceW2.data());
  }
  for (int b = 0; b < Model.NBins; ++b) {
    MC[b] += DeviceMC[static_cast<std::size_t>(b)];
    if (FillW2) W2[b] += DeviceW2[static_cast<std::size_t>(b)];
  }
  LastFillOnGPU = true;
}

// ************************************************
void CudaReweight::AttachDerivatives(const std::vector<WeightDerivativeSource>& Sources, const int NPars_) {
// ************************************************
  CpuReweight::AttachDerivatives(Sources, NPars_);
  std::vector<int> ParStart, TermSig, TermDeriv, TermPool;
  ParameterTerms(ParStart, TermSig, TermDeriv, TermPool);
  GPU->SetTerms(ParStart, TermSig, TermDeriv, TermPool, PoolDerivative.size());
}

// ************************************************
void CudaReweight::Gradient(const std::vector<double>& BinFactor, std::vector<double>& Out) {
// ************************************************
  if (!LastFillOnGPU) { CpuReweight::Gradient(BinFactor, Out); return; }
  // After a full pass on the GPU its gradient needs the groups of the new bins, which the GPU's
  // gradient kernels read; the signature products are already there
  EnsureGroups();
  if (UploadedGroups != GroupsVersion) {
    std::vector<int> GroupBin, SigGroupStart, SigGroup;
    IndexGroups(GroupBin, SigGroupStart, SigGroup);
    GPU->SetGroups(Groups.BinStart, Groups.Sig, Groups.Const, Groups.Const2, SigGroupStart, SigGroup, GroupBin);
    UploadedGroups = GroupsVersion;
  }
  double* Deriv = GPU->DerivBuffer();
  const int NDeriv = static_cast<int>(PoolDerivative.size());
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int j = 0; j < NDeriv; ++j) {
    const M3::float_t* D = PoolDerivative[static_cast<std::size_t>(j)];
    Deriv[j] = (D == nullptr) ? 1. : static_cast<double>(*D);
  }
  GPU->Gradient(BinFactor, Out);
}
