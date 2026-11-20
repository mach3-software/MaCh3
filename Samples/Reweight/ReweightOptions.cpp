#include "Samples/Reweight/ReweightOptions.h"
#include "Samples/Reweight/CpuReweight.h"
#ifdef MaCh3_METAL
#include "Samples/Reweight/GPU/MetalReweight.h"
#endif
#ifdef MaCh3_CUDA
#include "Samples/Reweight/GPU/CudaReweight.h"
#endif

// ************************************************
ReweightOptions ReweightOptions::Read(const YAML::Node& Config, const std::string& SampleName) {
// ************************************************
  // Each key is optional; the struct's own values are the defaults
  ReweightOptions O;
  if (CheckNodeExists(Config, "Reweight", "Backend")) O.Backend = Get<std::string>(Config["Reweight"]["Backend"], __FILE__ , __LINE__);
  if (CheckNodeExists(Config, "Reweight", "SelfCheckPoints")) O.SelfCheckPoints = Get<int>(Config["Reweight"]["SelfCheckPoints"], __FILE__ , __LINE__);
  if (CheckNodeExists(Config, "Reweight", "FastAccess")) O.FastAccess = Get<bool>(Config["Reweight"]["FastAccess"], __FILE__ , __LINE__);
  if (CheckNodeExists(Config, "Reweight", "LinearShifts")) O.LinearShifts = Get<bool>(Config["Reweight"]["LinearShifts"], __FILE__ , __LINE__);

  // Refuse a backend that doesn't exist, or wasn't built into this MaCh3
  if (O.Backend != "auto" && O.Backend != "legacy" && O.Backend != "cpu" && O.Backend != "metal" && O.Backend != "cuda") {
    MACH3LOG_ERROR("{}: Reweight.Backend is '{}', but has to be one of auto, legacy, cpu, cuda or metal", SampleName, O.Backend);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  #ifndef MaCh3_CUDA
  if (O.Backend == "cuda") {
    MACH3LOG_ERROR("{}: Reweight.Backend is cuda, but MaCh3 was built without CUDA", SampleName);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  #endif
  #ifndef MaCh3_METAL
  if (O.Backend == "metal") {
    MACH3LOG_ERROR("{}: Reweight.Backend is metal, but MaCh3 was built without it (MaCh3_METAL_ENABLED, Apple only)", SampleName);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  #endif
  return O;
}

// ************************************************
std::unique_ptr<ReweightBackend> MakeReweightBackend([[maybe_unused]] const std::string& Choice, ReweightModel Model, ReweightEventAccess& Access) {
// ************************************************
  #ifdef MaCh3_METAL
  if (Choice == "metal") return std::make_unique<MetalReweight>(std::move(Model), Access);
  #endif
  #ifdef MaCh3_CUDA
  //HW auto takes the GPU when there is one: it agrees with the CPU to rounding
  if (Choice == "cuda" || (Choice == "auto" && CudaReweight::Available())) {
    return std::make_unique<CudaReweight>(std::move(Model), Access);
  }
  #endif
  // cpu, or auto without a GPU
  return std::make_unique<CpuReweight>(std::move(Model), Access);
}
