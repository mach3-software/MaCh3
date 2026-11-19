#include "Algorithms/HybridSettings.h"
#include "Manager/MaCh3Logger.h"

#include <algorithm>

namespace {
  // ***************
  /// @brief Throw with an error if a setting isn't strictly between 0 and 1
  void RequireOpenUnitInterval(const double Value, const std::string& Key, const std::string& Why = "") {
  // ***************
    if (Value > 0. && Value < 1.) return;
    MACH3LOG_ERROR("{} has to be strictly between 0 and 1, you asked for {}{}", Key, Value, Why);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // ***************
  /// @brief The systematic types named in RandomWalkTypes, throwing on a name that isn't one
  std::vector<SystType> ParseSystTypes(const std::vector<std::string>& TypeNames) {
  // ***************
    std::vector<SystType> Types;
    for (const auto& Name : TypeNames) {
      bool Found = false;
      for (int t = 0; t < SystType::kSystTypes; ++t) {
        const auto Type = static_cast<SystType>(t);
        if (SystType_ToString(Type) != Name) continue;
        Types.push_back(Type);
        Found = true;
        break;
      }
      if (!Found) {
        MACH3LOG_ERROR("RandomWalkTypes has '{}', which isn't a systematic type. Valid ones are:", Name);
        for (int t = 0; t < SystType::kSystTypes; ++t) {
          MACH3LOG_ERROR("  {}", SystType_ToString(static_cast<SystType>(t)));
        }
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    return Types;
  }
}

// ***************
HybridSettings HybridSettings::Read(const YAML::Node& Hybrid, const bool AutoMode, const unsigned int WarmupSteps) {
// ***************
  HybridSettings S;

  const auto Scheme = GetFromManager<std::string>(Hybrid["Scheme"], "composition", __FILE__, __LINE__);
  if (Scheme == "composition") S.Compose = true;
  else if (Scheme == "mixture") S.Compose = false;
  else {
    MACH3LOG_ERROR("Hybrid Scheme is '{}', which isn't a scheme. Use 'composition' to run both "
                   "kernels every step, or 'mixture' to toss a coin for one of them", Scheme);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  S.RandomWalkProb = GetFromManager<double>(Hybrid["RandomWalkProb"], 0.5, __FILE__, __LINE__);
  S.RandomWalkProbGiven = static_cast<bool>(Hybrid["RandomWalkProb"]);

  //HW "auto" (and Auto's default) chooses the number of sweeps during warmup, see HybridMCMC::TuneSweepCount.
  //It starts from 10
  S.AutoSweeps = Hybrid["RandomWalkSweeps"] ? (Hybrid["RandomWalkSweeps"].as<std::string>() == "auto") : AutoMode;
  S.Sweeps = S.AutoSweeps ? 10 : GetFromManager<int>(Hybrid["RandomWalkSweeps"], 1, __FILE__, __LINE__);
  S.TimeShare = GetFromManager<double>(Hybrid["RandomWalkTimeShare"], 0.4, __FILE__, __LINE__);
  RequireOpenUnitInterval(S.TimeShare, "RandomWalkTimeShare");

  S.Names = GetFromManager<std::vector<std::string>>(Hybrid["RandomWalkNames"], {}, __FILE__, __LINE__);
  S.FlipMoves = GetFromManager<int>(Hybrid["FlipMoves"], -1, __FILE__, __LINE__);

  //HW Auto scales the tuning interval to the warmup, so a short warmup still gets tuned
  const int DefaultCheckEvery = AutoMode ? std::max(10, std::min(250, static_cast<int>(WarmupSteps) / 10)) : 1000;
  S.CheckEvery = GetFromManager<int>(Hybrid["RandomWalkCheckEvery"], DefaultCheckEvery, __FILE__, __LINE__);
  S.SurrogateStride = GetFromManager<int>(Hybrid["RandomWalkSurrogateStride"], 0, __FILE__, __LINE__);
  if (S.SurrogateStride < 0) {
    MACH3LOG_ERROR("RandomWalkSurrogateStride has to be 0 (off) or positive, you asked for {}", S.SurrogateStride);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  S.ScaleFactor = GetFromManager<double>(Hybrid["RandomWalkScale"], 1.0, __FILE__, __LINE__);
  S.TuneScale = GetFromManager<bool>(Hybrid["RandomWalkTuneScale"], true, __FILE__, __LINE__);
  S.TargetAcceptance = GetFromManager<double>(Hybrid["RandomWalkTargetAcceptance"], 0.234, __FILE__, __LINE__);
  S.PenaltyExponent = GetFromManager<double>(Hybrid["RandomWalkPenaltyExponent"], 1.0, __FILE__, __LINE__);
  if (S.PenaltyExponent <= 0.) {
    MACH3LOG_ERROR("RandomWalkPenaltyExponent has to be positive, you asked for {}. It is the "
                   "p in dlogL proportional to scale^p: 1 for a V shaped log posterior, 2 for "
                   "a smooth quadratic one", S.PenaltyExponent);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  RequireOpenUnitInterval(S.TargetAcceptance, "RandomWalkTargetAcceptance");
  if (S.TuneScale && S.CheckEvery <= 0) {
    MACH3LOG_ERROR("RandomWalkTuneScale is on but RandomWalkCheckEvery is {}, so there is no "
                   "batch to measure an acceptance over. Give it a positive number of steps, "
                   "or turn the tuning off", S.CheckEvery);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (S.ScaleFactor <= 0.) {
    MACH3LOG_ERROR("RandomWalkScale has to be positive, you asked for {}", S.ScaleFactor);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Functional parameters are the ones the hybrid sampler exists for. Oscillation parameters are in
  // by default for the same reason: NuOscillator hands back no derivative, so SampleHandlerBase
  // drops their gradient and HMC would sample them from their prior
  //HW Without RandomWalkTypes the blocks are chosen automatically: every parameter some sample
  //can't differentiate with respect to goes to the random walk, see HybridBlocks
  S.AutoTypes = !Hybrid["RandomWalkTypes"];
  S.Types = ParseSystTypes(GetFromManager<std::vector<std::string>>(Hybrid["RandomWalkTypes"], {}, __FILE__, __LINE__));

  if (S.Compose) {
    if (S.Sweeps < 1) {
      MACH3LOG_ERROR("RandomWalkSweeps has to be at least 1, you asked for {}. At 0 the "
                     "gradient free parameters would never move", S.Sweeps);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  } else {
    RequireOpenUnitInterval(S.RandomWalkProb, "RandomWalkProb",
                            ". At 0 the gradient free parameters never move and at 1 nothing else does");
    if (S.Sweeps != 1) {
      MACH3LOG_ERROR("RandomWalkSweeps is {} but the scheme is mixture, which runs one kernel "
                     "per step. Set Scheme to composition to use sweeps", S.Sweeps);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  return S;
}

// ***************
void HybridSettings::LogScheme() const {
// ***************
  if (Compose) {
    if (RandomWalkProbGiven) {
      MACH3LOG_WARN("RandomWalkProb is set but the scheme is composition, which runs both "
                    "kernels every step and so never consults it. Ignoring it");
    }
    MACH3LOG_INFO("Hybrid sampler, composition: {} random walk sweep{} on the gradient free "
                  "block then one NUTS step on the rest, every step", Sweeps, Sweeps == 1 ? "" : "s");
    return;
  }
  MACH3LOG_INFO("Hybrid sampler, mixture: {:.0f}% random walk on the gradient free block, "
                "{:.0f}% NUTS on the rest", 100. * RandomWalkProb, 100. * (1. - RandomWalkProb));
  MACH3LOG_INFO("Each block is frozen on the steps the other one got, so expect an "
                "autocorrelation time in chain steps inflated by about {:.1f}x. Scheme "
                "'composition' avoids that", 1. / (1. - RandomWalkProb));
}
