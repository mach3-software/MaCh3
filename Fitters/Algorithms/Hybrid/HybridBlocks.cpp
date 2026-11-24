#include "Algorithms/Hybrid/HybridBlocks.h"
#include "Parameters/ParameterHandlerGeneric.h"
#include "Samples/SampleHandlerInterface.h"

#include <cctype>
#include <regex>

namespace {
  // ***************
  /// @brief Shell style wildcard match, so a config can ask for "TotalEScale*"
  bool NameMatches(const std::string& Pattern, const std::string& Name) {
  // ***************
    if (Pattern == Name) return true;
    if (Pattern.find('*') == std::string::npos) return false;

    std::string Expr;
    Expr.reserve(Pattern.size() * 2);
    for (const char C : Pattern) {
      if (C == '*') {
        Expr += ".*";
      } else if (std::isalnum(static_cast<unsigned char>(C)) || C == '_') {
        Expr += C;
      } else {
        // Anything else could be a regex metacharacter, so escape it
        Expr += '\\';
        Expr += C;
      }
    }
    return std::regex_match(Name, std::regex(Expr));
  }
}

// *******************
void HybridBlocks::Build(const std::vector<ParameterHandlerBase*>& Systematics,
                         const std::vector<SampleHandlerInterface*>& Samples, const HybridSettings& Settings) {
// *******************
  //HW Automatic blocks: a parameter any sample can't differentiate with respect to has no force
  //to integrate, so the random walk takes it. The samples' gradients are on by now, so they know
  //their coverage
  if (Settings.AutoTypes) {
    FindParametersWithoutGradient(Systematics, Samples);
    MACH3LOG_INFO("Hybrid blocks chosen automatically: every parameter some sample has no gradient "
                  "for goes to the random walk. Set RandomWalkTypes to choose by type instead");
  }

  // One mask per handler for each kernel. Every parameter is set in exactly one of the two masks,
  // fixed ones included; that is masking, not ownership, and the counts below treat a fixed
  // parameter as belonging to neither block
  RandomWalkMask.assign(Systematics.size(), {});
  NutsMask.assign(Systematics.size(), {});
  NRandomWalk_ = NNuts_ = NFixed = 0;
  for (size_t s = 0; s < Systematics.size(); ++s) {
    const int NPars = Systematics[s]->GetNumParams();
    RandomWalkMask[s].assign(static_cast<size_t>(NPars), 0);
    NutsMask[s].assign(static_cast<size_t>(NPars), 0);
    for (int i = 0; i < NPars; ++i) {
      const bool RandomWalk = BelongsToRandomWalk(*Systematics[s], s, i, Settings);
      RandomWalkMask[s][static_cast<size_t>(i)] = RandomWalk ? 1 : 0;
      NutsMask[s][static_cast<size_t>(i)] = RandomWalk ? 0 : 1;

      // A fixed parameter is in neither block's working set. Counting it would overstate what each
      // kernel is actually moving, which matters when reading the rates later
      if (Systematics[s]->IsParameterFixed(i)) ++NFixed;
      else if (RandomWalk) ++NRandomWalk_;
      else ++NNuts_;
    }
  }

  if (NFixed > 0) {
    MACH3LOG_INFO("Hybrid split: {} parameters on the random walk, {} on NUTS, {} fixed and in "
                  "neither", NRandomWalk_, NNuts_, NFixed);
  } else {
    MACH3LOG_INFO("Hybrid split: {} parameters on the random walk, {} on NUTS", NRandomWalk_, NNuts_);
  }
}

// *******************
void HybridBlocks::FindParametersWithoutGradient(const std::vector<ParameterHandlerBase*>& Systematics,
                                                 const std::vector<SampleHandlerInterface*>& Samples) {
// *******************
  NoGradient.assign(Systematics.size(), {});
  for (size_t s = 0; s < Systematics.size(); ++s) {
    NoGradient[s].assign(static_cast<size_t>(Systematics[s]->GetNumParams()), 0);
    // A parameter lacks a gradient if any sample driven by this handler can't provide one
    for (const auto* Sample : Samples) {
      if (Sample->GetGradParHandler() != Systematics[s]) continue;
      const std::vector<bool> Covered = Sample->GetGradientCoverage();
      for (size_t i = 0; i < Covered.size() && i < NoGradient[s].size(); ++i) {
        if (!Covered[i]) NoGradient[s][i] = 1;
      }
    }
  }
}

// *******************
bool HybridBlocks::BelongsToRandomWalk(const ParameterHandlerBase& Handler, const size_t s, const int i,
                                       const HybridSettings& Settings) const {
// *******************
  const std::string Name = Handler.GetParFancyName(i);
  for (const auto& Pattern : Settings.Names) {
    if (NameMatches(Pattern, Name)) return true;
  }

  if (Settings.AutoTypes) return NoGradient[s][static_cast<size_t>(i)] != 0;
  if (Settings.Types.empty()) return false;

  // Only a generic handler knows what type each of its parameters is. Anything else can still be
  // selected by name, which is why that comes first
  const auto* Generic = dynamic_cast<const ParameterHandlerGeneric*>(&Handler);
  if (Generic == nullptr) return false;

  const SystType Type = Generic->GetParamType(i);
  for (const SystType Wanted : Settings.Types) {
    if (Type == Wanted) return true;
  }
  return false;
}

// *******************
bool HybridBlocks::HandlerHasRandomWalk(const size_t s) const {
// *******************
  for (const char Flag : RandomWalkMask[s]) {
    if (Flag != 0) return true;
  }
  return false;
}

// *******************
void HybridBlocks::RequireBothBlocks(const double RandomWalkProb) const {
// *******************
  if (NRandomWalk_ == 0) {
    MACH3LOG_ERROR("No parameter matched RandomWalkTypes or RandomWalkNames, so the random walk "
                   "kernel has nothing to sample and {:.0f}% of steps would do nothing. Run "
                   "plain NUTS instead, or widen the selection", 100. * RandomWalkProb);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (NNuts_ == 0) {
    MACH3LOG_ERROR("Every parameter was put on the random walk, so NUTS has nothing to sample. "
                   "Run plain MR2T2 instead, which is the same sampler without the overhead");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// *******************
void HybridBlocks::LogRandomWalkMembers(const std::vector<ParameterHandlerBase*>& Systematics) const {
// *******************
  for (size_t s = 0; s < Systematics.size(); ++s) {
    int NThisHandler = 0;
    for (const char Flag : RandomWalkMask[s]) NThisHandler += (Flag != 0);
    if (NThisHandler == 0) continue;
    MACH3LOG_INFO("  {} contributes {} parameters to the random walk", Systematics[s]->GetName(), NThisHandler);
    for (int i = 0; i < Systematics[s]->GetNumParams(); ++i) {
      if (!InRandomWalk(s, i)) continue;
      MACH3LOG_DEBUG("    {}", Systematics[s]->GetParFancyName(i));
    }
  }
}
