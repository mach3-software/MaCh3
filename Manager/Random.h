#pragma once

#include <cstdint>
#include <vector>

#include "Manager/Core.h"
#include "Manager/Monitor.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT include
#include "TRandom3.h"
_MaCh3_Safe_Include_End_ //}

namespace M3::rand
{
/// @brief Universal random number generator for MaCh3.
///
/// @note STL random-number generators were also benchmarked, but ROOT's
///       TRandom3 was faster for the tested workloads. The implementation
///       can be switched to an STL engine if needed.
///
/// @author Kamil Skwarczynski
class Random
{
 public:
  /// @brief Get the global Random instance.
  static Random& Instance()
  {
    static Random instance;
    return instance;
  }

  /// @brief Set the global seed and initialise all RNG streams.
  ///
  /// This should be called once at the beginning of the executable,
  /// before any multithreaded work is started.
  void SetSeed(const std::uint64_t seed)
  {
    fSeed = seed;
    const auto nThreads = GetNThreads();
    fEngines.clear();
    fEngines.reserve(nThreads);

    for (int thread = 0; thread < nThreads; ++thread) {
      fEngines.emplace_back(
        std::make_unique<TRandom3>(MakeThreadSeed(fSeed, thread)));
    }
  }

  /// @brief Get the global seed.
  std::uint64_t GetSeed() const {
    return fSeed;
  }

  /// @brief Get the random generator associated with the current thread.
  TRandom3* Engine() {
    return fEngines[ThreadID()].get();
  }

  /// @brief Generate a uniform random number in [min, max).
  double Uniform(const double min = 0.0,
                 const double max = 1.0) {
    return Engine()->Uniform(min, max);
  }

  /// @brief Generate a uniform random integer in [min, max].
  int UniformInt(const int min, const int max) {
    return Engine()->Integer(max - min + 1) + min;
  }

  /// @brief Generate a Gaussian random number.
  double Gaus(const double mean = 0.0,
              const double sigma = 1.0) {
    return Engine()->Gaus(mean, sigma);
  }
 private:
  /// @brief Constructor.
  Random() {
    SetSeed(0);
  }

  /// @brief Get the current OpenMP thread ID.
  static std::size_t ThreadID() {
    #ifdef MULTITHREAD
    return static_cast<std::size_t>(omp_get_thread_num());
    #else
    return 0;
    #endif
  }

  /// @brief Generate a deterministic, unique seed for each thread.
  static std::uint64_t MakeThreadSeed(const std::uint64_t seed,
                                      const std::size_t thread) {
    std::uint64_t x =
    seed + 0x9e3779b97f4a7c15ULL * (thread + 1);

    // SplitMix64 mixing.
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);

    return x;
  }

  std::uint64_t fSeed = 0;

  /// @brief One random generator per thread.
  std::vector<std::unique_ptr<TRandom3>> fEngines;
};

  /// @brief Set the global random seed.
  inline void SetSeed(const std::uint64_t seed) {
    Random::Instance().SetSeed(seed);
  }

  /// @brief Generate a uniform random number in [min, max).
  inline double Uniform(const double min = 0.0,
                        const double max = 1.0) {
    return Random::Instance().Uniform(min, max);
  }

  /// @brief Generate a uniform random integer in [min, max].
  inline int UniformInt(const int min, const int max) {
    return Random::Instance().UniformInt(min, max);
  }

  /// @brief Generate a Gaussian random number.
  inline double Gaus(const double mean = 0.0,
                     const double sigma = 1.0) {
    return Random::Instance().Gaus(mean, sigma);
  }
} // namespace M3
