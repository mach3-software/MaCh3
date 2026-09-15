/// @file SmearModule.hpp
/// @brief Allows you to smear contour. For example after performing sets of study one finds out that used sets of uncertainty doesn't fully cover analysis need. Then one can smear additionally contour.

#pragma once
#include "CLI/API/plugin.hpp"

namespace M3{

  /// @class SmearModule
  /// @brief Allows you to smear contour. For example after performing sets of study one finds out that used sets of uncertainty doesn't fully cover analysis need. Then one can smear additionally contour.
  class SmearModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~SmearModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;


      void SmearChain(const std::string& inputFile, const std::string& config);
      /// @brief Execute the MCMC diagnostics
      /// @return Exit code (0 on success)
      int Run() override;
  };
}
