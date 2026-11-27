/// @file ConfigModule.hpp
/// @brief MaCh3 installation and build configuration subcommand

#pragma once
#include "CLI/API/plugin.hpp"

namespace M3 {
    class ConfigModule final : public PluginBase {
        public:
            MaCh3ArgumentParser* get_parser() override;
            int Run() override;
    };
}
