/// @file ConfigModule.cpp
/// @brief Reports MaCh3 installation and build configuration

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "CLI/ConfigData.hpp"
#include "CLI/Modules/ConfigModule.hpp"

namespace M3 {
    MaCh3ArgumentParser* ConfigModule::get_parser() {
        m_parser = std::make_unique<MaCh3ArgumentParser>(
            "config", "1.0", argparse::default_arguments::help);
        m_parser->add_description("Print MaCh3 installation paths, build flags, and enabled features.");
        m_parser->add_argument("--incdir").help("Print the installed header directory.").flag();
        m_parser->add_argument("--cflags").help("Print compiler flags for MaCh3 headers.").flag();
        m_parser->add_argument("--libdir").help("Print the installed library directory.").flag();
        m_parser->add_argument("--libs").help("Print MaCh3 library linker flags.").flag();
        m_parser->add_argument("--features").help("Print enabled MaCh3 features.").flag();
        m_parser->add_argument("--has-feature")
            .help("Return success if FEATURE is enabled.")
            .metavar("FEATURE");
        m_parser->add_argument("--date").help("Print the build configuration date.").flag();
        return m_parser.get();
    }

    int ConfigModule::Run() {
        std::vector<std::string> output;
        if (m_parser->get<bool>("--incdir")) output.emplace_back(M3::CLIConfig::IncludeDir);
        if (m_parser->get<bool>("--cflags")) output.emplace_back(M3::CLIConfig::Cflags);
        if (m_parser->get<bool>("--libdir")) output.emplace_back(M3::CLIConfig::LibraryDir);
        if (m_parser->get<bool>("--libs")) output.emplace_back(M3::CLIConfig::Libraries);
        if (m_parser->get<bool>("--features")) output.emplace_back(M3::CLIConfig::Features);
        if (m_parser->get<bool>("--date")) output.emplace_back(M3::CLIConfig::Date);

        const auto feature = m_parser->present<std::string>("--has-feature");
        bool feature_enabled = true;
        if (feature) {
            std::istringstream features(M3::CLIConfig::Features);
            std::string enabled_feature;
            feature_enabled = false;
            while (features >> enabled_feature) {
                if (enabled_feature == *feature) {
                    feature_enabled = true;
                    break;
                }
            }
        }

        if (output.empty() && !feature) {
            std::cout << m_parser->help().str();
            return 0;
        }

        for (std::size_t i = 0; i < output.size(); ++i) {
            if (i != 0) std::cout << ' ';
            std::cout << output[i];
        }
        if (!output.empty()) std::cout << '\n';
        return feature_enabled ? 0 : 1;
    }
}
