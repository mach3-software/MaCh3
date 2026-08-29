/// @file m3argparse.hpp
/// @brief MaCh3 wrapper for argparse library with additional functionality

#pragma once
#include <argparse/argparse.hpp>

using argparse::ArgumentParser;

namespace M3 {

    /// @class MaCh3ArgumentParser
    /// @brief Extended ArgumentParser with MaCh3-specific functionality
    ///
    /// This class extends the standard ArgumentParser to provide additional
    /// methods for accessing parser name, subparsers, and tracking which
    /// subcommand was used.
    class MaCh3ArgumentParser: public ArgumentParser{
        public:
            using ArgumentParser::ArgumentParser;
            virtual ~MaCh3ArgumentParser() = default;

            /// @brief Get the name of this parser/subcommand
            /// @return The program or subcommand name
            const std::string name() const{
                return this->m_program_name;
            }

            /// @brief Get the list of registered subparsers
            /// @return Reference to the list of subparsers
            const std::list<std::reference_wrapper<ArgumentParser>>& subparsers() const{
                return this->m_subparsers;
            }

            /// @brief Get the subcommand that was used in parsing
            ///
            /// Recursively traverses the subparser hierarchy to find the
            /// deepest subcommand that was actually invoked.
            ///
            /// @return Reference to the MaCh3ArgumentParser of the used subcommand
            const MaCh3ArgumentParser& get_subcommand_used() const{
                for (const std::reference_wrapper<ArgumentParser>& subparser : this->m_subparsers) {
                    if (this->is_subcommand_used(subparser.get())) {
                        return static_cast<MaCh3ArgumentParser*>(&(subparser.get()))->get_subcommand_used();
                    }
                }
                return (*this);
            }

            /// @brief Return all long option names (--xxx) registered with this parser
            std::vector<std::string> get_long_option_names() const {
                std::vector<std::string> names;
                for (const auto& [name, _] : this->m_argument_map) {
                    if (name.size() >= 2 && name[0] == '-' && name[1] == '-')
                        names.push_back(name);
                }
                return names;
            }

            /// @brief Return true if the named option expects a value (not a flag)
            bool option_takes_value(const std::string& name) const {
                auto it = this->m_argument_map.find(name);
                if (it == this->m_argument_map.end()) return false;
                // get_usage_full() appends a metavar when max nargs > 0
                return it->second->get_usage_full().find(' ') != std::string::npos;
            }
    };
}
