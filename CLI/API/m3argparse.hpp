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

            /// @brief Get the description of this parser/subcommand
            /// @return The description string
            const std::string& description() const{
                return this->m_description;
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
                        return static_cast<MaCh3ArgumentParser&>(subparser.get()).get_subcommand_used();
                    }
                }
                return (*this);
            }

            /// @brief Set a subcommand to use when none of the registered subcommands is given
            /// @param name Name of the default subcommand
            void set_default_subcommand(const std::string& name){
                m_default_subcommand = name;
            }

            /// @brief Insert default subcommands into args where no explicit subcommand was given
            /// @param args Full argument list
            /// @param pos Index in args of this parser's own name (0 for the program itself)
            void insert_default_subcommands(std::vector<std::string>& args, std::size_t pos = 0) const{
                const std::size_t next = pos + 1;
                if (next < args.size()) {
                    auto it = this->m_subparser_map.find(args[next]);
                    if (it != this->m_subparser_map.end()) {
                        auto& sub = static_cast<MaCh3ArgumentParser&>(it->second->get());
                        // argparse only prefixes the immediate parent's name, so rebuild the full path
                        sub.m_parser_path = this->m_parser_path + " " + sub.m_program_name;
                        sub.insert_default_subcommands(args, next);
                        return;
                    }
                }
                if (!m_default_subcommand.empty()) {
                    auto it = this->m_subparser_map.find(m_default_subcommand);
                    if (it != this->m_subparser_map.end()) {
                        // hide the default subcommand name in its usage/help line
                        static_cast<MaCh3ArgumentParser&>(it->second->get()).m_parser_path = this->m_parser_path;
                    }
                    args.insert(args.begin() + static_cast<std::ptrdiff_t>(next), m_default_subcommand);
                }
            }

        private:
            std::string m_default_subcommand;
    };
}
