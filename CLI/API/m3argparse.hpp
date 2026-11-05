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

            /// @brief Complete subcommands or long options at any depth without parsing incomplete arguments
            /// @return Matching names, or an empty list to allow shell file completion
            std::vector<std::string> get_completions(const std::string& prefix,
                                                    const std::vector<std::string>& context) const {
                const MaCh3ArgumentParser* parser = this;
                bool positional = false;
                for (std::size_t i = 0; i < context.size(); ++i) {
                    const auto& word = context[i];
                    if (word == "--") return {};

                    if (!positional) {
                        auto sub = parser->m_subparser_map.find(word);
                        if (sub != parser->m_subparser_map.end()) {
                            parser = &static_cast<const MaCh3ArgumentParser&>(sub->second->get());
                            continue;
                        }
                        auto fallback = parser->m_subparser_map.find(parser->m_default_subcommand);
                        if (fallback != parser->m_subparser_map.end())
                            parser = &static_cast<const MaCh3ArgumentParser&>(fallback->second->get());
                    }

                    if (!word.empty() && word[0] == '-') {
                        const auto equals = word.find('=');
                        if (equals == std::string::npos && parser->option_takes_value(word)) {
                            if (++i == context.size()) return {};
                        }
                    } else {
                        // Once positional arguments start, later filenames are not subcommands.
                        positional = true;
                    }
                }

                std::vector<std::string> names;
                if (!prefix.empty() && prefix[0] == '-') {
                    auto fallback = parser->m_subparser_map.find(parser->m_default_subcommand);
                    if (fallback != parser->m_subparser_map.end())
                        parser = &static_cast<const MaCh3ArgumentParser&>(fallback->second->get());
                    for (const auto& name : parser->get_long_option_names()) {
                        if (name.rfind(prefix, 0) == 0) names.push_back(name);
                    }
                } else if (!positional) {
                    for (const auto& [name, sub] : parser->m_subparser_map) {
                        const auto& child = static_cast<const MaCh3ArgumentParser&>(sub->get());
                        if (!child.m_suppress && name.rfind(prefix, 0) == 0)
                            names.push_back(name);
                    }
                }
                return names;
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
