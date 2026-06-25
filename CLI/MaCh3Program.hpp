/// @file MaCh3Program.hpp
/// @brief Core program class for the MaCh3 command-line interface

#pragma once
#include <map>
#include <string>
#include <vector>
#include <filesystem>
#include "CLI/API/plugin.hpp"
#include "CLI/DynamicPlugin.hpp"

namespace fs = std::filesystem;

namespace M3{

    /// @class MaCh3Program
    /// @brief Main program class that manages modules, plugins, and command-line parsing
    ///
    /// This class extends MaCh3ArgumentParser to provide functionality for:
    /// - Managing core modules and dynamic plugins
    /// - Installing and generating shell completions
    /// - Running selected subcommands
    class MaCh3Program: public MaCh3ArgumentParser{
        public:
            using MaCh3ArgumentParser::MaCh3ArgumentParser;

            /// @brief Destructor that unloads all dynamic plugins
            virtual ~MaCh3Program(){
                //std::cout<<"UNLOADED"<< std::endl;
                this->unload_dynamic_plugins();
            }

            /// @brief Add a core module to the program
            /// @param module Reference to the module to be added
            void add_core_module(IModule& module);

            /// @brief Load dynamic plugins from the MACH3_PLUGINS environment variable
            void load_dynamic_plugins();

            /// @brief Install shell completion scripts for bash and zsh
            const void install_completions() const;

            /// @brief Generate completion suggestions for the given prefix
            /// @param prefix The prefix string to match against available subcommands
            const void completions(const std::string& prefix) const;

            /// @brief Unload all dynamically loaded plugins
            void unload_dynamic_plugins();

            /// @brief Run the selected subcommand
            /// @return Exit code from the executed module
            int Run();

        private:
            /// @brief Detect the user's shell from the SHELL environment variable
            /// @return Shell name ("bash", "zsh", or empty string if not detected)
            const std::string detect_shell() const;

            /// @brief Write content to a file, creating parent directories if needed
            /// @param path Filesystem path to write to
            /// @param content Content to write to the file
            /// @return true if successful, false otherwise
            const bool write_file(const fs::path& path, std::string_view content) const;

            /// @brief Expand a plugin path (file or directory) into a list of .so files
            /// @param path Path to a plugin file or directory containing plugins
            /// @return Vector of plugin file paths
            std::vector<std::string> expand_plugin_path(const std::string& path) const;

        private:
            std::vector<std::string> m_subcommands;                                     ///< List of registered subcommand names
            std::map<const MaCh3ArgumentParser*, IModule*> m_module_map;                ///< Map of parsers to core modules
            std::map<const MaCh3ArgumentParser*, DynamicPlugin*> m_dynamic_plugin_map;  ///< Map of parsers to dynamic plugins

            static constexpr std::string_view BASH_COMPLETION = R"(
_mach3_complete() {
    local cur prev words cword
    _init_completion || return
    COMPREPLY=( $(mach3 --complete "$cur" "${words[@]:1:$cword}") )
}
complete -F _mach3_complete mach3
)";
            static constexpr std::string_view ZSH_COMPLETION = R"(
#compdef mach3

local -a completions
local cur="$words[-1]"
completions=("${(@f)$(mach3 --complete "$cur" "${words[@]:1}")}")
compadd -a completions
)";
    };
};
