#pragma once
#include <map>
#include <string>
#include <vector>
#include <filesystem>
#include "cli/api/plugin.hpp"
#include "cli/DynamicPlugin.hpp"

namespace fs = std::filesystem;

namespace M3{

    // class NoArgsException: public std::exception{};

    class MaCh3Program: public MaCh3ArgumentParser{
        public:
            using MaCh3ArgumentParser::MaCh3ArgumentParser;
            virtual ~MaCh3Program(){
                //std::cout<<"UNLOADED"<< std::endl;
                this->unload_dynamic_plugins();
            }          
            // void parse_args(int argc, const char *const argv[]);
            void add_core_module(IModule& module);
            void load_dynamic_plugins();
            const void install_completions() const;
            const void completions(const std::string& prefix) const;
            void unload_dynamic_plugins();
            int run();

        private:
            const std::string detect_shell() const;
            const bool write_file(const fs::path& path, std::string_view content) const;
            std::vector<std::string> expand_plugin_path(const std::string& path) const;

        private:
            std::vector<std::string> m_subcommands;
            std::map<const MaCh3ArgumentParser*, IModule*> m_module_map;
            std::map<const MaCh3ArgumentParser*, DynamicPlugin*> m_dynamic_plugin_map;

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