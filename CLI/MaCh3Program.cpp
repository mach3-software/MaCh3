/// @file MaCh3Program.cpp
/// @brief Implementation of the MaCh3Program class

#include <dlfcn.h>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include "Manager/MaCh3Logger.h"
#include "CLI/MaCh3Program.hpp"

using std::vector, std::string, std::map;
namespace fs = std::filesystem;

namespace M3{
    bool MaCh3Program::write_file(const fs::path& path, std::string_view content) const{
        try {
            fs::create_directories(path.parent_path());
            std::ofstream out(path);
            if (!out) return false;
            out << content;
            return true;
        } catch (...) {
            return false;
        }
    }

    std::string MaCh3Program::detect_shell() const{
        const char* shell = std::getenv("SHELL");
        if (!shell) return "";

        std::string s(shell);
        if (s.find("bash") != std::string::npos) return "bash";
        if (s.find("zsh")  != std::string::npos) return "zsh";

        return "";
    }

    void MaCh3Program::install_completions() const{
        const char* home = std::getenv("HOME");
        if (!home) {
            std::cerr << "Cannot determine HOME directory\n";
            return;
        }

        std::string shell = detect_shell();
        if (shell.empty()) {
            std::cerr << "Could not detect your shell. Supported: bash, zsh\n";
            return;
        }

        // Resolve the absolute path of this binary so the completion script
        // works even when mach3 is not on PATH.
        std::string binary_path = "mach3";
        try {
            binary_path = fs::read_symlink("/proc/self/exe").string();
        } catch (...) {}

        // Replace the generic "mach3" invocation in the template with the
        // absolute binary path so the script is self-contained.
        auto embed_path = [&](std::string_view tmpl) -> std::string {
            std::string s(tmpl);
            const std::string placeholder = "mach3 --complete";
            const std::string replacement = binary_path + " --complete";
            for (std::size_t pos = s.find(placeholder); pos != std::string::npos;
                 pos = s.find(placeholder, pos + replacement.size()))
                s.replace(pos, placeholder.size(), replacement);
            return s;
        };

        fs::path home_path(home);

        if (shell == "bash") {
            fs::path path = home_path / ".local/share/bash-completion/completions/mach3";
            if (write_file(path, embed_path(MaCh3Program::BASH_COMPLETION)))
                std::cout << "Installed bash completions → " << path << "\n";
            return;
        }

        if (shell == "zsh") {
            fs::path path = home_path / ".local/share/zsh/site-functions/_mach3";
            if (write_file(path, embed_path(MaCh3Program::ZSH_COMPLETION)))
                std::cout << "Installed zsh completions → " << path << "\n";
            return;
        }
    }

    void MaCh3Program::completions(const std::string& prefix, const std::vector<std::string>& context) const {
        // Find the active subcommand from previously typed context words
        std::string active_subcommand;
        for (const auto& word : context) {
            for (const auto& sub : m_subcommands) {
                if (word == sub) {
                    active_subcommand = sub;
                    break;
                }
            }
            if (!active_subcommand.empty()) break;
        }

        if (active_subcommand.empty()) {
            // No subcommand typed yet — complete top-level subcommand names
            for (const auto& cmd : m_subcommands) {
                if (cmd.rfind(prefix, 0) == 0)
                    std::cout << cmd << "\n";
            }
            return;
        }

        // Find the parser for the active subcommand
        const MaCh3ArgumentParser* sub_parser = nullptr;
        for (const auto& [parser, mod] : m_module_map) {
            if (parser->name() == active_subcommand) { sub_parser = parser; break; }
        }
        if (!sub_parser) {
            for (const auto& [parser, plug] : m_dynamic_plugin_map) {
                if (parser->name() == active_subcommand) { sub_parser = parser; break; }
            }
        }
        if (!sub_parser) return; // unknown subcommand — fall back to files

        // If the last context word is an option that expects a value argument,
        // output nothing so the shell falls back to file completion.
        if (!context.empty()) {
            const std::string& prev = context.back();
            if (prev.size() >= 2 && prev[0] == '-') {
                if (sub_parser->option_takes_value(prev)) {
                    return; // option takes a value — let the shell complete files
                }
            }
        }

        // Complete long option names when the prefix starts with '-'
        if (!prefix.empty() && prefix[0] == '-') {
            for (const auto& name : sub_parser->get_long_option_names()) {
                if (name.rfind(prefix, 0) == 0)
                    std::cout << name << "\n";
            }
            return;
        }
        // Non-option prefix (or empty) — output nothing, shell falls back to files
    }

    // void MaCh3Program::parse_args(int argc, const char *const argv[]) {
    //     try{
    //         MaCh3ArgumentParser::parse_args(argc, argv);

    //     }
    //     catch (const std::exception& err) {
    //         if (auto prefix = this->present<std::string>("--complete")) {
    //             this->completions(*prefix);
    //             exit(0);
    //         }

    //         if (this->get<bool>("--install-completions")) {
    //             this->install_completions();
    //             exit(0);
    //         }
    //         std::cerr << err.what() << std::endl;
    //         std::cerr << this->get_subcommand_used();
    //         MACH3LOG_ERROR(err.what());
    //         throw;
    //     }

    //     if (!this->get_subcommand_used()){
    //         std::cerr << this->get_subcommand_used();
    //         throw NoArgsException();
    //     }
    // }

    void MaCh3Program::add_core_module(IModule& module){
        MaCh3ArgumentParser* parser = module.get_parser();
        this->add_subparser(*parser);
        m_module_map[parser] = &module;
        m_subcommands.push_back(parser->name());
    }

    /// @brief Load dynamic plugins from paths specified in MACH3_PLUGINS environment variable
    ///
    /// Reads the MACH3_PLUGINS environment variable (colon-separated paths), loads each .so file,
    /// and registers the plugins with the program. Each plugin must export create_plugin and
    /// destroy_plugin functions.
    void MaCh3Program::load_dynamic_plugins(){
        const char* env = std::getenv("MACH3_PLUGINS");
        if (!env) {
            //std::cerr << "No plugins specified\n";
            return;
        }

        std::stringstream ss(env);
        std::string path;

        while (std::getline(ss, path, ':')) {
            for (const auto& sofile : this->expand_plugin_path(path)) {

                try{
                    std::unique_ptr<DynamicPlugin> dlplugin = std::make_unique<DynamicPlugin>(sofile);
                    MaCh3ArgumentParser* parser = dlplugin->get_parser();
                    if (!parser) {
                        std::cerr << "Plugin " << sofile << " did not provide a valid parser.\n";
                        continue;
                    }

                    this->add_subparser(*parser);
                    m_dynamic_plugin_map[parser] = std::move(dlplugin);
                    m_subcommands.push_back(parser->name());
                }
                catch (const std::exception& e) {
                    std::cerr << "Failed to load plugin from " << sofile << ": " << e.what() << "\n";
                    continue;
                }
                catch (...) {
                    std::cerr << "Failed to load plugin from " << sofile << ": unknown error\n";
                    continue;
                }
            }
        }
    }

    /// @brief Unload all dynamic plugins and clean up resources
    ///
    /// Clears the plugin map. Smart pointers automatically handle cleanup.
    void MaCh3Program::unload_dynamic_plugins(){
        m_dynamic_plugin_map.clear();
    }

    /// @brief Execute the selected subcommand
    ///
    /// Determines which subcommand was invoked and runs the corresponding module or plugin.
    ///
    /// @return Exit code from the executed module (0 on success)
    int MaCh3Program::Run(){
        
        const MaCh3ArgumentParser& sub_parser = this->get_subcommand_used();
        
        if (sub_parser){
            auto plugin_itr = m_module_map.find(&sub_parser);
            if (plugin_itr != m_module_map.end()){
                return plugin_itr->second->Run();
            }
            auto dplugin_itr = m_dynamic_plugin_map.find(&sub_parser);
            if (dplugin_itr != m_dynamic_plugin_map.end()){
                return dplugin_itr->second->Run();
            }
        }
        return 0;        
    }


    /// @brief Expand a plugin path to a list of shared library files
    ///
    /// If the path is a file, returns it directly. If it's a directory,
    /// returns all .so files found in that directory.
    ///
    /// @param path Path to a plugin file or directory
    /// @return Vector of absolute paths to .so files
    std::vector<std::string> MaCh3Program::expand_plugin_path(const std::string& path) const {
        std::vector<std::string> result;

        fs::path p(path);

        if (!fs::exists(p)) {
            std::cerr << "Plugin path does not exist: " << path << "\n";
            return result;
        }

        if (fs::is_regular_file(p)) {
            // Single .so file
            result.push_back(path);
        }
        else if (fs::is_directory(p)) {
            // Load all .so files in directory
            for (auto& entry : fs::directory_iterator(p)) {
                if (entry.is_regular_file() && entry.path().extension() == ".so") {
                    result.push_back(entry.path().string());
                }
            }
        }
        else {
            std::cerr << "Invalid plugin path (not file or directory): " << path << "\n";
        }

        return result;
    }

};
