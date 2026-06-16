#include <dlfcn.h>
#include <filesystem>
#include <iostream>
#include <sstream>
#include "Manager/MaCh3Logger.h"
#include "CLI/MaCh3Program.hpp"

using std::vector, std::string, std::map;
namespace fs = std::filesystem;

namespace M3{

    const bool MaCh3Program::write_file(const fs::path& path, std::string_view content) const{
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

    const std::string MaCh3Program::detect_shell() const{
        const char* shell = std::getenv("SHELL");
        if (!shell) return "";

        std::string s(shell);
        if (s.find("bash") != std::string::npos) return "bash";
        if (s.find("zsh")  != std::string::npos) return "zsh";

        return "";
    }

    const void MaCh3Program::install_completions() const{
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

        fs::path home_path(home);

        if (shell == "bash") {
            fs::path path = home_path / ".local/share/bash-completion/completions/mach3";
            if (write_file(path, MaCh3Program::BASH_COMPLETION))
                std::cout << "Installed bash completions → " << path << "\n";
            return;
        }

        if (shell == "zsh") {
            fs::path path = home_path / ".local/share/zsh/site-functions/_mach3";
            if (write_file(path, MaCh3Program::ZSH_COMPLETION))
                std::cout << "Installed zsh completions → " << path << "\n";
            return;
        }

    }
                
    const void MaCh3Program::completions(const std::string& prefix) const{
        for (const std::string& cmd : m_subcommands) {
            if (cmd.rfind(prefix, 0) == 0)
                std::cout << cmd << "\n";
        }
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
                void* handle = dlopen(sofile.c_str(), RTLD_NOW);
                if (!handle) {
                    std::cerr << "dlopen failed - shared library '"<< sofile << "' not loaded: " << dlerror() << "\n";
                    continue;
                }

                auto create = reinterpret_cast<create_plugin_t>(dlsym(handle, "create_plugin"));
                auto destroy = reinterpret_cast<destroy_plugin_t>(dlsym(handle, "destroy_plugin"));

                if (!create || !destroy) {
                    std::cerr << "Invalid plugin: " << sofile << "\n";
                    dlclose(handle);
                    continue;
                }

                IPlugin* plugin = nullptr;
                try{
                    plugin = create();
                    if (!plugin){
                        throw std::runtime_error("null pointer");
                    }
                }
                catch (...){
                    std::cerr << "Error instantiating plugin: " << sofile << "\n";
                    dlclose(handle);
                    continue;
                }

                MaCh3ArgumentParser* parser = nullptr;
                try{
                    parser = plugin->get_parser();
                    if (!parser){
                        throw std::runtime_error("null pointer");
                    }
                }
                catch(...){
                    std::cerr<< "Error retrieving parser" << "\n";
                    destroy(plugin);
                    dlclose(handle);
                    continue;
                }


                try{
                    this->add_subparser(*parser);
                    m_dynamic_plugin_map[parser] = new DynamicPlugin(handle, plugin, destroy);
                    m_subcommands.push_back(parser->name());
                }
                catch(...){
                    std::cerr<< "Error finalising loading of plugin." << "\n";
                    destroy(plugin);
                    dlclose(handle);
                    continue;    
                }
            }
        }
    }

    void MaCh3Program::unload_dynamic_plugins(){
        for (auto& [_, plugin] : m_dynamic_plugin_map){
            plugin->destroy();
        }
        m_dynamic_plugin_map.clear();
    }

    int MaCh3Program::run(){
        
        const MaCh3ArgumentParser& sub_parser = this->get_subcommand_used();
        
        if (sub_parser){
            auto plugin_itr = m_module_map.find(&sub_parser);
            if (plugin_itr != m_module_map.end()){
                return plugin_itr->second->run();
            }
            auto dplugin_itr = m_dynamic_plugin_map.find(&sub_parser);
            if (dplugin_itr != m_dynamic_plugin_map.end()){
                return dplugin_itr->second->run();
            }
        }
        return 0;        
    }


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