#include <dlfcn.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <vector>
#include "api/plugin.hpp"
#include "Diagnostics/ProcessMCMCPlugin.hpp"
#include "Diagnostics/DiagMCMCPlugin.hpp"
#include "Diagnostics/GetPenaltyTermPlugin.hpp"

using namespace std;
namespace fs = std::filesystem;

namespace mach3{
    class DynamicPlugin: public IPlugin{
        public:
            DynamicPlugin(void* handle,
                          IPlugin* plugin,
                          destroy_plugin_t destroy_func):m_handle(handle),
                                                         m_plugin(plugin),
                                                         m_destroy_func(destroy_func){}

            int run() {
                return m_plugin->run();
            }
            MaCh3ArgumentParser* get_parser(){
                return m_plugin->get_parser();
            }

            void destroy(){
                m_destroy_func(m_plugin);
                m_plugin = 0;
                dlclose(m_handle);
                m_handle = 0;
                m_destroy_func = 0;
            }

        private:
            void* m_handle;
            IPlugin* m_plugin;
            destroy_plugin_t m_destroy_func;
    };

    class NoArgsException: public std::exception{};


    class MaCh3Program: public MaCh3ArgumentParser{
        public:
            using MaCh3ArgumentParser::MaCh3ArgumentParser;
            virtual ~MaCh3Program(){
                //std::cout<<"UNLOADED"<< std::endl;
                this->unload_dynamic_plugins();
            }

            
            void parse_args(int argc, const char *const argv[]) {
                try{
                    MaCh3ArgumentParser::parse_args(argc, argv);
                }
                catch (const std::exception& err) {
                    std::cerr << err.what() << std::endl;
                    std::cerr << this->get_subcommand_used();
                    MACH3LOG_ERROR(err.what());
                    throw;
                }

                if (!this->get_subcommand_used()){
                    std::cerr << this->get_subcommand_used();
                    throw NoArgsException();
                }
            }

            void add_core_plugin(IPlugin& plugin){
                MaCh3ArgumentParser* parser = plugin.get_parser();
                this->add_subparser(*parser);
                m_plugin_map[parser] = &plugin;
            }

            void load_dynamic_plugins(){
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

            void unload_dynamic_plugins(){
                for (auto& [_, plugin] : m_dynamic_plugin_map){
                    plugin->destroy();
                }
                m_dynamic_plugin_map.clear();
            }

            int run(){
                
                const MaCh3ArgumentParser& sub_parser = this->get_subcommand_used();
                
                if (sub_parser){
                    auto plugin_itr = m_plugin_map.find(&sub_parser);
                    if (plugin_itr != m_plugin_map.end()){
                        return plugin_itr->second->run();
                    }
                    auto dplugin_itr = m_dynamic_plugin_map.find(&sub_parser);
                    if (dplugin_itr != m_dynamic_plugin_map.end()){
                        return dplugin_itr->second->run();
                    }
                }
                return 0;        
            }

        private:
            std::vector<std::string> expand_plugin_path(const std::string& path) const {
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

        private:
            std::map<const MaCh3ArgumentParser*, IPlugin*> m_plugin_map;
            std::map<const MaCh3ArgumentParser*, DynamicPlugin*> m_dynamic_plugin_map;
    };
};


int main(int argc, char *argv[]) {
    mach3::MaCh3Program program;

    mach3::ProcessMCMCPlugin proc;
    mach3::DiagMCMCPlugin diag;
    mach3::GetPenaltyTermPlugin penterm;

    program.add_core_plugin(proc);
    program.add_core_plugin(diag);
    program.add_core_plugin(penterm);
    program.load_dynamic_plugins();

    try {
        program.parse_args(argc, argv);
    }
    catch(...){
        return 1;
    }

    try{
        program.run();
    }
    catch(...){
        return 2;
    }

    return 0;
}
