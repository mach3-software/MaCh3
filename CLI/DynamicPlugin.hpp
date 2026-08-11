/// @file DynamicPlugin.hpp
/// @brief Wrapper class for dynamically loaded plugins

#pragma once
#include <dlfcn.h>
#include <memory>
#include <functional>
#include "CLI/API/plugin.hpp"

namespace M3{
    /// @class DynamicPlugin
    /// @brief Manages the lifecycle of a dynamically loaded plugin
    ///
    /// This class wraps a plugin loaded from a shared library, maintaining the
    /// library handle and providing cleanup functionality using smart pointers.
    class DynamicPlugin: public IPlugin{
        public:
            /// @brief Constructor that loads a plugin from a shared library file
            /// @param sofile Path to the shared library file (.so)
            DynamicPlugin(const std::string& sofile)
                : DynamicPlugin(sofile.c_str()) {}

            /// @brief Constructor that loads a plugin from a shared library file
            /// @param sofile Path to the shared library file (.so)
            DynamicPlugin(const char* sofile)
                : m_handle(nullptr, &DynamicPlugin::dlcloser),
                  m_plugin(nullptr, nullptr){

                
                m_handle.reset(dlopen(sofile, RTLD_NOW));
                if (!m_handle) {
                    std::cerr << "dlopen failed - shared library '"<< sofile << "' not loaded: " << dlerror() << "\n";
                    throw std::runtime_error("dlopen failed");
                }

                auto create = reinterpret_cast<create_plugin_t>(dlsym(m_handle.get(), "create_plugin"));
                auto destroy = reinterpret_cast<destroy_plugin_t>(dlsym(m_handle.get(), "destroy_plugin"));

                if (!create || !destroy) {
                    std::cerr << "Invalid plugin: " << sofile << "\n";
                    m_handle.reset(nullptr);
                    throw std::runtime_error("Invalid plugin - missing create_plugin or destroy_plugin");
                }

                try{
                    m_plugin = std::unique_ptr<IPlugin, destroy_plugin_t>(create(), destroy);
                    if (!m_plugin){
                        throw std::runtime_error("null pointer");
                    }
                }
                catch (...){
                    std::cerr << "Error instantiating plugin: " << sofile << "\n";
                    m_handle.reset(nullptr);
                    throw std::runtime_error("Error instantiating plugin");
                }

                try{
                    m_parser = m_plugin->get_parser();
                    if (!m_parser){
                        throw std::runtime_error("null pointer");
                    }
                }
                catch(...){
                    std::cerr<< "Error retrieving parser" << "\n";
                    m_plugin.reset(nullptr);
                    m_handle.reset(nullptr);
                    throw std::runtime_error("Error retrieving parser");
                }
            }

            /// @brief Run the plugin's main functionality
            /// @return Exit code from the plugin
            int Run() override {
                return m_plugin->Run();
            }

            /// @brief Get the argument parser for this plugin
            /// @return Pointer to the plugin's argument parser
            ///
            /// pointer is owned by the shared library and will be deleted when the library is unloaded
            MaCh3ArgumentParser* get_parser() override {
                return m_parser;
            }

            /// @brief Destructor closes the shared library handle
            virtual ~DynamicPlugin() = default;

        private:

            static void dlcloser(void* h) {
                if (h) dlclose(h);
            }
            MaCh3ArgumentParser* m_parser = nullptr;  ///< Pointer to the plugin's argument parser, owned by the shared library
            // using dlcloser rather than built in dlclose to avoid gcc warning about casting function pointer to void pointer
            std::unique_ptr<void, decltype(&dlcloser)> m_handle;  ///< Handle to the loaded shared library
            std::unique_ptr<IPlugin, destroy_plugin_t> m_plugin;  ///< Smart pointer to plugin with custom deleter
    };
}
