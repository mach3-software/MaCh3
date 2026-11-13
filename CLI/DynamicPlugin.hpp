/// @file DynamicPlugin.hpp
/// @brief Wrapper class for dynamically loaded plugins

#pragma once
#include <dlfcn.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
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

                try{
                    // Fix API compatibility: reject legacy/mismatched plugins before calling their factory.
                    auto version = reinterpret_cast<plugin_api_version_t>(dlsym(m_handle.get(), "plugin_api_version"));
                    if (!version) {
                        throw std::runtime_error("Missing plugin_api_version function; rebuild the plugin with the current MaCh3 SDK");
                    }
                    const auto plugin_version = version();
                    if (plugin_version != PLUGIN_API_VERSION) {
                        throw std::runtime_error("Plugin API version " + std::to_string(plugin_version) +
                                                 " does not match host version " + std::to_string(PLUGIN_API_VERSION) +
                                                 "; rebuild the plugin with the current MaCh3 SDK");
                    }
                }
                catch (const std::exception& err){
                    // Fix plugin lifetime: translate exceptions before unloading plugin-defined exception code.
                    throw std::runtime_error("Cannot initialise plugin '" + std::string(sofile) + "': " + err.what());
                }
                catch (...){
                    throw std::runtime_error("Cannot initialise plugin '" + std::string(sofile) + "': unknown exception");
                }

                auto create = reinterpret_cast<create_plugin_t>(dlsym(m_handle.get(), "create_plugin"));
                auto destroy = reinterpret_cast<destroy_plugin_t>(dlsym(m_handle.get(), "destroy_plugin"));

                if (!create || !destroy) {
                    std::cerr << "Invalid plugin: " << sofile << "\n";
                    throw std::runtime_error("Invalid plugin - missing create_plugin or destroy_plugin");
                }

                try{
                    m_plugin = std::unique_ptr<IPlugin, destroy_plugin_t>(create(), destroy);
                    if (!m_plugin){
                        throw std::runtime_error("Plugin factory returned a null pointer");
                    }
                }
                catch (const std::exception& err){
                    throw std::runtime_error("Error instantiating plugin '" + std::string(sofile) + "': " + err.what());
                }
                catch (...){
                    std::cerr << "Error instantiating plugin: " << sofile << "\n";
                    throw std::runtime_error("Error instantiating plugin '" + std::string(sofile) + "': unknown exception");
                }

                try{
                    m_parser = m_plugin->get_parser();
                    if (!m_parser){
                        throw std::runtime_error("Plugin returned a null parser");
                    }
                }
                catch (const std::exception& err){
                    // Fix plugin lifetime: translate exceptions before unloading plugin-defined exception code.
                    throw std::runtime_error("Error retrieving parser '" + std::string(sofile) + "': " + err.what());
                }
                catch (...){
                    std::cerr<< "Error retrieving parser" << "\n";
                    throw std::runtime_error("Error retrieving parser '" + std::string(sofile) + "': unknown exception");
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
