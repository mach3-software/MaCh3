/// @file plugin.hpp
/// @brief Plugin interface and registration macros for MaCh3

#pragma once
#include <memory>
#include "CLI/API/argparse.hpp"


namespace M3 {

    /// @class IPlugin
    /// @brief Interface for MaCh3 plugins and modules
    ///
    /// All plugins and modules must implement this interface to provide
    /// argument parsing and execution functionality.
    class IPlugin {
        public:
            virtual ~IPlugin() = default;

            /// @brief Execute the plugin's main functionality
            /// @return Exit code (0 on success, non-zero on error)
            virtual int Run() = 0;

            /// @brief Get the argument parser for this plugin
            /// @return Pointer to the plugin's MaCh3ArgumentParser
            virtual MaCh3ArgumentParser* get_parser() = 0;
    };

    /// @typedef IModule
    /// @brief Alias for IPlugin, used for core modules
    typedef IPlugin IModule;

    /// @class IModuleBase
    /// @brief Base class for core modules in MaCh3
    ///
    /// This class provides a common interface for core modules, allowing them
    /// to be managed and executed by the MaCh3 program.
    class PluginBase : public IPlugin {
        public:
            virtual ~PluginBase() = default;

        protected:
            /// @brief Argument parser for this plugin
            std::unique_ptr<MaCh3ArgumentParser> m_parser;  ///< Argument parser for this plugin
    };



    /// @typedef IModuleBase
    /// @brief Alias for IPluginBase, used for core modules
    typedef PluginBase IModuleBase;

};

/// @typedef create_plugin_t
/// @brief Function pointer type for plugin factory function
extern "C" {
    typedef M3::IPlugin* (*create_plugin_t)();

    /// @typedef destroy_plugin_t
    /// @brief Function pointer type for plugin destructor function
    typedef void (*destroy_plugin_t)(M3::IPlugin*);
}

/// @def MACH3_REGISTER_PLUGIN
/// @brief Macro to register a plugin class with MaCh3
///
/// This macro generates the required factory and destructor functions
/// for a plugin class. Use it in your plugin implementation file.
///
/// @param PluginClass The plugin class that implements M3::IPlugin
///
/// Example usage:
/// @code
/// class MyPlugin : public M3::IPlugin {
///     // ... implementation
/// };
/// MACH3_REGISTER_PLUGIN(MyPlugin)
/// @endcode
#define MACH3_REGISTER_PLUGIN(PluginClass)                     \
static_assert(std::is_base_of<M3::IPlugin, PluginClass>::value, "PluginClass must derive from M3::IPlugin");  \
extern "C" M3::IPlugin* create_plugin() {                   \
    return new PluginClass();                                  \
}                                                              \
extern "C" void destroy_plugin(M3::IPlugin* p) {            \
    delete p;                                                  \
}
