/**
 * @file DynamicPlugin.hpp
 * @brief Wrapper class for dynamically loaded plugins
 */

#pragma once
#include <dlfcn.h>
#include "CLI/API/plugin.hpp"

namespace M3{
    /**
     * @class DynamicPlugin
     * @brief Manages the lifecycle of a dynamically loaded plugin
     *
     * This class wraps a plugin loaded from a shared library, maintaining the
     * library handle and providing cleanup functionality.
     */
    class DynamicPlugin: public IPlugin{
        public:
            /**
             * @brief Constructor for a dynamic plugin
             * @param handle Handle to the loaded shared library
             * @param plugin Pointer to the instantiated plugin
             * @param destroy_func Function pointer to destroy the plugin
             */
            DynamicPlugin(void* handle,
                          IPlugin* plugin,
                          destroy_plugin_t destroy_func):m_handle(handle),
                                                         m_plugin(plugin),
                                                         m_destroy_func(destroy_func){}

            /**
             * @brief Run the plugin's main functionality
             * @return Exit code from the plugin
             */
            int run() {
                return m_plugin->run();
            }

            /**
             * @brief Get the argument parser for this plugin
             * @return Pointer to the plugin's argument parser
             */
            MaCh3ArgumentParser* get_parser(){
                return m_plugin->get_parser();
            }

            /**
             * @brief Destroy the plugin and close the shared library
             *
             * Calls the plugin's destroy function and closes the library handle.
             */
            void destroy(){
                m_destroy_func(m_plugin);
                m_plugin = 0;
                dlclose(m_handle);
                m_handle = 0;
                m_destroy_func = 0;
            }

        private:
            void* m_handle;                      ///< Handle to the loaded shared library
            IPlugin* m_plugin;                   ///< Pointer to the plugin instance
            destroy_plugin_t m_destroy_func;     ///< Function pointer to destroy the plugin
    };
};