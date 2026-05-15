#pragma once
#include <dlfcn.h>
#include "api/plugin.hpp"

namespace M3{
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
};