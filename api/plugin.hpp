#pragma once
#include "api/argparse.hpp"


namespace M3 {

    class IPlugin {
        public:
            virtual ~IPlugin() = default;
            virtual int run() = 0;
            virtual MaCh3ArgumentParser* get_parser() = 0;
    };

};

// Factory function typedefs
extern "C" {
    typedef M3::IPlugin* (*create_plugin_t)();
    typedef void (*destroy_plugin_t)(M3::IPlugin*);
}

#define MACH3_REGISTER_PLUGIN(PluginClass)                     \
static_assert(std::is_base_of<M3::IPlugin, PluginClass>::value, "PluginClass must derive from M3::IPlugin");  \
extern "C" M3::IPlugin* create_plugin() {                   \
    return new PluginClass();                                  \
}                                                              \
extern "C" void destroy_plugin(M3::IPlugin* p) {            \
    delete p;                                                  \
}
