#pragma once
#include "api/argparse.hpp"


namespace mach3 {

    class IPlugin {
        public:
            virtual ~IPlugin() = default;
            virtual int run() = 0;
            virtual MaCh3ArgumentParser* get_parser() = 0;
    };

};

// Factory function typedefs
extern "C" {
    typedef mach3::IPlugin* (*create_plugin_t)();
    typedef void (*destroy_plugin_t)(mach3::IPlugin*);
}

#define MACH3_REGISTER_PLUGIN(PluginClass)                     \
static_assert(std::is_base_of<mach3::IPlugin, PluginClass>::value, "PluginClass must derive from mach3::IPlugin");  \
extern "C" mach3::IPlugin* create_plugin() {                   \
    return new PluginClass();                                  \
}                                                              \
extern "C" void destroy_plugin(mach3::IPlugin* p) {            \
    delete p;                                                  \
}
