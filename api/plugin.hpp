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

// Factory function signature
extern "C" {
    typedef mach3::IPlugin* (*create_plugin_t)();
    typedef void (*destroy_plugin_t)(mach3::IPlugin*);
}
