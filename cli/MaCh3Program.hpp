#pragma once
#include <map>
#include <string>
#include <vector>
#include "api/plugin.hpp"
#include "cli/DynamicPlugin.hpp"


namespace M3{

    class NoArgsException: public std::exception{};

    class MaCh3Program: public MaCh3ArgumentParser{
        public:
            using MaCh3ArgumentParser::MaCh3ArgumentParser;
            virtual ~MaCh3Program(){
                //std::cout<<"UNLOADED"<< std::endl;
                this->unload_dynamic_plugins();
            }          
            void parse_args(int argc, const char *const argv[]);
            void add_core_plugin(IPlugin& plugin);
            void load_dynamic_plugins();
            void unload_dynamic_plugins();
            int run();

        private:
            std::vector<std::string> expand_plugin_path(const std::string& path) const;

        private:
            std::map<const MaCh3ArgumentParser*, IPlugin*> m_plugin_map;
            std::map<const MaCh3ArgumentParser*, DynamicPlugin*> m_dynamic_plugin_map;
    };
};