#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/RHatModule.hpp"

int main(int argc, char const* argv[]){
    MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
    MACH3LOG_WARN("                   : you can use 'mach3 rhat' as a direct replacement instead.");
    // Create a new argv array with an additional element for "--high-mem"
    std::vector<const char*> new_argv;
    new_argv.push_back("rhat");
    new_argv.push_back("--high-mem");
    for (int i = 1; i < argc; ++i) {
        new_argv.push_back(argv[i]);
    }
    argc = static_cast<int>(new_argv.size());
    argv = new_argv.data();
    M3::RHatModule proc;
    ArgumentParser* parser = proc.get_parser();
    parser->parse_args(argc, argv);
    proc.Run();
    MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
    MACH3LOG_WARN("                   : you can use 'mach3 rhat' as a direct replacement instead.");
    return 0;
}
