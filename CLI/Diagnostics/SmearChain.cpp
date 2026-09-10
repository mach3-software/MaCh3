#include "Manager/MaCh3Logger.h"
#include "CLI/Modules/SmearModule.hpp"

int main(int argc, char const* argv[]){
    MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
    MACH3LOG_WARN("                   : you can use 'mach3 smear' as a direct replacement instead.");
    argv[0] = "smear";
    M3::SmearModule proc;
    ArgumentParser* parser = proc.get_parser();
    parser->parse_args(argc, argv);
    proc.Run();
    MACH3LOG_WARN("Deprecation Warning: Use of the standalone executable will be deprecated in future releases.");
    MACH3LOG_WARN("                   : you can use 'mach3 smear' as a direct replacement instead.");
    return 0;
}
