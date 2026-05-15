#include "cli/MaCh3Program.hpp"
#include "cli/modules/ProcessMCMCPlugin.hpp"
#include "cli/modules/DiagMCMCPlugin.hpp"
#include "cli/modules/GetPenaltyTermPlugin.hpp"


int main(int argc, char *argv[]) {
    M3::MaCh3Program program;

    M3::ProcessMCMCPlugin proc;
    M3::DiagMCMCPlugin diag;
    M3::GetPenaltyTermPlugin penterm;

    program.add_core_plugin(proc);
    program.add_core_plugin(diag);
    program.add_core_plugin(penterm);
    program.load_dynamic_plugins();

    try {
        program.parse_args(argc, argv);
    }
    catch(...){
        return 1;
    }

    try{
        program.run();
    }
    catch(...){
        return 2;
    }

    return 0;
}
