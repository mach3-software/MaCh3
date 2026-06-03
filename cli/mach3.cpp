#include "cli/MaCh3Program.hpp"
#include "cli/modules/ProcessMCMCPlugin.hpp"
#include "cli/modules/DiagMCMCPlugin.hpp"
#include "cli/modules/GetPenaltyTermPlugin.hpp"


int main(int argc, char *argv[]) {
    M3::MaCh3Program program("mach3");

    // Hidden completion option
    program.add_argument("--complete")
        .hidden()
        .nargs(1);

    // Hidden installer option
    program.add_argument("--install-completions")
        .help("")
        .flag();

    M3::ProcessMCMCPlugin proc;
    M3::DiagMCMCPlugin diag;
    M3::GetPenaltyTermPlugin penterm;

    program.add_core_module(proc);
    program.add_core_module(diag);
    program.add_core_module(penterm);
    program.load_dynamic_plugins();

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        // completions count as parsing error
        if (auto prefix = program.present<std::string>("--complete")) {
            program.completions(*prefix);
            return 0;
        }

        std::cerr << err.what() << std::endl;
        std::cerr << program.get_subcommand_used();
        MACH3LOG_ERROR(err.what());
        return 1;
    }
    catch(...){
        MACH3LOG_ERROR("Unknown parsing error.");
        return 1;
    }

    // no subcommand args
    if (!program.get_subcommand_used()){
        std::cerr << program.get_subcommand_used();
        return 2;
    }


    if (program.get<bool>("--install-completions")) {
        program.install_completions();
        return 0;
    }

    try{
        program.run();
    }
    catch(...){
        return 3;
    }

    return 0;
}
