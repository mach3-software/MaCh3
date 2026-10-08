/// @file mach3.cpp
/// @brief Main entry point for the MaCh3 command-line interface
///
/// This file contains the main function that initializes the MaCh3 program,
/// registers core modules, loads dynamic plugins, and processes command-line arguments.

#include "CLI/Modules/ProcessMCMCModule.hpp"
#include "CLI/Modules/DiagMCMCModule.hpp"
#include "CLI/Modules/GetPenaltyTermModule.hpp"
#include "CLI/MaCh3Program.hpp"

/// @brief Main entry point for the MaCh3 application
///
/// Initializes the MaCh3 program with core modules and dynamic plugins,
/// parses command-line arguments, and executes the selected subcommand.
///
/// @param argc Number of command-line arguments
/// @param argv Array of command-line argument strings
/// @return Exit code (0 on success, non-zero on error)
int main(int argc, char *argv[]) {
    M3::MaCh3Program program("mach3");

    // Hidden completion option — accepts prefix + preceding words as context
    program.add_argument("--complete")
        .hidden()
        .nargs(argparse::nargs_pattern::any);

    // Hidden installer option
    program.add_argument("--install-completions")
        .help("")
        .flag();

    M3::ProcessMCMCModule proc;
    M3::DiagMCMCModule diag;
    M3::GetPenaltyTermModule penterm;

    program.add_core_module(proc);
    program.add_core_module(diag);
    program.add_core_module(penterm);
    program.load_dynamic_plugins();

    // Pre-scan for --complete before argparse processes the full command line,
    // so incomplete commands (missing required args) still produce completions.
    for (int i = 1; i < argc; i++) {
        if (std::string_view(argv[i]) == "--complete") {
            std::string prefix = (i + 1 < argc) ? argv[i + 1] : "";
            std::vector<std::string> context;
            for (int j = i + 2; j < argc; j++) {
                context.emplace_back(argv[j]);
            }
            program.completions(prefix, context);
            return 0;
        }
    }

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {

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
        program.Run();
    }
    catch(...){
        return 3;
    }

    return 0;
}
