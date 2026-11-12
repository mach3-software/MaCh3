/// @file mach3.cpp
/// @brief Main entry point for the MaCh3 command-line interface
///
/// This file contains the main function that initializes the MaCh3 program,
/// registers core modules, loads dynamic plugins, and processes command-line arguments.

#include "CLI/Modules/ProcessMCMCModule.hpp"
#include "CLI/Modules/DiagModule.hpp"
#include "CLI/Modules/GetPenaltyTermModule.hpp"
#include "CLI/Modules/RHatModule.hpp"
#include "CLI/Modules/SmearModule.hpp"
#include "CLI/Modules/ReweightModule.hpp"
#include "CLI/Modules/CombineModule.hpp"
#include "CLI/Modules/Plot/PlotModule.hpp"
#include "CLI/Modules/UmbrellaModule.hpp"
#include "CLI/MaCh3Program.hpp"

/// @brief Main entry point for the MaCh3 application
///
/// Initializes the MaCh3 program with core modules and dynamic plugins,
/// parses command-line arguments, and executes the selected subcommand.
///
/// @param argc Number of command-line arguments
/// @param argv Array of command-line argument strings
/// @return Exit code (0 on success, non-zero on error)
// Fix startup errors: the outer boundary also covers construction, registration and discovery.
int main(int argc, char *argv[]) try {
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
    M3::DiagModule diag;
    M3::GetPenaltyTermModule penterm;
    M3::RHatModule rhat;
    M3::SmearModule smear;
    M3::ReweightModule reweight;
    M3::CombineModule combine;
    M3::PlotModule plot;
    M3::UmbrellaModule umbrella;

    program.add_core_module(proc);
    program.add_core_module(diag);
    program.add_core_module(penterm);
    program.add_core_module(rhat);
    program.add_core_module(smear);
    program.add_core_module(reweight);
    program.add_core_module(combine);
    program.add_core_module(plot);
    program.add_core_module(umbrella);
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
        // Fix error diagnostics: exception text is data, not a logger format string.
        MACH3LOG_ERROR("{}", err.what());
        return 1;
    }
    catch(...){
        MACH3LOG_ERROR("Unknown parsing error.");
        return 1;
    }

    // Fix handling of no-argument commands: allow leaf defaults; bare groups still need a child.
    const auto& selected_parser = program.get_subcommand_used();
    if (!selected_parser && !selected_parser.subparsers().empty()){
        std::cerr << selected_parser;
        return 2;
    }


    if (program.get<bool>("--install-completions")) {
        // Fix handling of completion installation: do not turn an installer failure into process success.
        return program.install_completions();
    }

    try{
        return program.Run();
    }
    catch(const std::exception& err){
        MACH3LOG_ERROR("{}", err.what());
        std::cerr << "An exception occurred during program execution: " << err.what() << std::endl;
        return 3;
    }
    catch(...){
        MACH3LOG_ERROR("Unknown exception occurred during program execution.");
        std::cerr << "An unknown exception occurred during program execution." << std::endl;
        return 3;
    }

}
catch (const std::exception& err) {
    MACH3LOG_ERROR("CLI initialisation failed: {}", err.what());
    std::cerr << "CLI initialisation failed: " << err.what() << "\n";
    return 1;
}
catch (...) {
    MACH3LOG_ERROR("Unknown exception during CLI initialisation.");
    std::cerr << "Unknown exception during CLI initialisation.\n";
    return 1;
}
