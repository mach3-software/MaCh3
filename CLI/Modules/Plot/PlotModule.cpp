#include "CLI/Modules/Plot/PlotModule.hpp"

namespace M3{

    PlotModule::~PlotModule() = default;

    MaCh3ArgumentParser* PlotModule::get_parser(){
        m_parser = std::make_unique<MaCh3ArgumentParser>("plot", "1.0", argparse::default_arguments::help);
        m_parser->add_description("Plotting subcommand.");
        m_parser->add_subparser(*m_llhMapModule.get_parser());
        m_parser->add_subparser(*m_llhModule.get_parser());
        m_parser->add_subparser(*m_postFitParamModule.get_parser());
        m_parser->add_subparser(*m_matrixModule.get_parser());
        m_parser->add_subparser(*m_sigmaVariation.get_parser());
        m_parser->add_subparser(*m_predictiveModule.get_parser());
        return m_parser.get();
    }

    int PlotModule::Run() {
        // if (m_parser->is_subcommand_used("llhmap")) return m_llhMapModule.Run();
        // if (m_parser->is_subcommand_used("llh")) return m_llhModule.Run();
        // if (m_parser->is_subcommand_used("postfit")) return m_postFitParamModule.Run();
        // if (m_parser->is_subcommand_used("matrix")) return m_matrixModule.Run();
        // if (m_parser->is_subcommand_used("sigmavar")) return m_sigmaVariation.Run();
        // if (m_parser->is_subcommand_used("predictive")) return m_predictiveModule.Run();
        return 0;
    }
} // namespace M3