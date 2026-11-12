#include "CLI/Modules/DiagModule.hpp"

namespace M3{

DiagModule::~DiagModule() = default;

MaCh3ArgumentParser* DiagModule::get_parser(){
  m_parser = std::make_unique<MaCh3ArgumentParser>("diagnostics", "1.0", argparse::default_arguments::help);
  m_parser->add_description("MCMC diagnostics subcommand.");
  m_parser->add_subparser(*m_llhMapModule.get_parser());
  m_parser->add_subparser(*m_diagMCMCModule.get_parser());
  return m_parser.get();
}

int DiagModule::Run() {
  if (m_parser->is_subcommand_used("llhmap")) return m_llhMapModule.Run();
  if (m_parser->is_subcommand_used("mcmc")) return m_diagMCMCModule.Run();
  return 0;
}
} // namespace M3
