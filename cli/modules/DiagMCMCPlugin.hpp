#pragma once
#include "cli/api/plugin.hpp"

namespace M3{

  class DiagMCMCPlugin: public IModule{

    public:
      virtual ~DiagMCMCPlugin();
      MaCh3ArgumentParser* get_parser() override;
      int run() override;


    private:
      MaCh3ArgumentParser* m_parser;
  };
};
