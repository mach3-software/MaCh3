#pragma once
#include "cli/api/plugin.hpp"

namespace M3{

  class DiagMCMCModule: public IModule{

    public:
      virtual ~DiagMCMCModule();
      MaCh3ArgumentParser* get_parser() override;
      int run() override;


    private:
      MaCh3ArgumentParser* m_parser;
  };
};
