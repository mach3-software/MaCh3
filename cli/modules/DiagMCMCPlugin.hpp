#pragma once
#include "api/plugin.hpp"

namespace M3{

  class DiagMCMCPlugin: public IPlugin{

    public:
      virtual ~DiagMCMCPlugin();
      MaCh3ArgumentParser* get_parser() override;
      int run() override;


    private:
      MaCh3ArgumentParser* m_parser;
  };
};
