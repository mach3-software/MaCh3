#pragma once
#include "api/plugin.hpp"

namespace mach3{

  class DiagMCMCPlugin: public IPlugin{

    public:
      virtual ~DiagMCMCPlugin(){
        if (m_parser) { delete m_parser; } 
      }
      MaCh3ArgumentParser* get_parser() override;
      int run() override;


    private:
      MaCh3ArgumentParser* m_parser;
  };
};
