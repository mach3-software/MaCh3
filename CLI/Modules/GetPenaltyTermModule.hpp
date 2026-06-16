#pragma once
#include "yaml-cpp/yaml.h"
#include "CLI/API/plugin.hpp"

namespace M3{

  class GetPenaltyTermModule: public IModule{

    public:
      virtual ~GetPenaltyTermModule();
      MaCh3ArgumentParser* get_parser() override;
      int run() override;

    private:
      void ReadCovFile(const std::string& inputFile,
                       std::vector <double>& Prior,
                       std::vector <bool>& isFlat,
                       std::vector<std::string>& ParamNames,
                       std::vector<std::vector<double>>& invCovMatrix,
                       int& nParams);
      void LoadSettings(YAML::Node& Settings,
                        std::vector<std::string>& SetsNames,
                        std::vector<std::string>& FancyTitle,
                        std::vector<std::vector<bool>>& isRelevantParam,
                        const std::vector<std::string>& ParamNames,
                        const int nParams);
      void GetPenaltyTerm(const std::string& inputFile, const std::string& configFile);

    private:
      MaCh3ArgumentParser* m_parser;
  };
};
