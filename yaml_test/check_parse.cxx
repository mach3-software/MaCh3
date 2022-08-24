
#include "yaml-cpp/yaml.h"
#include "covarianceParserYAML.h"

#include <iostream>
#include "/scratch2/akaboth/MaCh3_YAML/build/_deps/eigen-src/Eigen/Dense"

int main(int argc, char const *argv[]) {

   auto doc = YAML::LoadFile("/scratch2/akaboth/MaCh3_YAML/yaml_test/xsec_test.yaml");
   
   auto keys = doc["keys"].as<std::vector<std::string>>();
   
   std::cout << "Keys: " << std::endl;
   for (auto k : keys) {
      std::cout << "\t" << k << std::endl;
   }
   std::cout << "\n\n" << std::endl;

   std::cout << doc["parameters"].size() << std::endl;

   std::cout << "Parameters: " << std::endl;
   for (auto const &param : doc["parameters"]) {
      std::cout << "\tname: " << param["name"].as<std::string>() << std::endl;
      if (param["correlation"]) {
	 std::cout << "\t\thas " << param["correlation"].size() << " correlations:" << std::endl;
	 for (auto const &corr : param["correlation"]) {
	    std::cout << "\t\t\t" << corr.first.as<std::string>() << ": " << corr.second.as<double>() << std::endl;
	 }
      }
   }

   covarianceParserYAML parser = covarianceParserYAML("/scratch2/akaboth/MaCh3_YAML/yaml_test/xsec_test.yaml");
   std::cout << parser.GetNParameters() << std::endl;
   parser.GetCovMatrix()->Print();

}
