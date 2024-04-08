#include "manager/YamlHelper.h"

// **********************
// Function to convert a YAML string to a YAML node
YAML::Node STRINGtoYAML(const std::string& yaml_string) {
// **********************

  try {
    return YAML::Load(yaml_string);
  } catch (const YAML::ParserException& e) {
    std::cerr << "Error parsing YAML string: " << e.what() << std::endl;
    return YAML::Node();
  }
}

// **********************
// Function to convert a YAML node to a YAML string
std::string YAMLtoSTRING(const YAML::Node& node) {
// **********************

  YAML::Emitter emitter;
  emitter << node;
  return emitter.c_str();
}

// **********************
//KS: Converts ROOT TMacro to string
std::string TMacroToString(const TMacro& macro) {
// **********************

  std::stringstream ss;

  // Retrieve lines from TMacro
  TList* linesList = macro.GetListOfLines();
  if (!linesList) {
    std::cerr << "Error: Failed to retrieve lines from TMacro." << std::endl;
    return "";
  }

  TIter nextLine(linesList);
  TObject *obj = nullptr;
  while ((obj = nextLine())) {
    TObjString* line = dynamic_cast<TObjString*>(obj);
    if (!line) {
      std::cerr << "Error: Failed to cast object to TObjString." << std::endl;
      continue;
    }
    ss << line->GetString() << std::endl;
  }

  return ss.str();
}

// **********************
//KS: Converts ROOT TMacro to yaml node
YAML::Node TMacroToYAML(const TMacro& macro)  {
// **********************

  std::string yaml_string = TMacroToString(macro);

  // Convert the YAML string to a YAML node
  YAML::Node yaml_node = STRINGtoYAML(yaml_string);

  return yaml_node;
}
