#pragma once

// C++ Includes
#include <iostream>
#include <fstream>
#include <string>

// ROOT Includes
#include "TMacro.h"
#include "TList.h"
#include "TObjString.h"

// yaml Includes
#include "yaml-cpp/yaml.h"

// **********************
/// @brief Get content of config file if node is not found take default value specified
/// @param node Yaml node
/// @param defval Default value which will be used in case node doesn't exist
template<typename Type>
Type GetFromManager(const YAML::Node& node, Type defval) {
// **********************
  YAML::Node tmpNode = node;

  if(!tmpNode){
    return defval;
  }
  return tmpNode.as<Type>();
}

// **********************
/// @brief Use this like this CheckNodeExists(config, "LikelihoodOptions", "TestStatistic");
/// KS: Base case for recursion
template<typename T>
bool CheckNodeExistsHelper(const T& node) {
// **********************
  (void)node;
  return true; // Node exists, return true
}

// **********************
/// @brief KS: Recursive function to traverse YAML nodes
template<typename T, typename... Args>
bool CheckNodeExistsHelper(const T& node, const std::string& key, Args... args) {\
// **********************
  if (!node[key]) {
    //std::cerr << "Node " << key << " doesn't exist." << std::endl;
    return false;
  }
  return CheckNodeExistsHelper(node[key], args...);
}

// **********************
/// @brief KS: Wrapper function to call the recursive helper
template<typename... Args>
bool CheckNodeExists(const YAML::Node& node, Args... args) {
// **********************
  return CheckNodeExistsHelper(node, args...);
}

// **********************
/// Use this like this FindFromManager<std::string>(config, "LikelihoodOptions", "TestStatistic");
/// Base case for recursion
template<typename T>
T FindFromManagerHelper(const YAML::Node& node) {
// **********************
  return node.as<T>(); // Convert YAML node to the specified type
}

// **********************
/// @brief Recursive function to traverse YAML nodes
template<typename T, typename... Args>
T FindFromManagerHelper(const YAML::Node& node, const std::string& key, Args... args) {
  // **********************
  if (!node[key]) {
    std::cerr << "Node " << key << " doesn't exist." << std::endl;
    throw;
    return T();
  }
  return FindFromManagerHelper<T>(node[key], args...); // Recursive call
}
// **********************
/// @brief Wrapper function to call the recursive helper
template<typename T, typename... Args>
T FindFromManager(const YAML::Node& node, Args... args) {
// **********************
  return FindFromManagerHelper<T>(node, args...);
}

// **********************
/// @brief Function to convert a YAML string to a YAML node
/// @param yaml_string String which will be converted to yaml node
inline YAML::Node STRINGtoYAML(const std::string& yaml_string){
// **********************
  try {
    return YAML::Load(yaml_string);
  } catch (const YAML::ParserException& e) {
    std::cerr << "Error parsing YAML string: " << e.what() << std::endl;
    return YAML::Node();
  }
}

// **********************
/// @brief KS: Convert a YAML node to a string representation.
/// @param node The YAML node to convert to a string.
/// @return std::string The string representation of the YAML node.
inline std::string YAMLtoSTRING(const YAML::Node& node) {
// **********************
  YAML::Emitter emitter;
  emitter << node;
  return emitter.c_str();
}

// **********************
/// @brief KS: Convert a ROOT TMacro object to a string representation.
/// @param macro The ROOT TMacro object to convert.
/// @return std::string The string representation of the TMacro object.
inline std::string TMacroToString(const TMacro& macro) {
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
/// @brief KS: Convert a ROOT TMacro object to a YAML node.
/// @param macro The ROOT TMacro object to convert.
/// @return YAML::Node The YAML node representing the TMacro object.
inline YAML::Node TMacroToYAML(const TMacro& macro) {
// **********************
  std::string yaml_string = TMacroToString(macro);

  // Convert the YAML string to a YAML node
  YAML::Node yaml_node = STRINGtoYAML(yaml_string);

  return yaml_node;
}

// **********************
/// @brief Convert a YAML node to a ROOT TMacro object.
/// @param yaml_node The YAML node to convert to a TMacro.
/// @param name Name of TMacro that will be saved
/// @return TMacro The TMacro object constructed from the YAML node.
inline TMacro YAMLtoTMacro(const YAML::Node& yaml_node, const std::string& name) {
// **********************
  // Convert the YAML node to a string representation
  std::string macro_string = YAMLtoSTRING(yaml_node);

  // Create a TMacro object with the collected lines
  TMacro macro(name.c_str(), name.c_str());
  macro.AddLine(macro_string.c_str());

  return macro;
}

// **********************
/// @brief Compare if yaml nodes are identical
/// @param node1 The first YAML node to compare.
/// @param node2 The second YAML node to compare.
/// @return true If the two nodes are equivalent in type and content.
/// @return false If the two nodes differ in structure or content.
inline bool compareYAMLNodes(const YAML::Node& node1, const YAML::Node& node2) {
// **********************
  // Check if the types of the nodes match
  if (node1.Type() != node2.Type()) {
    return false;
  }

  // Compare scalar types (like strings, numbers)
  if (node1.IsScalar() && node2.IsScalar()) {
    return node1.as<std::string>() == node2.as<std::string>();
  }

  // Compare sequences (like YAML lists)
  if (node1.IsSequence() && node2.IsSequence()) {
    if (node1.size() != node2.size()) {
      return false;
    }
    for (std::size_t i = 0; i < node1.size(); ++i) {
      if (!compareYAMLNodes(node1[i], node2[i])) {
        return false;
      }
    }
    return true;
  }

  // Compare maps (like YAML dictionaries)
  if (node1.IsMap() && node2.IsMap()) {
    if (node1.size() != node2.size()) {
      return false;
    }
    for (auto it1 = node1.begin(); it1 != node1.end(); ++it1) {
      auto key = it1->first.as<std::string>();
      if (!node2[key] || !compareYAMLNodes(it1->second, node2[key])) {
        return false;
      }
    }
    return true;
  }

  // Default case: if it's neither scalar, sequence, nor map, consider it unequal
  return false;
}


// **********************
/// @brief Overrides the configuration settings based on provided arguments.
///
/// This function allows you to set configuration options in a nested YAML node.
/// @param node YAML node that will be modified
/// @param args The arguments to override the configuration. The last argument
///             will be used as the value
///
/// @note Example usage:
/// @code
/// OverrideConfig(config, "General", "OutputFile", "Wooimbouttamakeanameformyselfere.root");
/// OverrideConfig(config, "General", "MyDouble", 5.3);
/// @endcode
template <typename TValue>
void OverrideConfig(YAML::Node &node, std::string const &key, TValue val) {
// **********************
  node[key] = val;
}
template <typename... Args>
void OverrideConfig(YAML::Node &node, std::string const &key, Args... args) {
// **********************
  OverrideConfig(node[key], args...);
}
