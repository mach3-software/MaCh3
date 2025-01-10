#pragma once

// C++ Includes
#include <iostream>
#include <fstream>
#include <string>
#include <cxxabi.h>

// ROOT Includes
#include "TMacro.h"
#include "TList.h"
#include "TObjString.h"

// yaml Includes
#include "yaml-cpp/yaml.h"

// MaCh3 Includes
#include "manager/MaCh3Exception.h"

/// @file YamlHelper.h
/// @brief Utility functions for handling YAML nodes
/// @author Kamil Skwarczynski
/// @author Luke Pickering

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
    MACH3LOG_ERROR("Node {} doesn't exist.", key);
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
    MACH3LOG_ERROR("Error parsing YAML string: {}", e.what());
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
    MACH3LOG_ERROR("Failed to retrieve lines from TMacro.");
    return "";
  }

  TIter nextLine(linesList);
  TObject *obj = nullptr;
  while ((obj = nextLine())) {
    TObjString* line = dynamic_cast<TObjString*>(obj);
    if (!line) {
      MACH3LOG_ERROR("Failed to cast object to TObjString.");
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
void OverrideConfig(YAML::Node node, std::string const &key, TValue val) {
// **********************
  node[key] = val;
}
template <typename... Args>
void OverrideConfig(YAML::Node node, std::string const &key, Args... args) {
// **********************
  OverrideConfig(node[key], args...);
}

// **********************
/// @brief Function to demangle type names
inline std::string DemangleTypeName(const std::string& mangledName) {
// **********************
  int status = 0;
  char* demangledName = abi::__cxa_demangle(mangledName.c_str(), nullptr, nullptr, &status);
  std::string result = (status == 0) ? demangledName : mangledName;
  free(demangledName);
  return result;
}

// **********************
/// @brief Get content of config file
/// @param node Yaml node
template<typename Type>
Type Get(const YAML::Node& node, const std::string File, const int Line) {
// **********************
  if (!node) {
    MACH3LOG_ERROR("Empty Yaml node");
    throw MaCh3Exception(File , Line );
  }
  try {
    // Attempt to convert the node to the expected type
    return node.as<Type>();
  } catch (const YAML::BadConversion& e) {
    const std::string nodeAsString = YAMLtoSTRING(node);
    MACH3LOG_ERROR("YAML type mismatch: {}", e.what());
    MACH3LOG_ERROR("While trying to access variable {}", nodeAsString);
    throw MaCh3Exception(File , Line );
  }
}

// **********************
/// @brief Get content of config file if node is not found take default value specified
/// @param node Yaml node
/// @param defval Default value which will be used in case node doesn't exist
template<typename Type>
Type GetFromManager(const YAML::Node& node, Type defval, const std::string File = "", const int Line = 1) {
// **********************
  if (!node) {
    return defval;
  }
  try {
    // Attempt to convert the node to the expected type
    return node.as<Type>();
  } catch (const YAML::BadConversion& e) {
    const std::string nodeAsString = YAMLtoSTRING(node);
    MACH3LOG_ERROR("YAML type mismatch: {}", e.what());
    MACH3LOG_ERROR("While trying to access variable {}", nodeAsString);
    //const std::string expectedType = DemangleTypeName(typeid(Type).name());
    //MACH3LOG_ERROR("Expected argument is {}", expectedType);
    if(File == "") {
      throw MaCh3Exception(__FILE__ , __LINE__);
    } else {
      throw MaCh3Exception(File , Line );
    }
  }
}

// **********************
/// @brief Open YAML file
/// @param filename name of filename to open
/// @param File name of file where function is called
/// @param Line number where function is called
inline YAML::Node LoadYamlConfig(const std::string& filename, const std::string& File, const int Line) {
// **********************
  try {
    return YAML::LoadFile(filename);
  } catch (const std::exception& e) {
    MACH3LOG_ERROR("{}", e.what());
    MACH3LOG_ERROR("Can't open file {}", filename);
    throw MaCh3Exception(File, Line);
  }
}

/// Macro to simplify calling LoadYaml with file and line info
#define M3OpenConfig(filename) LoadYamlConfig((filename), __FILE__, __LINE__)
