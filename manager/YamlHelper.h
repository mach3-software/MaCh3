#pragma once

#include <iostream>
#include <string>

#include "yaml-cpp/yaml.h"

template<typename Type>
Type GetFromManager(const YAML::Node& node, Type defval)
{
  YAML::Node tmpNode = node;

  if(!tmpNode){
    return defval;
  }
  return tmpNode.as<Type>();
}

// Use this like this CheckNodeExists(config, "LikelihoodOptions", "TestStatistic");
//KS: Base case for recursion
template<typename T>
bool CheckNodeExistsHelper(const T& node) {
  (void)node;
  return true; // Node exists, return true
}

//KS: Recursive function to traverse YAML nodes
template<typename T, typename... Args>
bool CheckNodeExistsHelper(const T& node, const std::string& key, Args... args) {
  if (!node[key]) {
    //std::cerr << "Node " << key << " doesn't exist." << std::endl;
    return false;
  }
  return CheckNodeExistsHelper(node[key], args...);
}

//KS: Wrapper function to call the recursive helper
template<typename... Args>
bool CheckNodeExists(const YAML::Node& node, Args... args) {
  return CheckNodeExistsHelper(node, args...);
}



/// Use this like this FindFromManager<std::string>(config, "LikelihoodOptions", "TestStatistic");

// Base case for recursion
template<typename T>
T FindFromManagerHelper(const YAML::Node& node) {
  return node.as<T>(); // Convert YAML node to the specified type
}

// Recursive function to traverse YAML nodes
template<typename T, typename... Args>
T FindFromManagerHelper(const YAML::Node& node, const std::string& key, Args... args) {
  if (!node[key]) {
    std::cerr << "Node " << key << " doesn't exist." << std::endl;
    throw;
    return T();
  }
  return FindFromManagerHelper<T>(node[key], args...); // Recursive call
}

// Wrapper function to call the recursive helper
template<typename T, typename... Args>
T FindFromManager(const YAML::Node& node, Args... args) {
  return FindFromManagerHelper<T>(node, args...);
}

