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

//KS: Base case for recursion
template<typename T>
bool checkPathExistsHelper(const T& node) {
  (void)node;
  return true; // Node exists, return true
}

//KS: Recursive function to traverse YAML nodes
template<typename T, typename... Args>
bool checkPathExistsHelper(const T& node, const std::string& key, Args... args) {
  if (!node[key]) {
    //std::cerr << "Node " << key << " doesn't exist." << std::endl;
    return false;
  }
  return checkPathExistsHelper(node[key], args...);
}

//KS: Wrapper function to call the recursive helper
template<typename... Args>
bool checkPathExists(const YAML::Node& node, Args... args) {
  return checkPathExistsHelper(node, args...);
}

//KS: WARNING Don't use this. It works but tmpNode.as<Type>() will delete YAML node... need to figure out how to fix this
template<typename Type, typename... Args>
Type FindFromManager(const YAML::Node& node, Args... args) {
  YAML::Node tmpNode = node;
  ((tmpNode = tmpNode[args]), ...);
  return tmpNode.as<Type>();
}
