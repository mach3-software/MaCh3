#pragma once
#include "nudock.hpp"
#include <Manager/Manager.h>

void InitialiseNuDockObj(manager *man, std::unique_ptr<NuDock> &nudock_ptr);

void FormatOscParsForNuDock(const std::string &param_name, double &param_value);

void FormatOscParsForMaCh3(const std::string &param_name, double &param_value);

extern const std::unordered_map<std::string, std::string> NuDockOscNameMap;
extern const std::unordered_map<std::string, std::string> NuDockOscNameMap_r;
