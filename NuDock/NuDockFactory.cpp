#include "NuDockFactory.h"
#include "nudock.hpp"

const std::unordered_map<std::string, std::string> NuDockOscNameMap = {
  {"Theta12", "sin2th_12"},
  {"Theta13", "sin2th_13"},
  {"Theta23", "sin2th_23"},
  {"DeltaCP", "delta_cp"},
  {"Deltam2_21", "delm2_12"},
  {"Deltam2_32", "delm2_23"},
};

const std::unordered_map<std::string, std::string> NuDockOscNameMap_r = {
  {"sin2th_12", "Theta12"},
  {"sin2th_13", "Theta13"},
  {"sin2th_23", "Theta23"},
  {"delta_cp", "DeltaCP"},
  {"delm2_12", "Deltam2_21"},
  {"delm2_23", "Deltam2_32"},
};

const std::unordered_map<std::string, CommunicationType> CommunicationTypeMap = {
  {"LOCALHOST", CommunicationType::LOCALHOST},
  {"UNIX_DOMAIN_SOCKET", CommunicationType::UNIX_DOMAIN_SOCKET},
  {"TCP", CommunicationType::TCP}
};

const std::unordered_map<std::string, VerbosityLevel> VerbosityLevelMap = {
  {"DEBUG", VerbosityLevel::DEBUG},
  {"INFO", VerbosityLevel::INFO},
  {"WARNING", VerbosityLevel::WARNING},
  {"ERROR", VerbosityLevel::ERROR}
};

void InitialiseNuDockObj(manager *man,
                                std::unique_ptr<NuDock> &nudock_ptr) {
  // Here man can be either the sample manager or the fit manager, as long as it
  // has the NuDock config
  if (nudock_ptr) {
    MACH3LOG_INFO("NuDock object has already been created so I am not "
                  "re-initialising the object");
    return;
  }
  std::string nudock_conf_name = "NuDock";
  if (!man->raw()["NuDock"]) {
    if (!man->raw()["NuDockClient"]) {
      MACH3LOG_ERROR(
          "No NuDock configuration found in the provided manager object");
      throw MaCh3Exception(__FILE__, __LINE__);
    } else {
      MACH3LOG_INFO("Setting up NuDock client.");
      nudock_conf_name = "NuDockClient";
    }
  } else {
    MACH3LOG_INFO("Setting up NuDock server.");
  }

  bool useDebug = GetFromManager(man->raw()[nudock_conf_name]["Debug"], false);
  std::string schemaLocation = GetFromManager<std::string>(
      man->raw()[nudock_conf_name]["SchemaLocation"], "");
  std::string commTypeStr = GetFromManager<std::string>(
      man->raw()[nudock_conf_name]["CommunicationType"], "LOCALHOST");
  int port = GetFromManager<int>(man->raw()[nudock_conf_name]["Port"], 1234);
  std::string verbosity =
      GetFromManager<std::string>(man->raw()[nudock_conf_name]["NuDockVerbosity"], "INFO");

  CommunicationType commType;
  if (CommunicationTypeMap.find(commTypeStr) != CommunicationTypeMap.end()) {
    commType = CommunicationTypeMap.at(commTypeStr);
  } else {
    MACH3LOG_ERROR(
        "Unsupported communication type for NuDock: {}. Supported types "
        "are LOCALHOST, UNIX_DOMAIN_SOCKET, and TCP.",
        commTypeStr);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  VerbosityLevel verbLevel;
  if (VerbosityLevelMap.find(verbosity) != VerbosityLevelMap.end()) {
    verbLevel = VerbosityLevelMap.at(verbosity);
  } else {
    MACH3LOG_ERROR(
        "Unsupported verbosity level for NuDock: {}. Supported levels are "
        "DEBUG, INFO, WARN, and ERROR.",
        verbosity);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  nudock_ptr =
      std::make_unique<NuDock>(useDebug, schemaLocation, commType, port, verbLevel);
  MACH3LOG_INFO("NuDock object created with communication type: {} on port: {}", commTypeStr, port);
  MACH3LOG_INFO("NuDock schema location: {}", schemaLocation);
  MACH3LOG_INFO("NuDock debug mode: {}", useDebug);
  MACH3LOG_INFO("NuDock verbosity level: {}", verbosity);
}

void FormatOscParsForNuDock(const std::string &param_name, double &param_value) {
  if (param_name == "Theta12" || param_name == "Theta13" || param_name == "Theta23") {
    param_value = asin(sqrt(param_value));
  }
}

void FormatOscParsForMaCh3(const std::string &param_name, double &param_value) {
  if (param_name == "Theta12" || param_name == "Theta13" || param_name == "Theta23") {
    param_value = sin(param_value) * sin(param_value);
  }
}
