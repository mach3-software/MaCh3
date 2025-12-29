#include "nudock.hpp"
#include <Manager/Manager.h>

inline void InitialiseNuDockObj(manager *man,
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
  CommunicationType commType;
  if (commTypeStr == "LOCALHOST") {
    commType = CommunicationType::LOCALHOST;
  } else if (commTypeStr == "UNIX_DOMAIN_SOCKET") {
    commType = CommunicationType::UNIX_DOMAIN_SOCKET;
  } else if (commTypeStr == "TCP") {
    commType = CommunicationType::TCP;
  } else {
    MACH3LOG_ERROR(
        "Unsupported communication type for NuDock: {}. Supported types "
        "are LOCALHOST, UNIX_DOMAIN_SOCKET, and TCP.",
        commTypeStr);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  int port = GetFromManager<int>(man->raw()[nudock_conf_name]["Port"], 1234);
  nudock_ptr =
      std::make_unique<NuDock>(useDebug, schemaLocation, commType, port);
  MACH3LOG_INFO("NuDock object created with communication type: {} on port: {}", commTypeStr, port);
  MACH3LOG_INFO("NuDock schema location: {}", schemaLocation);
  MACH3LOG_INFO("NuDock debug mode: {}", useDebug);
}