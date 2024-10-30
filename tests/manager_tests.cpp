#include "catch2/catch_test_macros.hpp"

#include "manager/YamlHelper.h"

TEST_CASE("OverrideConfig", "[Yamlhelper]") {
  YAML::Node lineup = YAML::Load("{1B: Prince Fielder, 2B: Rickie Weeks, LF: Ryan Braun}");

  REQUIRE(lineup.size() == 3);
  REQUIRE(lineup["1B"].as<std::string>() == "Prince Fielder");

  OverrideConfig(lineup, "1B", "Fielder formerly know as Prince");

  REQUIRE(lineup["1B"].as<std::string>() == "Fielder formerly know as Prince");

  OverrideConfig(lineup, "1B", 123);

  REQUIRE(lineup["1B"].as<unsigned>() == 123u);
}