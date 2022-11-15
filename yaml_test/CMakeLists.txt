add_executable(check_parse check_parse.cxx covarianceParserYAML.cpp)
target_include_directories(check_parse PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(check_parse yaml-cpp ROOT::ROOT)

install(TARGETS check_parse DESTINATION bin)
