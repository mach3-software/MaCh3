add_executable(check_parse check_parse.cxx)
target_link_libraries(check_parse yaml-cpp)

install(TARGETS check_parse DESTINATION bin)