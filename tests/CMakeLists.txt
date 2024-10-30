add_executable(manager_tests manager_tests.cpp)
target_link_libraries(manager_tests PRIVATE Catch2::Catch2WithMain Manager)

catch_discover_tests(manager_tests)
