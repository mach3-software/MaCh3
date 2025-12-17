CPMAddPackage(
    NAME Catch2
    GITHUB_REPOSITORY catchorg/Catch2
    VERSION 3.3.2
)

LIST(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

include(CTest)
include(Catch)

add_executable(ParameterListTest ParameterListTest.cpp)
target_link_libraries(ParameterListTest PRIVATE Catch2::Catch2WithMain Parameters)
catch_discover_tests(ParameterListTest)

add_executable(StepProposerTest StepProposerTest.cpp)
target_link_libraries(StepProposerTest PRIVATE Catch2::Catch2WithMain Parameters)
catch_discover_tests(StepProposerTest)
