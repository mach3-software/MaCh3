include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

set(HEADERS
    plugin.hpp
    argparse.hpp
)


add_library(CLIApi INTERFACE)

target_include_directories(CLIApi INTERFACE ${argparse_SOURCE_DIR}/include)

target_include_directories(CLIApi
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/cli/api>
)

set_target_properties(CLIApi PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME CLIApi)


install(TARGETS CLIApi
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/cli/api)

#install(FILES ${HEADERS} DESTINATION include/api)
add_library(MaCh3::CLIApi ALIAS CLIApi)