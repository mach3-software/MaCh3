include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG d924b84eba1f0f0adf38b20b7b4829f6f65b6570
)
FetchContent_MakeAvailable(argparse)

set(HEADERS
    plugin.hpp
    m3argparse.hpp
)
list(TRANSFORM HEADERS PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

install(DIRECTORY ${argparse_SOURCE_DIR}/include/
        DESTINATION include/CLI/API/argparse)
target_include_directories(argparse INTERFACE
    $<BUILD_INTERFACE:${argparse_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/CLI/API/argparse>
)
install(TARGETS argparse
        EXPORT MaCh3-targets
        INCLUDES DESTINATION include/CLI/API/argparse)


add_library(CLIApi INTERFACE)

target_link_libraries(CLIApi INTERFACE argparse)
target_include_directories(CLIApi INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/CLI/API>
)

set_target_properties(CLIApi PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME CLIApi)


install(TARGETS CLIApi
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/CLI/API)

#install(FILES ${HEADERS} DESTINATION include/api)
add_library(MaCh3::CLIApi ALIAS CLIApi)