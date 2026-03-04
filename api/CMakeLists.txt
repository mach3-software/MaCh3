include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

set(HEADERS
    plugin.hpp
)


add_library(api INTERFACE)

target_include_directories(api INTERFACE ${argparse_SOURCE_DIR}/include)

target_include_directories(api
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/api>
)


set_target_properties(api PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME api)


install(TARGETS api
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/api)

#install(FILES ${HEADERS} DESTINATION include/api)
add_library(MaCh3::api ALIAS api)