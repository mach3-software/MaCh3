set(HEADERS
    NuDockServer.h
    SampleHandlerNuDock.h
    NuDockFactory.h
)

add_library(MaCh3NuDock SHARED
    NuDockServer.cpp
    SampleHandlerNuDock.cpp
    NuDockFactory.cpp
)

set_target_properties(MaCh3NuDock PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME MaCh3NuDock)

# nlohmann_json is found at the parent CMakeLists level
# Create an alias if needed for consistency
if(NOT TARGET nlohmann_json::nlohmann_json)
    add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
endif()

CPMAddPackage(
    NAME NuDock
    GITHUB_REPOSITORY "NuDock/nudock"
    GIT_TAG feature/verbosity
    GIT_SHALLOW NO
    GIT_SUBMODULES_RECURSE ON
    # Apply patch to fix CMakeLists.txt until merged upstream
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy 
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/patches/NuDock_CMakeLists.txt
        <SOURCE_DIR>/CMakeLists.txt
)

if(TARGET NuDock::nudock)
    set(NUDOCK_TARGET NuDock::nudock)
elseif(TARGET nudock)
    # Create an alias if the namespaced target doesn't exist
    add_library(NuDock::nudock ALIAS nudock)
    set(NUDOCK_TARGET NuDock::nudock)
else()
    cmessage(FATAL_ERROR "MaCh3 Expected dependency target: NuDock::nudock or nudock")
endif()

target_link_libraries(MaCh3NuDock PUBLIC MaCh3::All MaCh3CompilerOptions ${NUDOCK_TARGET}) 

target_include_directories(MaCh3NuDock PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<BUILD_INTERFACE:${NuDock_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/cpp-httplib-src>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/nudock-src/externals/cpp-httplib>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/nudock-src/externals/json-schema-validator/src>
  $<INSTALL_INTERFACE:include>)

install(TARGETS nudock MaCh3NuDock
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/MaCh3NuDock)

add_library(MaCh3::MaCh3NuDock ALIAS MaCh3NuDock)