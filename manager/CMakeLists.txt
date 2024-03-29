set(HEADERS
    manager.h
    YamlHelper.h
)

add_library(Manager SHARED
    manager.cpp
)

set_target_properties(Manager PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Manager)

target_link_libraries(Manager yaml-cpp MaCh3CompilerOptions)
target_include_directories(Manager PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Manager
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/manager)

add_library(MaCh3::Manager ALIAS Manager)


install(TARGETS yaml-cpp
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/yaml-cpp
        )

install(DIRECTORY ${yaml-cpp_SOURCE_DIR}/include/
    DESTINATION include/
)
