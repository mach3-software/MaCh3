set(HEADERS manager.h

)

add_library(Manager SHARED
    manager.cpp
)

set_target_properties(Manager PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Manager)


target_link_libraries(Manager yaml-cpp)
target_include_directories(Manager PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS yaml-cpp
        EXPORT mach3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/yaml-cpp
        )

install(TARGETS Manager
        EXPORT mach3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/manager)

add_library(MACH3::Manager ALIAS Manager)


