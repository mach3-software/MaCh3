set(manager_implementation_files
	manager.cpp)

set(manager_header_files
	manager.h)

add_library(Manager SHARED ${manager_implementation_files})
target_include_directories(Manager PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(Manager yaml-cpp)

set_target_properties(Manager PROPERTIES 
  PUBLIC_HEADER "${manager_header_files}"
  EXPORT_NAME Manager)

install(TARGETS Manager
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/manager)

