set(manager_implementation_files
	manager.cpp)

set(manager_header_files
	manager.h)

add_library(Manager SHARED ${manager_implementation_files})
target_include_directories(Manager PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(Manager)

install(TARGETS Manager DESTINATION lib)

