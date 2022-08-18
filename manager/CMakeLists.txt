set(manager_implementation_files
	manager.cpp)

set(manager_header_files
	manager.h)

add_library(Manager SHARED ${manager_implementation_files})
target_link_libraries(Manager yaml-cpp)

install(TARGETS Manager DESTINATION lib)

