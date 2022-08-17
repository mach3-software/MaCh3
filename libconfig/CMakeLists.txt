include(ExternalProject)

execute_process(COMMAND ./configure --disable-examples 
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
execute_process(COMMAND make 
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#file(GLOB_RECURSE LIBCONFIG_FILES lib/.libs/lib*.so*)
file(GLOB_RECURSE LIBCONFIG_FILES lib/.libs/libconfig++.so.9.1.3)

cmessage ( STATUS "[LIBCONFIG]: files : ${LIBCONFIG_FILES} ")

foreach(a ${LIBCONFIG_FILES})
	 get_filename_component(barename ${a} NAME)
         file(COPY ${a} DESTINATION ${CMAKE_BINARY_DIR}/lib/)
	 execute_process(COMMAND ln -sf ${barename} ${CMAKE_BINARY_DIR}/lib/libconfig++.so.9.1.3
	 	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#	  install(FILES ${a} DESTINATION lib)
endforeach()

