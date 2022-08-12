include(ExternalProject)

execute_process(COMMAND make -j4 shared
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE LIBPROB3_FILES ${CMAKE_CURRENT_SOURCE_DIR}/../lib/libThreeProb*.so*)

cmessage ( STATUS "[PROB3++]: files : ${LIBPROB3_FILES} ")

foreach(a ${LIBPROB3_FILES})
	  file(COPY ${a} DESTINATION ${CMAKE_BINARY_DIR}/lib)
endforeach()

