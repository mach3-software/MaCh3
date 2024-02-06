set(OscClass_implementation_files 
	OscClass_CUDAProb3.cpp 
	)

set(OscClass_header_files 
	OscClass_CUDAProb3.h
	)

add_library(OscClass SHARED ${OscClass_implementation_files})
target_include_directories(OscClass PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

if (USE_PROB3)
   message(FATAL_ERROR "OscClass only defined using CUDAProb3. USE_PROB3 is defined and indicates Prob3++/ProbGPU is being expected to be used.")
else()
   target_link_libraries(OscClass Manager CUDAProb3Beam CUDAProb3Atmos)
endif()

set_target_properties(OscClass PROPERTIES
  PUBLIC_HEADER "${OscClass_header_files}"
  EXPORT_NAME OscClass)

install(TARGETS OscClass
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/OscClass)


