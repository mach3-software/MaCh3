set(OscClass_implementation_files 
	OscClass_CUDAProb3.cpp 
	)

set(OscClass_header_files 
	OscClass_CUDAProb3.h
	)

add_library(OscClass SHARED ${OscClass_implementation_files})
target_include_directories(OscClass PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

if (USE_PROB3)
   target_link_libraries(OscClass Manager CUDAProb3Beam CUDAProb3Atmos)
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


