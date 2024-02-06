set(SamplePDF_implementation_files 
	samplePDFBase.cpp 
	samplePDFFDBase.cpp 
	interfacePDFEbE.cpp 
	Structs.cpp)

set(SamplePDF_header_files 
	samplePDFBase.h
	samplePDFFDBase.h
	interfacePDFEbE.h
	Structs.h )

if(NOT CPU_ONLY)
  list(APPEND SamplePDF_implementation_files probGpu.cu)
endif()

add_library(SamplePDF SHARED ${SamplePDF_implementation_files})
target_include_directories(SamplePDF PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

if (USE_PROB3)
  target_link_libraries(SamplePDF Splines Covariance Manager OscClass Prob3plusplus)
else()
  target_link_libraries(SamplePDF Splines Covariance Manager OscClass CUDAProb3Beam CUDAProb3Atmos)
endif()

if(NOT CPU_ONLY)
  add_library(probGpu probGpu.cu)
  set_target_properties(SamplePDF PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_property(TARGET probGpu SamplePDF PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)
  target_link_libraries(SamplePDF probGpu)
  target_compile_definitions(SamplePDF PUBLIC GPU_ON)
endif()

set_target_properties(SamplePDF PROPERTIES 
  PUBLIC_HEADER "${SamplePDF_header_files}"
  EXPORT_NAME SamplePDF)

install(TARGETS SamplePDF
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/samplePDF)
