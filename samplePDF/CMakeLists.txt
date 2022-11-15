set(SamplePDF_implementation_files 
	samplePDFBase.cpp 
	samplePDFFDBase.cpp 
	interfacePDFEbE.cpp 
	Structs.cpp
        probGpu.cu )
        

set(SamplePDF_header_files 
	samplePDFBase.h
	samplePDFFDBase.h
	interfacePDFEbE.h
	Structs.h )

add_library(SamplePDF SHARED ${SamplePDF_implementation_files})

target_include_directories(SamplePDF PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(SamplePDF Splines Covariance Manager Prob3plusplus)

if(NOT ${CPU_ONLY})
  add_library(probGpu probGpu.cu)
  set_target_properties(SamplePDF PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_property(TARGET probGpu SamplePDF PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)
  target_link_libraries(SamplePDF probGpu)
else()
endif()

set_target_properties(SamplePDF PROPERTIES 
  PUBLIC_HEADER "${SamplePDF_header_files}"
  EXPORT_NAME SamplePDF)

install(TARGETS SamplePDF
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/samplePDF)


