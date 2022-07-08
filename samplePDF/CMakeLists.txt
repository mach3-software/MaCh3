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

add_library(probGpu probGpu.cu)

target_include_directories(SamplePDF PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_compile_features(SamplePDF PUBLIC)

set_target_properties(SamplePDF PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

set_property(TARGET probGpu SamplePDF PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)

add_dependencies(SamplePDF Splines Covariance probGpu)

install(TARGETS SamplePDF DESTINATION lib)
