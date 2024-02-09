set(HEADERS
    samplePDFBase.h
    samplePDFFDBase.h
    interfacePDFEbE.h
    Structs.h
)

add_library(SamplePDF SHARED
	samplePDFBase.cpp
	samplePDFFDBase.cpp
	interfacePDFEbE.cpp
	Structs.cpp
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:probGpu.cu>
)

target_link_libraries(SamplePDF Splines Covariance Manager OscClass
    $<$<BOOL:${USE_PROB3}>:Prob3plusplus>
    $<$<NOT:$<BOOL:${USE_PROB3}>>:CUDAProb3Beam CUDAProb3Atmos>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:probGpu>
)
target_include_directories(SamplePDF PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

if(NOT ${CPU_ONLY})
    set_target_properties(SamplePDF PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET SamplePDF PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)
    target_compile_definitions(SamplePDF PUBLIC GPU_ON)
endif()

set_target_properties(SamplePDF PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME SamplePDF)

install(TARGETS SamplePDF
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/samplePDF)

add_library(MaCh3::SamplePDF ALIAS SamplePDF)


