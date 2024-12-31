set(HEADERS
    samplePDFBase.h
    samplePDFFDBase.h
    Structs.h
    HistogramUtils.h
    FarDetectorCoreInfoStruct.h
)

add_library(SamplePDF SHARED
    samplePDFBase.cpp
    samplePDFFDBase.cpp
    HistogramUtils.cpp
)

target_link_libraries(SamplePDF PUBLIC Splines NuOscillator Covariance)
target_link_libraries(SamplePDF PRIVATE MaCh3Warnings)

target_include_directories(SamplePDF PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

if(NOT CPU_ONLY)
    set_target_properties(SamplePDF PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET SamplePDF PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

set_target_properties(SamplePDF PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME SamplePDF)

install(TARGETS SamplePDF
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/samplePDF)

add_library(MaCh3::SamplePDF ALIAS SamplePDF)
