set(HEADERS
    SampleHandlerBase.h
    SampleHandlerFD.h
    Structs.h
    HistogramUtils.h
    FarDetectorCoreInfoStruct.h
)

add_library(SampleHandler SHARED
    SampleHandlerBase.cpp
    SampleHandlerFD.cpp
    HistogramUtils.cpp
)

target_link_libraries(SampleHandler PUBLIC Splines NuOscillator)
target_link_libraries(SampleHandler PRIVATE MaCh3Warnings)

target_include_directories(SampleHandler PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

if(MaCh3_GPU_ENABLED)
    set_target_properties(SampleHandler PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET SampleHandler PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

set_target_properties(SampleHandler PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME SampleHandler)

install(TARGETS SampleHandler
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/SampleHandler)

add_library(MaCh3::SampleHandler ALIAS SampleHandler)
