set(HEADERS
    SampleHandlerBase.h
    SampleHandlerFD.h
    Structs.h
    HistogramUtils.h
    FarDetectorCoreInfoStruct.h
    OscillationHandler.h
)

add_library(Samples SHARED
    SampleHandlerBase.cpp
    SampleHandlerFD.cpp
    HistogramUtils.cpp
    OscillationHandler.cpp
)

target_link_libraries(Samples PUBLIC Splines NuOscillator)
target_link_libraries(Samples PRIVATE MaCh3Warnings)

target_include_directories(Samples PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

if(MaCh3_GPU_ENABLED)
    set_target_properties(Samples PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET Samples PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

set_target_properties(Samples PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Samples)

install(TARGETS Samples
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Samples)

add_library(MaCh3::Samples ALIAS Samples)
