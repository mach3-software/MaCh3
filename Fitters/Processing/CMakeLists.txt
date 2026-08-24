set(HEADERS
    MCMCProcessor.h
    OscProcessor.h
    RHatCalculator.h
    StatisticalUtils.h
    gpuMCMCProcessorUtils.cuh
)

add_library(FitProcessing SHARED
    MCMCProcessor.cpp
    OscProcessor.cpp
    RHatCalculator.cpp
    StatisticalUtils.cpp
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:gpuMCMCProcessorUtils.cu>
)

if(MaCh3_GPU_ENABLED)
    set_target_properties(FitProcessing PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET FitProcessing PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

set_target_properties(FitProcessing PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME FitProcessing)

target_link_libraries(FitProcessing PUBLIC Samples ROOT::ROOTDataFrame)
target_link_libraries(FitProcessing PRIVATE MaCh3Warnings)

target_include_directories(FitProcessing PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS FitProcessing
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/FitProcessing)

add_library(MaCh3::FitProcessing ALIAS FitProcessing)
target_link_libraries(MaCh3 INTERFACE FitProcessing)
