set(HEADERS
    FitterBase.h
    mcmc.h
    LikelihoodFit.h
    MinuitFit.h
    PSO.h
    MCMCProcessor.h
    OscProcessor.h
    SampleSummary.h
    MaCh3Factory.h
    StatisticalUtils.h
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:gpuMCMCProcessorUtils.cuh>
)

add_library(Fitters SHARED
    FitterBase.cpp
    mcmc.cpp
    LikelihoodFit.cpp
    $<$<BOOL:${MaCh3_MINUIT2_ENABLED}>:MinuitFit.cpp>
    PSO.cpp
    MCMCProcessor.cpp
    OscProcessor.cpp
    SampleSummary.cpp
    MaCh3Factory.cpp
    StatisticalUtils.cpp
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:gpuMCMCProcessorUtils.cu>
)

if(MaCh3_GPU_ENABLED)
    set_target_properties(Fitters PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET Fitters PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

set_target_properties(Fitters PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Fitters)

target_link_libraries(Fitters PUBLIC Samples)
target_link_libraries(Fitters PRIVATE MaCh3Warnings)

target_include_directories(Fitters PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS Fitters
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Fitters)

add_library(MaCh3::Fitters ALIAS Fitters)
