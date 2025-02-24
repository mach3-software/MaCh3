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
    gpuFitterProcessorUtils.cuh
)

add_library(Fitter SHARED
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
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:gpuFitterProcessorUtils.cu>
)

if(MaCh3_GPU_ENABLED)
    set_target_properties(Fitter PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET Fitter PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

set_target_properties(Fitter PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Fitter)

target_link_libraries(Fitter PUBLIC Samples)
target_link_libraries(Fitter PRIVATE MaCh3Warnings)

target_include_directories(Fitter PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS Fitter
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/mcmc)

add_library(MaCh3::Fitter ALIAS Fitter)
