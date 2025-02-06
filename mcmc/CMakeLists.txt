set(HEADERS
    FitterBase.h
    mcmc.h
    LikelihoodFit.h
    MinuitFit.h
    PSO.h
    MCMCProcessor.h
    SampleSummary.h
    MaCh3Factory.h
    StatisticalUtils.h
    gpuMCMCProcessorUtils.cuh
)

add_library(MCMC SHARED
    FitterBase.cpp
    mcmc.cpp
    LikelihoodFit.cpp
    $<$<BOOL:${MaCh3_MINUIT2_ENABLED}>:MinuitFit.cpp>
    PSO.cpp
    MCMCProcessor.cpp
    SampleSummary.cpp
    MaCh3Factory.cpp
    StatisticalUtils.cpp
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:gpuMCMCProcessorUtils.cu>
)

if(NOT CPU_ONLY)
    set_target_properties(MCMC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET MCMC PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

set_target_properties(MCMC PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME MCMC)

target_link_libraries(MCMC PUBLIC SamplePDF)
target_link_libraries(MCMC PRIVATE MaCh3Warnings)

target_include_directories(MCMC PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS MCMC
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/mcmc)

add_library(MaCh3::MCMC ALIAS MCMC)
