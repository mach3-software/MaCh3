set(HEADERS
    FitterBase.h
    mcmc.h
    LikelihoodFit.h
    tune.h
    stretch.h
    MCMCProcessor.h
)

add_library(MCMC SHARED
    FitterBase.cpp
    mcmc.cpp
    LikelihoodFit.cpp
    MinuitFit.cpp
    tune.cpp
    stretch.cpp
    MCMCProcessor.cpp
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:gpuMCMCProcessorUtils.cu>
)


if(NOT ${CPU_ONLY})
    set_target_properties(MCMC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET MCMC PROPERTY CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES})
endif()

set_target_properties(MCMC PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME MCMC)

target_link_libraries(MCMC Manager ROOT::ROOT)
target_include_directories(MCMC PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:$ENV{CUDAPATH}/include>>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:$ENV{CUDAPATH}/common/inc>>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:$ENV{CUDAPATH}/samples/common/inc>>
)

install(TARGETS MCMC
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/mcmc)

add_library(MaCh3::MCMC ALIAS MCMC)
