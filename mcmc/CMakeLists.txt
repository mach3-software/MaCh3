set(mcmc_implementation_files
    FitterBase.cpp
    mcmc.cpp
    LikelihoodFit.cpp
    MinuitFit.cpp
    tune.cpp
    stretch.cpp
    MCMCProcessor.cpp
    gpuMCMCProcessorUtils.cu)

set(mcmc_implementation_files_CPU 
    FitterBase.cpp
    mcmc.cpp
    LikelihoodFit.cpp
    MinuitFit.cpp
    tune.cpp
    stretch.cpp
    MCMCProcessor.cpp)

set(mcmc_header_files
    FitterBase.h
    mcmc.h
    LikelihoodFit.h
    MinuitFit.h
    tune.h
    stretch.h
    MCMCProcessor.h)

if(NOT ${CPU_ONLY})
  add_library(MCMC SHARED ${mcmc_implementation_files})
  target_include_directories(MCMC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
  add_library(gpuMCMC gpuMCMCProcessorUtils.cu)
  set_target_properties(MCMC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_property(TARGET gpuMCMCProcessorUtils MCMC PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)
  target_link_libraries(gpuMCMCProcessorUtils)
else()
    add_library(MCMC SHARED ${mcmc_implementation_files_CPU})
    target_include_directories(MCMC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
	target_link_libraries(MCMC Manager ROOT::ROOT)
endif()

set_target_properties(MCMC PROPERTIES 
  PUBLIC_HEADER "${mcmc_header_files}"
  EXPORT_NAME MCMC)

install(TARGETS MCMC
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/mcmc)
