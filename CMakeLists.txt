set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.17 FATAL_ERROR)  

project(MaCh3 LANGUAGES CXX)

# creates a custom target to drive download, update/patch, configure, build, 
# install and test steps of an external project
include(ExternalProject)

# Changes default install path to be a subdirectory of the build dir.
# Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif()

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Use the compilers found in the path
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# To make output lines about environment vars
include(${CMAKE_SOURCE_DIR}/cmake/cmessage.cmake)

# Check the where to install 
cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")

# add cmake script files
if(NOT ${CMAKE_SCRIPT_SETUP})  #Check if setup by Experiment MaCh3
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

  #CUDA
  if(NOT ${CPU_ONLY})
    include(${CMAKE_SOURCE_DIR}/cmake/CUDASetup.cmake)
  endif()

  # ROOT  
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
  include(${CMAKE_SOURCE_DIR}/cmake/ROOTSetup.cmake)

endif()

#LW - Possibly use CPM to install dependencies...
#CPM
#set(CPM_DOWNLOAD_VERSION 0.34.0)

#if(CPM_SOURCE_CACHE)
#  set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
#elseif(DEFINED ENV{CPM_SOURCE_CACHE})
#  set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
#else()
#  set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
#endif()

#if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
#  message(STATUS "Downloading CPM.cmake to ${CPM_DOWNLOAD_LOCATION}")
#  file(DOWNLOAD
#       https://github.com/TheLartians/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
#       ${CPM_DOWNLOAD_LOCATION}
#  )
#endif()

#include(${CPM_DOWNLOAD_LOCATION})

##Prob3++ with CPM
#CPMFindPackage(
#   NAME prob3
#   URL "https://github.com/hyperk/Prob3.git"
#)


# Custom commands
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

#Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
add_compile_options(-I${CUDA_SAMPLES}/common/inc -L${ROOT_LIBDIR})
add_compile_options(-g -fopenmp -std=c++11 -O3 -Wextra)
add_compile_options(-I$ENV{MACH3CORE})
add_link_options(-fopenmp)


# Build components
#add_subdirectory(libconfig)
#add_subdirectory(Diagnostics)
add_subdirectory(Prob3++)
#add_subdirectory(CUDAProb3)
add_subdirectory(throwParms)
#add_subdirectory(manager)
add_subdirectory(covariance)
add_subdirectory(splines)
add_subdirectory(samplePDF)
#add_subdirectory(mcmc)
