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

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

include(CPM)

# Luke's helper cmake modules for the moment
CPMFindPackage(
    NAME CMakeModules
    GIT_TAG v0.0.11
    GITHUB_REPOSITORY NuHepMC/CMakeModules
    DOWNLOAD_ONLY
)
include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)

CPMFindPackage(
    NAME yaml-cpp
    VERSION 0.7.0
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG "yaml-cpp-0.7.0"
    OPTIONS
      "YAML_BUILD_SHARED_LIBS ON"
)

CPMFindPackage(
    NAME prob3pp
    VERSION 3.10.1
    GITHUB_REPOSITORY "mach3-software/Prob3plusplus"
    GIT_TAG v3.10.1
)

include(CMessage)

include(ROOT)

set(CMAKE_CXX_STANDARD 17)

# If ROOT sets the cxx standard, listen to it
if(DEFINED ROOT_CXX_STANDARD)
  if(NOT DEFINED CMAKE_CXX_STANDARD OR ROOT_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
  endif()
endif()

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

# Custom commands
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )


#If Single thread is chosen MultiThread flags are omitted
if(NOT ${SINGLE_THREAD_ONLY})

  #Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
  add_compile_options(-I${CUDA_SAMPLES}/common/inc)
  add_compile_options(-g -fopenmp -O3 -Wextra-DMULTITHREAD)
  add_link_options(-fopenmp -DMULTITHREAD)

else()

  #Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
  add_compile_options(-I${CUDA_SAMPLES}/common/inc)
  add_compile_options(-g -O3 -Wextra)

endif()

# Build components
# add_subdirectory(Diagnostics)
# add_subdirectory(CUDAProb3)
add_subdirectory(manager)
add_subdirectory(throwParms)
add_subdirectory(covariance)
add_subdirectory(splines)
add_subdirectory(samplePDF)
add_subdirectory(mcmc)
