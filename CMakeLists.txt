set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)  

SET(MaCh3_VERSION 0.0.2)
if(${CPU_ONLY})
  project(MaCh3 VERSION ${MaCh3_VERSION} LANGUAGES CXX)
else()
  project(MaCh3 VERSION ${MaCh3_VERSION} LANGUAGES CXX CUDA)
endif()
set(CMAKE_CXX_STANDARD 17)

# Would be good to have a check on CMAKE_CXX_STANDARD here
# Changes default install path to be a subdirectory of the build dir.
# Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif()

# Use the compilers found in the path
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

#include(CPM)
include(${CPM_DOWNLOAD_LOCATION})

#Luke's handing cmake modules which Neutrino hep experiments might want
CPMFindPackage(
    NAME CMakeModules
    GIT_TAG v0.0.11
    GITHUB_REPOSITORY NuHepMC/CMakeModules
    DOWNLOAD_ONLY
)
include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)

#YAML for reading in config files
CPMFindPackage(
    NAME yaml-cpp
    VERSION 0.7.0
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG "yaml-cpp-0.7.0"
    OPTIONS
      "YAML_BUILD_SHARED_LIBS ON"
)

# Oscillation calcualtion
# In the future which osc calc we use might be set with a flag

if (${USE_PROB3})

  CPMFindPackage(
      NAME Prob3plusplus
      VERSION 3.10.3
      GITHUB_REPOSITORY "mach3-software/Prob3plusplus"
      GIT_TAG v3.10.3
  )

else()

  if (${CPU_ONLY})

    CPMFindPackage(
        NAME CUDAProb3
        GITHUB_REPOSITORY "mach3-software/CUDAProb3"
        GIT_TAG "feature_cleanup"
    )

  else()
  
    CPMFindPackage(
        NAME CUDAProb3
        GITHUB_REPOSITORY "mach3-software/CUDAProb3"
        GIT_TAG "feature_cleanup"
	    OPTIONS "GPU_ON ON"
    )

  endif()

endif()

CPMAddPackage(
  NAME Eigen
  VERSION 3.2.8
  URL https://gitlab.com/libeigen/eigen/-/archive/3.2.8/eigen-3.2.8.tar.gz
  # Eigen's CMakelists are not intended for library use
  DOWNLOAD_ONLY YES 
)

if(Eigen_ADDED)
  add_library(Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
endif()

include(NuHepMCUtils)
dump_cmake_variables(Prob3plusplus)

#Nice coloured messaged from cmake
include(CMessage)

include(ROOT)
#Check to see if you successfully found ROOT or not
if(NOT ROOT_FOUND)
  cmessage(FATAL_ERROR " couldn't find ROOT dependency...")
endif()

# Would be good to have a check on CMAKE_CXX_STANDARD here
if(${ROOT_CXX_STANDARD} LESS 14)
  cmessage(WARNING "ROOT CXX STANDARD: ${ROOT_CXX_STANDARD}")
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
endif()

# Custom commands
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

#Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
add_compile_options(-I${CUDA_SAMPLES}/common/inc)
add_compile_options(-g -O3 -Wextra)

#Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
add_compile_options(-g -Wextra)

#If Single thread is chosen MultiThread flags are omitted
if(NOT ${SINGLE_THREAD_ONLY})
  add_compile_options(-DMULTITHREAD -fopenmp)
  add_link_options(-DMULTITHREAD -fopenmp)
endif()

if(${CPU_ONLY})
  add_compile_options(-DCPU_ONLY)
else()
  add_compile_options(-I${CUDA_SAMPLES}/common/inc)
  add_compile_options(-I$ENV{CUDAPATH}/targets/x86_64-linux/include)
endif()

if(${USE_PROB3})
  add_compile_options(-DUSE_PROB3)
endif()

# Build components
# add_subdirectory(Diagnostics)
add_subdirectory(manager)
add_subdirectory(throwParms)
add_subdirectory(covariance)
add_subdirectory(splines)
add_subdirectory(OscClass)
add_subdirectory(samplePDF)
add_subdirectory(mcmc)
#ETA - commenting out as not c++11 compatible...
#add_subdirectory(yaml_test)

#This is to export the target properties of MaCh3
#Anything that links to "MaCh3" will get all of these target properties
if(NOT TARGET MaCh3)
  add_library(MaCh3 INTERFACE)
  set_target_properties(MaCh3 PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}"
      INTERFACE_LINK_LIBRARIES "MCMC;OscClass;SamplePDF;Covariance;Splines;ThrowParms;Manager"
      EXPORT_NAME MaCh3
  )
  cmessage(STATUS "Added target MaCh3")
endif()
