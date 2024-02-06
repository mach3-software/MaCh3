set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)  

# Try to find CUDA
find_package(CUDAToolkit)
# Check if CUDA was found
if(CUDAToolkit_FOUND AND NOT(USE_CPU))
    message(STATUS "CUDA found. Adding CUDA support.")
    SET(MaCh3_GPU_ENABLED TRUE)
    SET(CPU_ONLY FALSE)
else()
    message(STATUS "CUDA not found. Proceeding without CUDA support.")
    SET(MaCh3_GPU_ENABLED FALSE)
    SET(CPU_ONLY TRUE)
endif()

SET(MaCh3_VERSION 0.0.2)
if(CPU_ONLY)
  project(MaCh3 VERSION ${MaCh3_VERSION} LANGUAGES CXX)
else()
  project(MaCh3 VERSION ${MaCh3_VERSION} LANGUAGES CXX CUDA)
endif()

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

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
################################## Dependencies ################################

include(CPM)

#Luke's handing cmake modules which Neutrino hep experiments might want
CPMFindPackage(
      NAME CMakeModules
      GIT_TAG stable
      GITHUB_REPOSITORY NuHepMC/CMakeModules
      DOWNLOAD_ONLY
  )
include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)
include(NuHepMCUtils)
include(ROOT)
if(NOT TARGET ROOT::ROOT)
  cmessage(FATAL_ERROR "MaCh3 Expected dependency target: ROOT::ROOT")
endif()

#YAML for reading in config files
CPMFindPackage(
    NAME yaml-cpp
    VERSION 0.7.0
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG "yaml-cpp-0.7.0"
    OPTIONS
      "YAML_BUILD_SHARED_LIBS ON"
)

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


# If we're doing GPU stuff, we need the CUDA helper module
if(MaCh3_GPU_ENABLED)
    CPMAddPackage(
        NAME cuda-samples
        GITHUB_REPOSITORY "NVIDIA/cuda-samples"
        GIT_TAG v12.3
	DOWNLOAD_ONLY YES
    )
# Now we add the library
    if(cuda-samples_ADDED)
         add_library(cuda-samples INTERFACE IMPORTED)
         target_include_directories(cuda-samples INTERFACE ${cuda-samples_SOURCE_DIR})
    endif()
endif()


################################## Oscillation ################################

#If USE_PROB3 not defined turn it off by default
if(NOT DEFINED USE_PROB3)
  SET(USE_PROB3 FALSE)
endif()

# Oscillation calcualtion
# In the future which osc calc we use might be set with a flag

if (USE_PROB3)
  CPMFindPackage(
      NAME Prob3plusplus
      VERSION 3.10.3
      GITHUB_REPOSITORY "mach3-software/Prob3plusplus"
      GIT_TAG v3.10.3
  )
else()
  if (CPU_ONLY)
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

dump_cmake_variables(Prob3plusplus)

############################  Setting Flags  ####################################
# Check the where to install 
cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")

# add cmake script files
if(NOT ${CMAKE_SCRIPT_SETUP})  #Check if setup by Experiment MaCh3
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
  #CUDA
  if(NOT CPU_ONLY)
    include(${CMAKE_SOURCE_DIR}/cmake/CUDASetup.cmake)
  endif()
endif()

# Custom commands
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

#Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
add_compile_options(-g -Wextra)

#KS: If Debug is not defined disable it by default
if(NOT DEFINED MaCh3_DEBUG_ENABLED)
  SET(MaCh3_DEBUG_ENABLED FALSE)
endif()

#If DEBUG_LEVEL was defined but MaCh3_DEBUG_ENABLED not, enable debug flag
if(DEFINED DEBUG_LEVEL)
  SET(MaCh3_DEBUG_ENABLED TRUE)
else()
  #If MaCh debug was enable but level not, set it to 1. In very rare cases we want to go beyond 1.
  if(MaCh3_DEBUG_ENABLED)
    SET(DEBUG_LEVEL 1)
  endif()
endif()

#KS: If Debug add debugging compile flag if not add optimisation for speed
if(MaCh3_DEBUG_ENABLED)
  add_compile_definitions(DEBUG=${DEBUG_LEVEL})
  cmessage(STATUS "Enabling DEBUG with Level: \"${DEBUG_LEVEL}\"")
else()
  add_compile_options(-O3)
endif()

#KS: If multithreading is not defined enable it by default
if(NOT DEFINED MaCh3_MULTITHREAD_ENABLED)
  SET(MaCh3_MULTITHREAD_ENABLED TRUE)
endif()

#Add MultiThread flags
if(MaCh3_MULTITHREAD_ENABLED)
  add_compile_options(-fopenmp)
  add_compile_definitions(MULTITHREAD)
endif()

if(CPU_ONLY)
  add_compile_definitions(CPU_ONLY)
else()
  add_compile_options(-I${CUDA_SAMPLES}/common/inc)
  add_compile_options(-I$ENV{CUDAPATH}/targets/x86_64-linux/include)
endif()

if(USE_PROB3)
  add_compile_definitions(USE_PROB3)
endif()

############################  C++ Compiler  ####################################
if (NOT DEFINED CMAKE_CXX_STANDARD OR "${CMAKE_CXX_STANDARD} " STREQUAL " ")
  SET(CMAKE_CXX_STANDARD 11)
endif()

if(DEFINED ROOT_CXX_STANDARD AND ROOT_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
endif()
cmessage(STATUS "CMAKE CXX Standard: ${CMAKE_CXX_STANDARD}")

if(${ROOT_CXX_STANDARD} LESS 14)
  cmessage(WARNING "ROOT CXX STANDARD: ${ROOT_CXX_STANDARD}")
endif()

################################# Features ##################################

LIST(APPEND ALL_FEATURES
  DEBUG
  MULTITHREAD
  GPU
  )
cmessage(STATUS "MaCh3 Features: ")
foreach(f ${ALL_FEATURES})
  cmessage(STATUS "     ${f}: ${MaCh3_${f}_ENABLED}")
endforeach()

################################# Build MaCh3 ##################################

# Build components
# add_subdirectory(Diagnostics)
add_subdirectory(manager)
add_subdirectory(throwParms)
add_subdirectory(covariance)
add_subdirectory(splines)
if (NOT USE_PROB3)
  add_subdirectory(OscClass)
endif()
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

############################  Install  ####################################

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/MaCh3Config.cmake.in ${CMAKE_BINARY_DIR}/MaCh3Config.cmake
  INSTALL_DESTINATION
    /this/is/ignored/for/some/reason/thanks/kitware
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

  write_basic_package_version_file(${CMAKE_BINARY_DIR}/MaCh3ConfigVersion.cmake
  VERSION ${MaCh3_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES
    ${CMAKE_BINARY_DIR}/MaCh3Config.cmake
    ${CMAKE_BINARY_DIR}/MaCh3ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/MaCh3)
