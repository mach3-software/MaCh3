set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

#KS: Enable language, necessary when using CUDA
enable_language(CXX)
SET(MaCh3_VERSION 1.1.0)

# Try to find CUDA
find_package(CUDAToolkit)
add_library(MaCh3CompilerOptions INTERFACE)
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

# Check if CUDA was found
if(MaCh3_GPU_ENABLED)
  include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/CUDASetup.cmake)
endif()

include(ROOT)
if(NOT TARGET ROOT::ROOT)
  cmessage(FATAL_ERROR "MaCh3 Expected dependency target: ROOT::ROOT")
endif()

if(ROOT_VERSION VERSION_LESS 6.18.0)
    message(FATAL_ERROR "Using ROOT version smaller than 6.18.0, this may lead to unexpected results")
endif()

#YAML for reading in config files
set(YAML_CPP_VERSION 0.7.0) #KS: We need it for version.h file also define this number olny once
set(YAML_CPP_GIT_TAG "yaml-cpp-${YAML_CPP_VERSION}")
CPMAddPackage(
    NAME yaml-cpp
    VERSION ${YAML_CPP_VERSION}
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG "${YAML_CPP_GIT_TAG}"
    OPTIONS
      "YAML_BUILD_SHARED_LIBS ON"
)

if(NOT TARGET yaml-cpp::yaml-cpp)
  cmessage(FATAL_ERROR "MaCh3 Expected dependency target: yaml-cpp::yaml-cpp")
endif()

#KS: Not being used right now so commenting it out.
#CPMAddPackage(
#  NAME Eigen
#  VERSION 3.2.8
#  URL https://gitlab.com/libeigen/eigen/-/archive/3.2.8/eigen-3.2.8.tar.gz
#  # Eigen's CMakelists are not intended for library use
#  DOWNLOAD_ONLY YES
#)

#if(Eigen_ADDED)
#  add_library(Eigen INTERFACE IMPORTED)
#  target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
#endif()


################################## Oscillation ################################
#KS: All of these should be moved to separate cmake and be handled by osc class, keep it for now
#If USE_PROB3 not defined turn it off by default
if(NOT DEFINED USE_PROB3)
  SET(USE_PROB3 FALSE)
endif()

# Oscillation calcualtion
# In the future which osc calc we use might be set with a flag
SET(MaCh3_Oscillator_ENABLED "")
if (USE_PROB3)
  CPMFindPackage(
    NAME Prob3plusplus
    VERSION 3.10.3
    GITHUB_REPOSITORY "mach3-software/Prob3plusplus"
    GIT_TAG v3.10.3
  )
  LIST(APPEND MaCh3_Oscillator_ENABLED "Prob3++")
else()
  CPMFindPackage(
    NAME CUDAProb3
    GITHUB_REPOSITORY "mach3-software/CUDAProb3"
    GIT_TAG "feature_cleanup"
    DOWNLOAD_ONLY YES
  )
  LIST(APPEND MaCh3_Oscillator_ENABLED "CUDAProb3")
endif()

#dump_cmake_variables(Prob3plusplus)

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

############################  Setting Flags  ####################################
# Check the where to install 
cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")

# add cmake script files
if(NOT ${CMAKE_SCRIPT_SETUP})  #Check if setup by Experiment MaCh3
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
endif()

# Custom commands
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

#Compile and link options
#KS: Add -Werror  and maybe -Wconversion in the future
target_compile_options(MaCh3CompilerOptions INTERFACE
                                                    -g
                                                    -Wextra
                                                    -Wall
                                                    -pedantic
                                                    -Wshadow
                                                    -Wuninitialized
                                                    -Wnon-virtual-dtor
                                                    -Woverloaded-virtual)
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
  target_compile_definitions(MaCh3CompilerOptions INTERFACE DEBUG=${DEBUG_LEVEL})
  cmessage(STATUS "Enabling DEBUG with Level: \"${DEBUG_LEVEL}\"")
else()
  #KS: Consider in future __attribute__((always_inline)) see https://indico.cern.ch/event/386232/sessions/159923/attachments/771039/1057534/always_inline_performance.pdf
  #https://gcc.gnu.org/onlinedocs/gcc-3.3.6/gcc/Optimize-Options.html
  target_compile_options(MaCh3CompilerOptions INTERFACE -O3 -funroll-loops --param=max-vartrack-size=100000000 -finline-limit=100000000)
  #KS: addomg Link-Time Optimization (LTO)
  target_link_options(MaCh3CompilerOptions INTERFACE -flto)
endif()

#KS: If multithreading is not defined enable it by default
if(NOT DEFINED MaCh3_MULTITHREAD_ENABLED)
  SET(MaCh3_MULTITHREAD_ENABLED TRUE)
endif()

#Add Multithread flags
if(MaCh3_MULTITHREAD_ENABLED)
  target_compile_options(MaCh3CompilerOptions INTERFACE -fopenmp)
  target_link_libraries(MaCh3CompilerOptions INTERFACE gomp)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE MULTITHREAD)
endif()

if(CPU_ONLY)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE CPU_ONLY)
endif()

if(USE_PROB3)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE USE_PROB3)
endif()

if(MaCh3_GPU_ENABLED)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE CUDA)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE GPU_ON)
endif()

if(NOT DEFINED MaCh3_LOW_MEMORY_STRUCTS_ENABLED)
  SET(MaCh3_LOW_MEMORY_STRUCTS_ENABLED FALSE)
endif()

if(MaCh3_LOW_MEMORY_STRUCTS_ENABLED)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE _LOW_MEMORY_STRUCTS_)
endif()

set_target_properties(MaCh3CompilerOptions PROPERTIES EXPORT_NAME CompilerOptions)
install(TARGETS MaCh3CompilerOptions
    EXPORT MaCh3-targets
    LIBRARY DESTINATION lib/)

#KS: Consider moving it somewhere else
SET(MaCh3_Fitter_ENABLED "MR2T2")
SET(MaCh3_MINUIT2_ENABLED FALSE)
if(ROOT_CXX_FLAGS MATCHES "-DMINUIT2_ENABLED")
  LIST(APPEND MaCh3_Fitter_ENABLED " Minuit2")
  SET(MaCh3_MINUIT2_ENABLED TRUE)
endif()

#Include logger
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/Logger.cmake)
################################# Features ##################################

LIST(APPEND ALL_FEATURES
  DEBUG
  MULTITHREAD
  GPU
  LOW_MEMORY_STRUCTS
  Oscillator
  Fitter
  )
cmessage(STATUS "MaCh3 Features: ")
foreach(f ${ALL_FEATURES})
  cmessage(STATUS "     ${f}: ${MaCh3_${f}_ENABLED}")
endforeach()

################################# Build MaCh3 ##################################

# Build components
add_subdirectory(manager)
add_subdirectory(covariance)
add_subdirectory(splines)
if (NOT USE_PROB3)
  add_subdirectory(OscClass)
endif()
add_subdirectory(samplePDF)
add_subdirectory(mcmc)
add_subdirectory(Diagnostics)

#This is to export the target properties of MaCh3
#Anything that links to "MaCh3" will get all of these target properties
add_library(MaCh3 INTERFACE)
target_link_libraries(MaCh3 INTERFACE MCMC OscClass SamplePDF Covariance Splines Manager MaCh3CompilerOptions)
set_target_properties(MaCh3 PROPERTIES EXPORT_NAME All)

install(TARGETS MaCh3
    EXPORT MaCh3-targets
    LIBRARY DESTINATION lib/)

add_library(MaCh3::All ALIAS MaCh3)

configure_file(cmake/Templates/setup.MaCh3.sh.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.MaCh3.sh" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.MaCh3.sh" DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


set(export_destinations
    ${CMAKE_INSTALL_PREFIX}/lib/cmake/
    ${CMAKE_INSTALL_PREFIX}/
)

foreach(dest ${export_destinations})
  install(EXPORT MaCh3-targets
          FILE MaCh3Targets.cmake
          NAMESPACE MaCh3::
          DESTINATION ${dest}
  )
endforeach()


#KS: Options to print dependency graph
if(NOT DEFINED MaCh3_DependancyGraph)
  set(MaCh3_DependancyGraph FALSE)
endif()

if(MaCh3_DependancyGraph)
  add_custom_target(graphviz ALL
                    COMMAND ${CMAKE_COMMAND} "--graphviz=foo.dot" .
                    COMMAND dot -Tpng foo.dot -o foo.png
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
                    COMMENT "Generating dependency graph"
  )
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

include(mach3-config)

############################ Version Info ##########################
## Based on Ewan Miller implementation in OAGenWeightsApp (T2K)
## Make a custom file named version.h that has a
## dump of the versions of the main dependencies of MaCh3
## We also check the git hash of each of the dependencies and info
## on whether or not they were using a tagged commit or not
set(GENERATED_VERSION_HEADER_DIR ${CMAKE_CURRENT_BINARY_DIR}/)

if(EXISTS ${GENERATED_VERSION_HEADER_DIR}/_version.h)
  message(FATAL_ERROR "File \"${GENERATED_VERSION_HEADER_DIR}/_version.h\" found, \
  this should never be created, remove!")
endif()

add_custom_target(
  version_header ALL
  DEPENDS ${GENERATED_VERSION_HEADER_DIR}/_version.h
)

ADD_CUSTOM_COMMAND(
  OUTPUT ${GENERATED_VERSION_HEADER_DIR}/version.h
         ${GENERATED_VERSION_HEADER_DIR}/_version.h

  COMMAND ${CMAKE_COMMAND}
          -DCMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}
          -DCMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
          -DROOT_VERSION=${ROOT_VERSION}
          -DYAML_CPP_VERSION=${YAML_CPP_VERSION}
          -DSPDLOG_VERSION=${SPDLOG_VERSION}
          -DMaCh3_GPU_ENABLED=${MaCh3_GPU_ENABLED}
          -P
          ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/version.cmake
)

set_target_properties(version_header PROPERTIES INCLUDE_DIRECTORIES ${GENERATED_VERSION_HEADER_DIR})
