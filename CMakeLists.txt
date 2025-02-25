set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(MaCh3 VERSION 1.4.5 LANGUAGES CXX)
set(MaCh3_VERSION ${PROJECT_VERSION})

option(MaCh3_PYTHON_ENABLED "Whether to build MaCh3 python bindings" OFF)
option(MaCh3_WERROR_ENABLED "Whether to build MaCh3 with heightened compiler pedancy" ON)
option(MaCh3_LOW_MEMORY_STRUCTS_ENABLED "This will use float/short string for many structures" OFF)
option(MaCh3_MULTITHREAD_ENABLED "This will enable multithreading (OMP) in whole framework" ON)
option(MaCh3_DEBUG_ENABLED "Enable special debugging mode, with additional printouts and data structure, increases RAM usage and low performance" OFF)
option(MaCh3_DependancyGraph "Will produce dependency graph" OFF)
option(MaCh3_GPU_BENCHMARK "Perform fancy CUDA Benchmarking" OFF)
option(MaCh3_NATIVE_ENABLED "Enable native CPU optimizations for improved performance during benchmarking)" OFF)
option(MaCh3_GPU_ENABLED "Use GPU acceleration or not" ON)

#KS: Load cmake function like DefineEnabledRequiredSwitch allowing to write more compact cmake
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules)
include(MaCh3Utils)

add_library(MaCh3CompileDefinitions INTERFACE)
#Create a separate target for GPU compiler options
#if GPU is not enabled nothing will be added to this and it will be empty
# this saves a few if statements
add_library(MaCh3GPUCompilerOptions INTERFACE)

find_package(CUDAToolkit QUIET)
# Check if CUDA was found
if(CUDAToolkit_FOUND AND NOT(MaCh3_GPU_ENABLED))
  cmessage(STATUS "CUDA found. Adding CUDA support.")
  set(MaCh3_GPU_ENABLED ON)
  enable_language(CUDA)
  set(CPU_ONLY FALSE)
else()
  cmessage(STATUS "CUDA not found. Proceeding without CUDA support.")
  set(MaCh3_GPU_ENABLED OFF)
  set(CPU_ONLY TRUE)
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

################################## Dependencies ################################
include(MaCh3Dependencies)

############################  C++ Compiler  ####################################
if (NOT DEFINED CMAKE_CXX_STANDARD OR "${CMAKE_CXX_STANDARD} " STREQUAL " ")
  set(CMAKE_CXX_STANDARD 14)
  cmessage(STATUS "Set default CXX standard: \"${CMAKE_CXX_STANDARD}\"")
endif()

if(DEFINED ROOT_CXX_STANDARD AND NOT ROOT_CXX_STANDARD EQUAL CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
  cmessage(STATUS "Set CXX standard due to ROOT: \"${ROOT_CXX_STANDARD}\"")
endif()
cmessage(STATUS "CMAKE_CXX_STANDARD: \"${CMAKE_CXX_STANDARD}\"")
set(MACH3_CXX_STANDARD ${CMAKE_CXX_STANDARD})
set_property(GLOBAL PROPERTY MACH3_CXX_STANDARD "${MACH3_CXX_STANDARD}")

set(MACH3_CXX_STANDARD ${CMAKE_CXX_STANDARD})
set_property(GLOBAL PROPERTY MACH3_CXX_STANDARD "${MACH3_CXX_STANDARD}")
############################  Setting Flags  ####################################
cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")

# add cmake script files
if(NOT ${CMAKE_SCRIPT_SETUP})  #Check if setup by Experiment MaCh3
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
endif()

# Custom commands
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# KS: Compile and link options for more see https://github.com/cpp-best-practices/cppbestpractices/tree/master
add_library(MaCh3Warnings INTERFACE)
target_compile_options(MaCh3Warnings INTERFACE
  -Wextra                 # Enable extra warning flags
  -Wall                   # Enable all standard warning flags
  -Wshadow                # Warn when a variable declaration shadows one from an outer scope
  -Wuninitialized         # Warn about uninitialized variables
  -Wnon-virtual-dtor      # Warn when a class with virtual functions has a non-virtual destructor
  -Woverloaded-virtual    # Warn when a function declaration hides a virtual function from a base class
  -Wformat=2              # Warn on security issues around functions that format output (ie printf)
  -Wunused                # Warn on anything being unused
  -Wredundant-decls       # Warn about multiple declarations of the same entity. Useful for code cleanup.
  -Wstrict-aliasing=2     # Helps detect potential aliasing issues that could lead to undefined behavior.
  -Wnull-dereference      # Warn if a null dereference is detected (only in GCC >= 6.0)
  -Wold-style-cast        # Warn for c-style casts
  -Wconversion            # Warn on type conversions that may lose data
  -Wformat-security       # Warn on functions that are potentially insecure for formatting
  -Walloca                # Warn if `alloca` is used, as it can lead to stack overflows
  -Wswitch-enum           # Warn if a `switch` statement on an enum does not cover all values
  #-Wfloat-equal          # Warn if floating-point values are compared directly
  #-Wpadded               # Warn when padding is added to a structure or class for alignment
)
# KS Some compiler options are only available in GCC, in case we move to other compilers we will have to expand this
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(MaCh3Warnings INTERFACE
    -Wlogical-op            # Warn about logical operations being used where bitwise were probably wanted (only in GCC)
    -Wduplicated-cond       # Warn if if / else chain has duplicated conditions (only in GCC >= 6.0)
    -Wduplicated-branches   # Warn if if / else branches have duplicated code (only in GCC >= 7.0)
    -Wuseless-cast          # Warn if you perform a cast to the same type (only in GCC >= 4.8)
  )
endif()

# KS: Clang is super picky, it is almost impossible to make it work with Werror
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(MaCh3_WERROR_ENABLED FALSE)
endif()

if(MaCh3_WERROR_ENABLED)
target_compile_options(MaCh3Warnings INTERFACE
  -Werror                # Treat Warnings as Errors
)
endif()

add_library(MaCh3CompilerOptions INTERFACE)
target_link_libraries(MaCh3CompilerOptions INTERFACE MaCh3CompileDefinitions)
target_compile_options(MaCh3CompilerOptions INTERFACE
    -g                      # Generate debug information
    -pedantic               # Enforce strict ISO compliance (all versions of GCC, Clang >= 3.2)
)

#If DEBUG_LEVEL was defined but MaCh3_DEBUG_ENABLED not, enable debug flag
if(DEFINED DEBUG_LEVEL)
  set(MaCh3_DEBUG_ENABLED TRUE)
else()
  #If MaCh3 debug was enable but level not, set it to 1. In very rare cases we want to go beyond 1.
  if(MaCh3_DEBUG_ENABLED)
    set(DEBUG_LEVEL 1)
  endif()
endif()

#KS: If Debug add debugging compile flag if not add optimisation for speed
if(MaCh3_DEBUG_ENABLED)
  target_compile_options(MaCh3CompilerOptions INTERFACE
    -O0                                  # Turn off any optimisation to have best debug experience
  )
  target_compile_definitions(MaCh3CompileDefinitions INTERFACE DEBUG=${DEBUG_LEVEL})
  cmessage(STATUS "Enabling DEBUG with Level: \"${DEBUG_LEVEL}\"")
else()
  #KS: Consider in future __attribute__((always_inline)) see https://indico.cern.ch/event/386232/sessions/159923/attachments/771039/1057534/always_inline_performance.pdf
  #https://gcc.gnu.org/onlinedocs/gcc-3.3.6/gcc/Optimize-Options.html
  target_compile_options(MaCh3CompilerOptions INTERFACE
    -O3                                   # Optimize code for maximum speed
    #-fstrict-aliasing                    # Assume that pointers of different types do not point to the same memory location
    # KS: After benchmarking below didn't in fact worse performance, leave it for future tests and documentation
    #-funroll-loops                       # Unroll loops where possible for performance
    #--param=max-vartrack-size=100000000  # Set maximum size of variable tracking data to avoid excessive memory usage
    #-flto                                # Enable link-time optimization (commented out for now, needs more testing)
  )
  # KS Some compiler options are only available in GCC, in case we move to other compilers we will have to expand this
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(MaCh3CompilerOptions INTERFACE
      -finline-limit=100000000              # Increase the limit for inlining functions to improve performance
    )
  endif()
  #KS: add Link-Time Optimization (LTO)
  #target_link_libraries(MaCh3CompilerOptions INTERFACE -flto)
endif()

if(MaCh3_NATIVE_ENABLED)
  target_compile_options(MaCh3CompilerOptions INTERFACE -march=native)
endif()

#Add Multithread flags
if(MaCh3_MULTITHREAD_ENABLED)
  target_compile_options(MaCh3CompilerOptions INTERFACE -fopenmp)
  # KS: Not a clue why clang need different...
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_libraries(MaCh3CompilerOptions INTERFACE omp)
  else()
    target_link_libraries(MaCh3CompilerOptions INTERFACE gomp)
  endif()
  target_compile_definitions(MaCh3CompileDefinitions INTERFACE MULTITHREAD)
endif()

if(MaCh3_GPU_ENABLED) 
  target_compile_definitions(MaCh3GPUCompilerOptions INTERFACE CUDA)
  target_compile_definitions(MaCh3GPUCompilerOptions INTERFACE GPU_ON)
endif()

if(MaCh3_LOW_MEMORY_STRUCTS_ENABLED)
  target_compile_definitions(MaCh3CompileDefinitions INTERFACE _LOW_MEMORY_STRUCTS_)
endif()

set_target_properties(MaCh3CompilerOptions PROPERTIES EXPORT_NAME CompilerOptions)
set_target_properties(MaCh3GPUCompilerOptions PROPERTIES EXPORT_NAME GPUCompilerOptions)
install(TARGETS MaCh3CompilerOptions MaCh3GPUCompilerOptions MaCh3CompileDefinitions
    EXPORT MaCh3-targets
    LIBRARY DESTINATION lib/)

#Include logger
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/Logger.cmake)
#Setup NuOscillator, this has to be here as it requires all flags to be setup
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/NuOscillatorSetup.cmake)

################################# Features ##################################
#KS: Retrieve and print the compile options and link libraries for maximal verbose
get_target_property(compile_options MaCh3CompilerOptions INTERFACE_COMPILE_OPTIONS)
get_target_property(link_libraries MaCh3CompilerOptions INTERFACE_LINK_LIBRARIES)

cmessage(STATUS "Compile options for MaCh3CompilerOptions: ${compile_options}")
cmessage(STATUS "Link libraries for MaCh3CompilerOptions: ${link_libraries}")

if(MaCh3_GPU_ENABLED)
  #KS: Retrieve and print the compile options and link libraries for maximal verbose for GPU
  get_target_property(gpu_compile_options MaCh3GPUCompilerOptions INTERFACE_COMPILE_OPTIONS)
  get_target_property(gpu_link_libraries MaCh3GPUCompilerOptions INTERFACE_LINK_LIBRARIES)
  #Nice print out of the GPU options if GPU is enabled
  cmessage(STATUS "Compile options for MaCh3GPUCompilerOptions: ${gpu_compile_options}")
  cmessage(STATUS "Link libraries for MaCh3GPUCompilerOptions: ${gpu_link_libraries}")
endif()

LIST(APPEND ALL_FEATURES
  DEBUG
  MULTITHREAD
  GPU
  PYTHON
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
add_subdirectory(samplePDF)
add_subdirectory(mcmc)
add_subdirectory(Diagnostics)
add_subdirectory(plotting)
if (MaCh3_PYTHON_ENABLED)
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
  add_subdirectory(python)
endif()

#This is to export the target properties of MaCh3
#Anything that links to "MaCh3" will get all of these target properties
add_library(MaCh3 INTERFACE)
target_link_libraries(MaCh3 INTERFACE MCMC SamplePDF Covariance Splines Manager MaCh3CompilerOptions MaCh3GPUCompilerOptions Plotting)
set_target_properties(MaCh3 PROPERTIES EXPORT_NAME All)

install(TARGETS MaCh3
    EXPORT MaCh3-targets
    LIBRARY DESTINATION lib/)

add_library(MaCh3::All ALIAS MaCh3)

configure_file(cmake/Templates/setup.MaCh3.sh.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.MaCh3.sh" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.MaCh3.sh" DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

install(EXPORT MaCh3-targets
        FILE MaCh3Targets.cmake
        NAMESPACE MaCh3::
        DESTINATION ${CMAKE_INSTALL_PREFIX}/
)
install(DIRECTORY cmake DESTINATION ${CMAKE_BINARY_DIR})

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
if(EXISTS "${CMAKE_INSTALL_PREFIX}/_version.h")
  cmessage(FATAL_ERROR "File \"${CMAKE_INSTALL_PREFIX}/_version.h\" found, this should never be created, remove!")
endif()

# Define custom target for version_header
add_custom_target(
  version_header ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/_version.h
)

# Add custom command to generate version.h and _version.h
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.h
         ${CMAKE_CURRENT_BINARY_DIR}/_version.h
  COMMAND ${CMAKE_COMMAND}
          -DCMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}
          -DCMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
          -DROOT_VERSION=${ROOT_VERSION}
          -DMaCh3_VERSION=${MaCh3_VERSION}
          -DYAML_CPP_VERSION=${YAML_CPP_VERSION}
          -DSPDLOG_VERSION=${SPDLOG_VERSION}
          -DMaCh3_GPU_ENABLED=${MaCh3_GPU_ENABLED}
          -P ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/version.cmake
  COMMENT "Generating version.h and _version.h"
)

# Set properties for the target
set_target_properties(version_header PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR})

# Install version.h to the install directory
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/version.h
    DESTINATION ${CMAKE_INSTALL_PREFIX}/)

# uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

# KS: Configure the Doxygen input file, this is to ensure whenever we update MaCh3 version Doxyfile will have same version.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/Doxyfile.in ${CMAKE_CURRENT_SOURCE_DIR}/Doc/Doxyfile @ONLY)
