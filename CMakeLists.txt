set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

#KS: Enable language, necessary when using CUDA
enable_language(CXX)
set(MaCh3_VERSION 1.1.7)

# Try to find CUDA
find_package(CUDAToolkit QUIET)
add_library(MaCh3CompilerOptions INTERFACE)
# Check if CUDA was found

message(STATUS "USE_CPU flag: ${USE_CPU}")
message(STATUS "USE_FPGA flag: ${USE_FPGA}")

if(CUDAToolkit_FOUND AND NOT(USE_CPU OR USE_FPGA))
  message(STATUS "CUDA found. Adding CUDA support.")
  project(MaCh3 VERSION ${MaCh3_VERSION} LANGUAGES CXX CUDA)
  set(MaCh3_GPU_ENABLED TRUE)
  set(CUDA TRUE)
  set(CPU_ONLY FALSE)
else()
  project(MaCh3 VERSION ${MaCh3_VERSION} LANGUAGES CXX)
  if(USE_CPU)  
    message(STATUS "CPU only mode, disabling CUDA")
    set(CPU_ONLY TRUE)
    set(USE_FPGA FALSE)
  elseif(USE_FPGA)
    message(STATUS "FPGA mode, disabling CUDA") 
    set(MaCh3_GPU_ENABLED FALSE)
    set(CUDA FALSE)
    set(CPU_ONLY FALSE)
  endif()
endif()

if(USE_FPGA)
  if(NOT DEFINED FPGA_DEVICE)
    set(FPGA_DEVICE "Agilex7")
    message(STATUS "FPGA_DEVICE was not specified.\
                    \nConfiguring the design to the default FPGA family target: ${FPGA_DEVICE}\
                    \nPlease refer to the README for information on target selection.")
  endif()

  set(USER_FPGA_FLAGS ${USER_FPGA_FLAGS})

  set(EMULATOR_TARGET fpga_emu)
  set(SIMULATOR_TARGET fpga_sim)
  set(REPORT_TARGET report)
  set(FPGA_TARGET fpga)

  # Set the names of the generated files per makefile target
  set(EMULATOR_OUTPUT_NAME ${TARGET_NAME}.${EMULATOR_TARGET})
  set(SIMULATOR_OUTPUT_NAME ${TARGET_NAME}.${SIMULATOR_TARGET})
  set(REPORT_OUTPUT_NAME ${TARGET_NAME}.${REPORT_TARGET})
  set(FPGA_OUTPUT_NAME ${TARGET_NAME}.${FPGA_TARGET})

  message(STATUS "Additional USER_FPGA_FLAGS=${USER_FPGA_FLAGS}")
  message(STATUS "Additional USER_FLAGS=${USER_FLAGS}")
  set(QACTYPES "-qactypes")

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
#KS: Store some handy cmake functions
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/MaCh3Utils.cmake)
#Loads all dependencies
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/MaCh3Dependencies.cmake)

############################  C++ Compiler  ####################################
if (NOT DEFINED CMAKE_CXX_STANDARD OR "${CMAKE_CXX_STANDARD} " STREQUAL " ")
  set(CMAKE_CXX_STANDARD 14)
endif()

# KS: If C++ standard is lower than C++ standard used for ROOT compilation things will go terribly wrong
if(DEFINED ROOT_CXX_STANDARD AND ROOT_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
endif()
cmessage(STATUS "CMAKE CXX Standard: ${CMAKE_CXX_STANDARD}")

# KS: ROOT changed cmake in 6.32, we should move avay from using Luke's hack, keep it for due to compatibility
if(ROOT_CXX_STANDARD LESS 14 AND ROOT_VERSION VERSION_LESS 6.32.00)
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

# KS: Compile and link options for more see https://github.com/cpp-best-practices/cppbestpractices/tree/master
target_compile_options(MaCh3CompilerOptions INTERFACE
    -g                      # Generate debug information
    -Wextra                 # Enable extra warning flags
    -Wall                   # Enable all standard warning flags
    -pedantic               # Enforce strict ISO compliance (all versions of GCC, Clang >= 3.2)
    -Wshadow                # Warn when a variable declaration shadows one from an outer scope
    -Wuninitialized         # Warn about uninitialized variables
    -Wnon-virtual-dtor      # Warn when a class with virtual functions has a non-virtual destructor
    -Woverloaded-virtual    # Warn when a function declaration hides a virtual function from a base class
    -Wformat=2              # Warn on security issues around functions that format output (ie printf)
    -Wunused                # Warn on anything being unused
    -Wredundant-decls       # Warn about multiple declarations of the same entity. Useful for code cleanup.
    #-Wstrict-aliasing=2    # Helps detect potential aliasing issues that could lead to undefined behavior.
    #-Wuseless-cast         # Warn if you perform a cast to the same type (only in GCC >= 4.8)
    #-Wpadded               # Warn when padding is added to a structure or class for alignment
    #-Wnull-dereference     # Warn if a null dereference is detected (only in GCC >= 6.0)
    #-Wold-style-cast       # Warn for c-style casts
    #-Wconversion           # Warn on type conversions that may lose data
    #-Werror                # Treat Warnings as Errors
)
# KS Some compiler options are only available in GCC, in case we move to other compilers we will have to expand this
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(MaCh3CompilerOptions INTERFACE
    -Wlogical-op            # Warn about logical operations being used where bitwise were probably wanted (only in GCC)
    -Wduplicated-cond       # Warn if if / else chain has duplicated conditions (only in GCC >= 6.0)
    -Wduplicated-branches   # Warn if if / else branches have duplicated code (only in GCC >= 7.0)
  )
endif()

if(USE_FPGA)
  target_compile_options(MaCh3CompilerOptions INTERFACE
    -DUSE_FPGA              # preprocessor directive to switch on FPGA code
  )
endif()

#KS: If Debug is not defined disable it by default
DefineEnabledRequiredSwitch(MaCh3_DEBUG_ENABLED FALSE)

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
  target_compile_definitions(MaCh3CompilerOptions INTERFACE DEBUG=${DEBUG_LEVEL})
  cmessage(STATUS "Enabling DEBUG with Level: \"${DEBUG_LEVEL}\"")
else()
  #KS: Consider in future __attribute__((always_inline)) see https://indico.cern.ch/event/386232/sessions/159923/attachments/771039/1057534/always_inline_performance.pdf
  #https://gcc.gnu.org/onlinedocs/gcc-3.3.6/gcc/Optimize-Options.html
  target_compile_options(MaCh3CompilerOptions INTERFACE
    -O3                                   # Optimize code for maximum speed
    -finline-limit=100000000              # Increase the limit for inlining functions to improve performance
    # KS: After benchmarking below didn't in fact worse performance, leave it for future tests and documentation
    #-funroll-loops                       # Unroll loops where possible for performance
    #--param=max-vartrack-size=100000000  # Set maximum size of variable tracking data to avoid excessive memory usage
    #-flto                                # Enable link-time optimization (commented out for now, needs more testing)
  )
  #KS: add Link-Time Optimization (LTO)
  #target_link_libraries(MaCh3CompilerOptions INTERFACE -flto)
endif()

#KS: If multithreading is not defined enable it by default
DefineEnabledRequiredSwitch(MaCh3_MULTITHREAD_ENABLED TRUE)

#Add Multithread flags
if(MaCh3_MULTITHREAD_ENABLED)
  if(USE_FPGA)
    target_compile_options(MaCh3CompilerOptions INTERFACE -qopenmp)
  else()
    target_compile_options(MaCh3CompilerOptions INTERFACE -fopenmp)
  endif()
  target_link_libraries(MaCh3CompilerOptions INTERFACE gomp)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE MULTITHREAD)
endif()

if(CPU_ONLY)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE CPU_ONLY)
endif()

if(USE_FPGA)
  target_link_libraries(MaCh3CompilerOptions INTERFACE sycl)
  target_compile_options(MaCh3CompilerOptions INTERFACE -Wall ${QACTYPES} ${USER_FLAGS})
endif()

if(USE_PROB3)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE USE_PROB3)
endif()

if(MaCh3_GPU_ENABLED)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE CUDA)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE GPU_ON)
endif()

if(NOT DEFINED MaCh3_LOW_MEMORY_STRUCTS_ENABLED)
  set(MaCh3_LOW_MEMORY_STRUCTS_ENABLED FALSE)
endif()

if(MaCh3_LOW_MEMORY_STRUCTS_ENABLED)
  target_compile_definitions(MaCh3CompilerOptions INTERFACE _LOW_MEMORY_STRUCTS_)
endif()

set_target_properties(MaCh3CompilerOptions PROPERTIES EXPORT_NAME CompilerOptions)
install(TARGETS MaCh3CompilerOptions
    EXPORT MaCh3-targets
    LIBRARY DESTINATION lib/)

#Include logger
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/Logger.cmake)
################################# Features ##################################
#KS: Retrieve and print the compile options and link libraries for maximal verbose
get_target_property(compile_options MaCh3CompilerOptions INTERFACE_COMPILE_OPTIONS)
get_target_property(link_libraries MaCh3CompilerOptions INTERFACE_LINK_LIBRARIES)

cmessage(STATUS "Compile options for MaCh3CompilerOptions: ${compile_options}")
cmessage(STATUS "Link libraries for MaCh3CompilerOptions: ${link_libraries}")

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
if (NOT USE_PROB3)
  add_subdirectory(OscClass)
endif()
add_subdirectory(samplePDF)
add_subdirectory(mcmc)
add_subdirectory(Diagnostics)
add_subdirectory(plotting)


################################# pybind11 stuff ##################################

if( MaCh3_PYTHON_ENABLED )
  ## EM: make a module target out of all the python*Module.cpp files (currently just one...)
  pybind11_add_module( 
    pyMaCh3 MODULE 
    pyMaCh3.cpp 
    plotting/plottingUtils/pythonPlottingModule.cpp
  )
  ## EM: only works with code compiled with -fPIC enabled.. I think this flag can things slightly slower
  ## so would be good to find a way around this. 
  set_property( TARGET pyMaCh3 PROPERTY POSITION_INDEPENDENT_CODE ON )
  target_link_libraries( pyMaCh3 PUBLIC Plotting )
  install( TARGETS pyMaCh3 DESTINATION pyMaCh3/)
endif()


#This is to export the target properties of MaCh3
#Anything that links to "MaCh3" will get all of these target properties
add_library(MaCh3 INTERFACE)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MaCh3 INTERFACE OpenMP::OpenMP_CXX)
endif()
target_link_libraries(MaCh3 INTERFACE MCMC OscClass SamplePDF Covariance Splines Manager MaCh3CompilerOptions Plotting)
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

#KS: Options to print dependency graph
DefineEnabledRequiredSwitch(MaCh3_DependancyGraph FALSE)

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
