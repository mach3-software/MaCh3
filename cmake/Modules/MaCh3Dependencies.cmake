# Load all MaCh3 dependencies
# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.38.8/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

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
  include(${CMAKE_CURRENT_LIST_DIR}/CUDASetup.cmake)
endif()

include(ROOT)
if(NOT TARGET ROOT::ROOT)
  cmessage(FATAL_ERROR "MaCh3 Expected dependency target: ROOT::ROOT")
endif()

if(ROOT_VERSION VERSION_LESS 6.18.00)
  cmessage(FATAL_ERROR "Using ROOT version smaller than 6.18.0, this may lead to unexpected results")
endif()

# KS: Since ROOT 6.32.0 Minuit is turned on by default
set(MaCh3_MINUIT2_ENABLED FALSE)
if(ROOT_VERSION GREATER_EQUAL 6.32.00 OR ROOT_CXX_FLAGS MATCHES "-DMINUIT2_ENABLED")
  set(MaCh3_MINUIT2_ENABLED TRUE)
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


set(MaCh3_Fitter_ENABLED "MR2T2")
LIST(APPEND MaCh3_Fitter_ENABLED " PSO")
if(MaCh3_MINUIT2_ENABLED)
  LIST(APPEND MaCh3_Fitter_ENABLED " Minuit2")
  target_compile_definitions(MaCh3CompilerOptions INTERFACE MaCh3_MINUIT2)
endif()


######################### python binding ##########################
# EM: If Debug is not defined disable it by default
DefineEnabledRequiredSwitch(MaCh3_PYTHON_ENABLED FALSE)

if( MaCh3_PYTHON_ENABLED )
  set(PYBIND11_FINDPYTHON ON)

  CPMFindPackage(
      NAME pybind11
      VERSION 2.13.5
      GITHUB_REPOSITORY "pybind/pybind11"
      GIT_TAG v2.13.5
    )
endif()
