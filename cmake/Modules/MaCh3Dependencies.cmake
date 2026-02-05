# Load all MaCh3 dependencies
# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.42.1/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

# Check if CUDA was found
if(MaCh3_GPU_ENABLED)
  include(${CMAKE_CURRENT_LIST_DIR}/CUDASetup.cmake)
endif()

### Begin ROOT setup
### KS: 6.18 is minimal version which has nasty bug-fix to SetBinError in TH2Poly
### https://github.com/root-project/root/commit/9ed733d77ca4b513ef68ba5343c9b2664f24ebc3
### KS: 6.20 is minimal due to ROOT_GENERATE_DICTIONARY
find_package(ROOT 6.20 REQUIRED)

STRING(STRIP "${ROOT_CXX_FLAGS}" ROOT_CXX_FLAGS_LIST)
STRING(REPLACE " " ";" ROOT_CXX_FLAGS_LIST ${ROOT_CXX_FLAGS_LIST})

list (FIND ROOT_CXX_FLAGS_LIST "-std=c++14" CPP14_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++1y" CPP1Y_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++17" CPP17_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++1z" CPP1Z_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++20" CPP20_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++23" CPP23_INDEX)
list (FIND ROOT_CXX_FLAGS_LIST "-std=c++2b" CPP2B_INDEX)

if (CPP14_INDEX GREATER -1 OR CPP1Y_INDEX GREATER -1)
  set(ROOT_CXX_STANDARD 14)
elseif (CPP17_INDEX GREATER -1 OR CPP1Z_INDEX GREATER -1)
  set(ROOT_CXX_STANDARD 17)
elseif (CPP20_INDEX GREATER -1)
  set(ROOT_CXX_STANDARD 20)
elseif (CPP23_INDEX GREATER -1 OR CPP2B_INDEX GREATER -1)
  set(ROOT_CXX_STANDARD 23)
endif()

cmessage(STATUS "ROOT_CXX_FLAGS: \"${ROOT_CXX_FLAGS}\" -> ROOT_CXX_STANDARD: ${ROOT_CXX_STANDARD}")

execute_process(COMMAND root-config --features
  OUTPUT_VARIABLE ROOT_CONFIG_FEATURES OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE " " ";" ROOT_FEATURES_LIST "${ROOT_CONFIG_FEATURES}")
# Check if "minuit2" is in the list of ROOT features
list(FIND ROOT_FEATURES_LIST "minuit2" ROOT_CONFIG_MINUIT2)

# KS: Since ROOT 6.32.0 Minuit is turned on by default
set(MaCh3_MINUIT2_ENABLED FALSE)
if(ROOT_VERSION GREATER_EQUAL 6.32.00 OR ROOT_CONFIG_MINUIT2 GREATER -1)
  set(MaCh3_MINUIT2_ENABLED TRUE)
endif()

# KS: Parts of MaCh3 uses Fast Fourier Transform, do some safeguard
list(FIND ROOT_FEATURES_LIST "fftw3" ROOT_CONFIG_FFT)
set(MaCh3_FFT_ENABLED FALSE)
if(ROOT_CONFIG_FFT GREATER -1)
  cmessage(STATUS "FFT Enabled")
  set(MaCh3_FFT_ENABLED TRUE)
  target_compile_definitions(MaCh3CompileDefinitions INTERFACE MaCh3_FFT)
endif()
### End ROOT setup

#YAML for reading in config files
set(YAML_CPP_VERSION 0.9.0) #KS: We need it for version.h file also define this number only once
CPMAddPackage(
    NAME yaml-cpp
    VERSION ${YAML_CPP_VERSION}
    GITHUB_REPOSITORY "jbeder/yaml-cpp"
    GIT_TAG yaml-cpp-${YAML_CPP_VERSION}
    GIT_SHALLOW NO
    OPTIONS
      "YAML_CPP_INSTALL ON"
      "YAML_CPP_BUILD_TESTS OFF"
      "YAML_CPP_BUILD_CONTRIB OFF"
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
LIST(APPEND MaCh3_Fitter_ENABLED " DelayedMR2T2")
LIST(APPEND MaCh3_Fitter_ENABLED " PSO")
if(MaCh3_MINUIT2_ENABLED)
  LIST(APPEND MaCh3_Fitter_ENABLED " Minuit2")
  target_compile_definitions(MaCh3CompileDefinitions INTERFACE MaCh3_MINUIT2)
endif()

######################### python binding ##########################
if( MaCh3_PYTHON_ENABLED )
  set(PYBIND11_FINDPYTHON ON)
  CPMFindPackage(
      NAME pybind11
      VERSION 2.13.5
      GITHUB_REPOSITORY "pybind/pybind11"
      GIT_SHALLOW YES
      GIT_TAG v2.13.5
    )
endif()
