set(MACH3_LIB_LIST "")
foreach(lib IN LISTS mach3_libs)
  string(APPEND MACH3_LIB_LIST " -l${lib}")
endforeach()
string(STRIP "${MACH3_LIB_LIST}" MACH3_LIB_LIST)

set(MACH3_FEATURES_LIST)
set(ALL_FEATURES_CONFIG ${ALL_FEATURES})
list(REMOVE_ITEM ALL_FEATURES_CONFIG Oscillator Fitter)
foreach(feature IN LISTS ALL_FEATURES_CONFIG)
  if(MaCh3_${feature}_ENABLED)
    list(APPEND MACH3_FEATURES_LIST "${feature}")
  endif()
endforeach()
list(APPEND MACH3_FEATURES_LIST ${MaCh3_Fitter_ENABLED})
list(APPEND MACH3_FEATURES_LIST ${MaCh3_Oscillator_ENABLED})
string(REPLACE ";" " " MACH3_FEATURES "${MACH3_FEATURES_LIST}")

set(MACH3_ROOT_INCLUDE_FLAGS "${ROOT_INCLUDE_DIRS}")
string(REPLACE ";" " -I" MACH3_ROOT_INCLUDE_FLAGS "-I${MACH3_ROOT_INCLUDE_FLAGS}")
set(MACH3_ROOT_DEFINITIONS "${ROOT_DEFINITIONS}")
string(REPLACE ";" " " MACH3_ROOT_DEFINITIONS "${MACH3_ROOT_DEFINITIONS}")
set(MACH3_CFLAGS "-I${CMAKE_INSTALL_PREFIX}/include ${MACH3_ROOT_INCLUDE_FLAGS} ${MACH3_ROOT_DEFINITIONS}")
set(MACH3_LIBS "-L${CMAKE_INSTALL_PREFIX}/lib ${MACH3_LIB_LIST}")
set(MACH3_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")
set(MACH3_INCDIR "${CMAKE_INSTALL_PREFIX}/include")
string(TIMESTAMP MACH3_CONFIG_DATE "%d-%m-%Y")

function(mach3_escape_cpp_string input_variable output_variable)
  set(value "${${input_variable}}")
  string(REPLACE "\\" "\\\\" value "${value}")
  string(REPLACE "\"" "\\\"" value "${value}")
  string(REPLACE "\n" "\\n" value "${value}")
  string(REPLACE "\r" "\\r" value "${value}")
  set(${output_variable} "${value}" PARENT_SCOPE)
endfunction()

foreach(variable IN ITEMS
    MACH3_INCDIR MACH3_CFLAGS MACH3_LIBDIR
    MACH3_LIBS MACH3_FEATURES MaCh3_VERSION MACH3_CONFIG_DATE)
  mach3_escape_cpp_string(${variable} MACH3_CPP_${variable})
endforeach()

set(MACH3_CONFIG_HEADER_DIR "${PROJECT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${MACH3_CONFIG_HEADER_DIR}/CLI")
configure_file(
  "${CMAKE_CURRENT_LIST_DIR}/../Templates/MaCh3CLIConfig.hpp.in"
  "${MACH3_CONFIG_HEADER_DIR}/CLI/ConfigData.hpp"
  @ONLY
)
target_include_directories(mach3 PRIVATE "${MACH3_CONFIG_HEADER_DIR}")
target_include_directories(CLIModules PRIVATE "${MACH3_CONFIG_HEADER_DIR}")
