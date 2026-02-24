
################################# pybind11 stuff ##################################
## EM: make a module target out of all the python*Module.cpp files (currently just one...)
pybind11_add_module(
  _pyMaCh3 MODULE
  pyMaCh3.cpp
  plotting.cpp
  fitters.cpp
  samples.cpp
  manager.cpp
  parameters.cpp
  splines.cpp
)
## EM: only works with code compiled with -fPIC enabled.. I think this flag can make things slightly slower
## so would be good to find a way around this.
set_property( TARGET _pyMaCh3 PROPERTY POSITION_INDEPENDENT_CODE ON )
target_link_libraries( _pyMaCh3 PRIVATE MaCh3::All NuOscillator yaml-cpp::yaml-cpp MaCh3Warnings )

## check if the current install dir is called <something>/pyMaCh3
## if it is then we are probably building with pip and so don't need to make 
## new pyMaCh3 directory below
## (a smidge hacky but should work)
cmake_path(GET CMAKE_INSTALL_PREFIX FILENAME INSTALL_DIR_NAME)

if (NOT ${INSTALL_DIR_NAME} STREQUAL "pyMaCh3")

  install( DIRECTORY pyMaCh3 DESTINATION ${CMAKE_INSTALL_PREFIX} )
  set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/pyMaCh3 )

else ()

  set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX})

endif()

message ( "INSTALLING pyMaCh3 to ${INSTALL_DIR}")

install( TARGETS _pyMaCh3 DESTINATION ${INSTALL_DIR}/ )

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/Templates/__init__.py.in ${INSTALL_DIR}/__init__.py
  INSTALL_DESTINATION
    /this/is/ignored/for/some/reason/thanks/kitware
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

## make a symlink to the lib dir so that __init__.py can find all the MaCh3.so files
install(CODE "execute_process( \
    COMMAND ${CMAKE_COMMAND} -E create_symlink \
    ${PROJECT_BINARY_DIR}/lib \
    ${INSTALL_DIR}/lib   \
    )"
)
