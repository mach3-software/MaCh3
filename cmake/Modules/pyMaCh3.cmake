function(setup_pyMaCh3)

  include(CMakePackageConfigHelpers)
  
  ## get pybind dependency
  set(PYBIND11_FINDPYTHON ON)
  CPMFindPackage(
      NAME pybind11
      VERSION 2.13.5
      GITHUB_REPOSITORY "pybind/pybind11"
      GIT_SHALLOW YES
      GIT_TAG v2.13.5
    )
  
  cmake_parse_arguments(
    ARGS
    "" 
    "INSTALL_DIR"
    "BINDING_FILES;LINK_TARGETS;EXTRA_MODULES"
    "${ARGN}"
  )

  if (${ARGS_INSTALL_DIR} STREQUAL "NONE")
    set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
  else ()
    set(INSTALL_DIR ${ARGS_INSTALL_DIR})
  endif()

  cmessage(STATUS "Making pyMaCh3 library")
  cmessage(STATUS "  - Install dir:              ${INSTALL_DIR}")
  cmessage(STATUS "  - Binding definition files: ${ARGS_BINDING_FILES}")
  cmessage(STATUS "  - Linking targets:          ${ARGS_LINK_TARGETS}")
  cmessage(STATUS "  - Extra modules:            ${ARGS_EXTRA_MODULES}")

  ################################# pybind11 stuff ##################################
  ## EM: make a module target out of all the python*Module.cpp files (currently just one...)
  pybind11_add_module(
    _pyMaCh3 MODULE
    ${ARGS_BINDING_FILES}
  )
  ## EM: only works with code compiled with -fPIC enabled.. I think this flag can make things slightly slower
  ## so would be good to find a way around this.
  set_property( TARGET _pyMaCh3 PROPERTY POSITION_INDEPENDENT_CODE ON )
  target_link_libraries( _pyMaCh3 PRIVATE MaCh3::All ${ARGS_LINK_TARGETS} )
  
  ## install our pybind11 object
  install( TARGETS _pyMaCh3 DESTINATION ${INSTALL_DIR}/ )

  ## create directory inside our python module to hold MaCh3 libraries
  file(MAKE_DIRECTORY ${INSTALL_DIR}/lib)       
  file(MAKE_DIRECTORY ${INSTALL_DIR}/lib/core)  ## <- where we will copy all the MaCh3 core libs
  file(MAKE_DIRECTORY ${INSTALL_DIR}/lib/other) ## <- where we will copy any experiment libs
  
  ## copy all the core libraries into our python module lib folder
  if ( MaCh3_LIB_DIR ) ## we are installing in experimet 
    add_custom_command(TARGET _pyMaCh3 POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${MaCh3_LIB_DIR}/* ${INSTALL_DIR}/lib/core 
      COMMAND_EXPAND_LISTS
    )

    set(MaCh3_PYTHON_INIT_TEMPLATE ${MaCh3_SOURCE_DIR}/cmake/Templates/__init__.py.in)
    
  else() ## we are installing core pyMaCh3 - so libs haven't been installed yet so can't just copy them
    install(
      TARGETS
        yaml-cpp
        OscProbCalcer
        Oscillator
        spdlog
        Manager
        SplineDict
        Parameters
        Splines
        Samples
        Fitters
        Plotting
      LIBRARY DESTINATION ${INSTALL_DIR}/lib/core)

      set(MaCh3_PYTHON_INIT_TEMPLATE ${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/__init__.py.in)
  endif()

  ## install any experiment specific libraries into python lib dir
  set(LINK_TARGET_LIB_LIST "")
  foreach(link_target ${ARGS_LINK_TARGETS})
    install( TARGETS ${ARGS_LINK_TARGETS} DESTINATION ${INSTALL_DIR}/lib/other )
    set(LINK_TARGET_LIB_LIST "${LINK_TARGET_LIB_LIST}\"${link_target}\", ")
  endforeach()

  set(PYMACH3_ADDITIONAL_MODULE_IMPORT "")
  set(PYMACH3_ADDITIONAL_MODULE_NAMES  "")

  foreach(module ${ARGS_EXTRA_MODULES})
    set(PYMACH3_ADDITIONAL_MODULE_IMPORT "${PYMACH3_ADDITIONAL_MODULE_IMPORT}, ${module}")
    set(PYMACH3_ADDITIONAL_MODULE_NAMES  "${PYMACH3_ADDITIONAL_MODULE_NAMES}, \"${module}\"")
  endforeach()

  ## set the MaCh3 root directory so can be used in the __init__.py 
  if (NOT MaCh3_PREFIX)
    set(MaCh3_BUILD_PATH ${MaCh3_BINARY_DIR})
  else()
    set(MaCh3_BUILD_PATH ${MaCh3_PREFIX})
  endif()

  ## generate our module __init__ file
  cmessage(STATUS "pyMaCh3 __init__.py template: ${MaCh3_PYTHON_INIT_TEMPLATE}")
  configure_package_config_file(
    ${MaCh3_PYTHON_INIT_TEMPLATE} ${INSTALL_DIR}/__init__.py
    INSTALL_DESTINATION
      /this/is/ignored/for/some/reason/thanks/kitware
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

endfunction()
