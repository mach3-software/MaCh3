function(setup_pyMaCh3)

  cmake_parse_arguments(
    ARGS
    "" 
    "TARGET_NAME;INSTALL_DIR"
    "BINDING_FILES;LINK_TARGETS;EXTRA_MODULES"
    "${ARGN}"
  )

  if (${ARGS_INSTALL_DIR} STREQUAL "NONE")
    set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
  else ()
    set(INSTALL_DIR ${ARGS_INSTALL_DIR})
  endif()

  message("Making pyMaCh3 library")
  message("Target name:              ${ARGS_TARGET_NAME}")
  message("install dir:              ${INSTALL_DIR}")
  message("Binding definition files: ${ARGS_BINDING_FILES}")
  message("linking targets:          ${ARGS_LINK_TARGETS}")
  message("extra modules:            ${ARGS_EXTRA_MODULES}")

  ################################# pybind11 stuff ##################################
  ## EM: make a module target out of all the python*Module.cpp files (currently just one...)
  pybind11_add_module(
    ${ARGS_TARGET_NAME} MODULE
    ${ARGS_BINDING_FILES}
  )
  ## EM: only works with code compiled with -fPIC enabled.. I think this flag can make things slightly slower
  ## so would be good to find a way around this.
  set_property( TARGET ${ARGS_TARGET_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON )
  target_link_libraries( ${ARGS_TARGET_NAME} PRIVATE MaCh3::All ${ARGS_LINK_TARGETS} )

  message ( "INSTALLING pyMaCh3 to ${INSTALL_DIR}" )
  
  ## install our pybind11 object
  install( TARGETS ${ARGS_TARGET_NAME} DESTINATION ${INSTALL_DIR}/ )

  ## create directory inside our python module to hold MaCh3 libraries
  file(MAKE_DIRECTORY ${INSTALL_DIR}/lib)       
  file(MAKE_DIRECTORY ${INSTALL_DIR}/lib/core)  ## <- where we will copy all the MaCh3 core libs
  file(MAKE_DIRECTORY ${INSTALL_DIR}/lib/other) ## <- where we will copy any experiment libs
  
  ## copy all the core libraries into our python module lib folder
  if ( NOT MaCh3_LIB_DIR )
    set(MaCh3_LIB_DIR ${PROJECT_BINARY_DIR}/lib)
  endif()

  add_custom_command(TARGET ${ARGS_TARGET_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${MaCh3_LIB_DIR}/* ${INSTALL_DIR}/lib/core
    COMMAND_EXPAND_LISTS
  )
  
  ## install any experiment specific libraries into python lib dir
  set(LINK_TARGET_LIB_LIST "")
  foreach(link_target ${ARGS_LINK_TARGETS})
    message("link_target: ${link_target}")
    install( TARGETS ${ARGS_LINK_TARGETS} DESTINATION ${INSTALL_DIR}/lib/other )
    set(LINK_TARGET_LIB_LIST "${LINK_TARGET_LIB_LIST}\"${link_target}\", ")
  endforeach()

  set(PYMACH3_ADDITIONAL_MODULE_IMPORT "")
  set(PYMACH3_ADDITIONAL_MODULE_NAMES  "")

  foreach(module ${ARGS_EXTRA_MODULES})
    set(PYMACH3_ADDITIONAL_MODULE_IMPORT "${PYMACH3_ADDITIONAL_MODULE_IMPORT}, ${module}")
    set(PYMACH3_ADDITIONAL_MODULE_NAMES  "${PYMACH3_ADDITIONAL_MODULE_NAMES}, \"${module}\"")
  endforeach()


  ## generate our module __init__ file
  message("__init__.py template file: ${MaCh3_PYTHON_INIT_TEMPLATE}")
  configure_package_config_file(
    ${MaCh3_PYTHON_INIT_TEMPLATE} ${INSTALL_DIR}/__init__.py
    INSTALL_DESTINATION
      /this/is/ignored/for/some/reason/thanks/kitware
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

endfunction()
