function(setup_pyMaCh3)

  cmake_parse_arguments(
    ARGS
    "" 
    "TARGET_NAME;INSTALL_DIR"
    "BINDING_FILES;TARGETS"
    "${ARGN}"
  )

  message("ARGS_TARGET_NAME:   ${ARGS_TARGET_NAME}")
  message("ARGS_INSTALL_DIR:   ${ARGS_INSTALL_DIR}")
  message("ARGS_BINDING_FILES: ${ARGS_BINDING_FILES}")
  message("ARGS_TARGETS:       ${ARGS_TARGETS}")

  ################################# pybind11 stuff ##################################
  ## EM: make a module target out of all the python*Module.cpp files (currently just one...)
  pybind11_add_module(
    ${ARGS_TARGET_NAME} MODULE
    ${ARGS_BINDING_FILES}
  )
  ## EM: only works with code compiled with -fPIC enabled.. I think this flag can make things slightly slower
  ## so would be good to find a way around this.
  set_property( TARGET ${ARGS_TARGET_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON )
  target_link_libraries( ${ARGS_TARGET_NAME} PRIVATE MaCh3::All NuOscillator MaCh3Warnings )

  if (${ARGS_INSTALL_DIR} STREQUAL "NONE")
    set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
  else ()
    set(INSTALL_DIR ${ARGS_INSTALL_DIR})
  endif()

  message ( "INSTALLING pyMaCh3 to ${INSTALL_DIR}" )

  install( TARGETS ${ARGS_TARGET_NAME} DESTINATION ${INSTALL_DIR}/ )

  message("__init__.py template file: ${MaCh3_PYTHON_INIT_TEMPLATE}")
  configure_package_config_file(
    ${MaCh3_PYTHON_INIT_TEMPLATE} ${INSTALL_DIR}/__init__.py
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

endfunction()
