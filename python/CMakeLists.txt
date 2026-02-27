set(HEADERS
  plotting.h
  fitters.h
  samples.h
  manager.h
  parameters.h
  splines.h
  pyMaCh3.h
)

add_library(pyMaCh3 INTERFACE)

set_target_properties(pyMaCh3 PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME pyMaCh3)

target_include_directories(
  pyMaCh3 
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>
)

install(TARGETS pyMaCh3
    EXPORT MaCh3-targets
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/python
)

add_library(MaCh3::pyMaCh3 ALIAS pyMaCh3)
