
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
install( DIRECTORY pyMaCh3 DESTINATION ./ )
install( TARGETS _pyMaCh3 DESTINATION pyMaCh3/ )
