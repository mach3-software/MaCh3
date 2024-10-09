
################################# pybind11 stuff ##################################

if( MaCh3_PYTHON_ENABLED )
  ## EM: make a module target out of all the python*Module.cpp files (currently just one...)
  pybind11_add_module( 
    pyMaCh3 MODULE 
    pyMaCh3.cpp 
    plotting.cpp
    fitter.cpp
    samplePDF.cpp
    manager.cpp
  )
  ## EM: only works with code compiled with -fPIC enabled.. I think this flag can make things slightly slower
  ## so would be good to find a way around this. 
  set_property( TARGET pyMaCh3 PROPERTY POSITION_INDEPENDENT_CODE ON )
  target_link_libraries( pyMaCh3 PUBLIC MaCh3::All )
  install( TARGETS pyMaCh3 DESTINATION pyMaCh3/)
endif()
