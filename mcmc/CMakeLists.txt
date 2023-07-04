set(mcmc_implementation_files
	mcmc.cpp
	tune.cpp 
	stretch.cpp 
	MCMCProcessor.cpp)

set(mcmc_header_files
	mcmc.h
	tune.h 
	stretch.h 
	MCMCProcessor.h)

add_library(MCMC SHARED ${mcmc_implementation_files})
target_include_directories(MCMC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(MCMC Manager ROOT::ROOT)

set_target_properties(MCMC PROPERTIES 
  PUBLIC_HEADER "${mcmc_header_files}"
  EXPORT_NAME MCMC)

install(TARGETS MCMC
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/mcmc)
