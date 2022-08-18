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
target_include_directories(MCMC PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(MCMC yaml-cpp)

install(TARGETS MCMC DESTINATION lib)



