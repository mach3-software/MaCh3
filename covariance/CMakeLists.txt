set(Covariance_implementation_files
	covarianceBase.cpp 
	covarianceOsc.cpp 
	covarianceXsec.cpp)


set(Covariance_header_files
	covarianceBase.cpp 
	covarianceOsc.cpp 
	covarianceXsec.cpp)

add_library(Covariance SHARED ${Covariance_implementation_files})
target_include_directories(Covariance PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(Covariance)
add_dependencies(Covariance throwParms)


install(TARGETS Covariance DESTINATION lib)

