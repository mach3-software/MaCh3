set(Covariance_implementation_files
	covarianceBase.cpp 
	covarianceOsc.cpp 
	covarianceXsec.cpp)


set(Covariance_header_files
	covarianceBase.cpp 
	covarianceOsc.cpp 
	covarianceXsec.cpp)

add_library(Covariance SHARED ${Covariance_implementation_files})
target_include_directories(Covariance INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(Covariance ThrowParms ROOT::ROOT Manager)

set_target_properties(Covariance PROPERTIES 
  PUBLIC_HEADER "${Covariance_header_files}"
  EXPORT_NAME Covariance)

install(TARGETS Covariance
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/covariance)
