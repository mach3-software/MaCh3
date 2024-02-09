set(HEADERS
	covarianceBase.h
	covarianceOsc.h
	covarianceXsec.h
)

add_library(Covariance SHARED
	covarianceBase.cpp
	covarianceOsc.cpp
	covarianceXsec.cpp
)

set_target_properties(Covariance PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Covariance)


target_link_libraries(Covariance ThrowParms ROOT::ROOT Manager)
target_include_directories(Covariance PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Covariance
        EXPORT mach3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/covariance)

add_library(MACH3::Covariance ALIAS Covariance)

