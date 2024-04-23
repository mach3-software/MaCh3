set(HEADERS
    ThrowParms.h
    CovarianceUtils.h
    covarianceBase.h
    covarianceOsc.h
    covarianceXsec.h
)

add_library(Covariance SHARED
    ThrowParms.cpp
    CovarianceUtils.cpp
    covarianceBase.cpp
    covarianceOsc.cpp
    covarianceXsec.cpp
)

set_target_properties(Covariance PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Covariance)


target_link_libraries(Covariance Manager MaCh3CompilerOptions)
target_include_directories(Covariance PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Covariance
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/covariance)

add_library(MaCh3::Covariance ALIAS Covariance)

