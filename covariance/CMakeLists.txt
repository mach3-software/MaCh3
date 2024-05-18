set(HEADERS
    ThrowParms.h
    CovarianceUtils.h
    covarianceBase.h
    covarianceOsc.h
    covarianceXsec.h
)

add_library(Covariance SHARED
    ThrowParms.cpp
    covarianceBase.cpp
    covarianceOsc.cpp
    covarianceXsec.cpp
)

set_target_properties(Covariance PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Covariance)


if(NOT CPU_ONLY)
    set_target_properties(Covariance PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET Covariance PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

target_link_libraries(Covariance Manager MaCh3CompilerOptions)
target_include_directories(Covariance PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Covariance
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/covariance)

add_library(MaCh3::Covariance ALIAS Covariance)

