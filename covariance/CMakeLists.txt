set(HEADERS
    CovarianceUtils.h
    covarianceBase.h
    covarianceOsc.h
    covarianceXsec.h
    AdaptiveMCMCHandler.h
    PCAHandler.h
)

add_library(Covariance SHARED
    covarianceBase.cpp
    covarianceOsc.cpp
    covarianceXsec.cpp
    AdaptiveMCMCHandler.cpp
    PCAHandler.cpp
)

set_target_properties(Covariance PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Covariance)


if(MaCh3_GPU_ENABLED)
    set_target_properties(Covariance PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET Covariance PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

target_link_libraries(Covariance PUBLIC Manager)
target_link_libraries(Covariance PRIVATE MaCh3Warnings)

target_include_directories(Covariance PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Covariance
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/covariance)

add_library(MaCh3::Covariance ALIAS Covariance)

