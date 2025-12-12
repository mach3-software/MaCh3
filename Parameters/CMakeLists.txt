set(HEADERS
    ParameterHandlerUtils.h
    ParameterHandlerBase.h
    ParameterHandlerGeneric.h
    AdaptiveMCMCHandler.h
    PCAHandler.h
    ParameterTunes.h
    ParameterStructs.h
)

add_library(Parameters SHARED
    ParameterHandlerBase.cpp
    ParameterHandlerGeneric.cpp
    AdaptiveMCMCHandler.cpp
    PCAHandler.cpp
    ParameterTunes.cpp
)

set_target_properties(Parameters PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Parameters)


if(MaCh3_GPU_ENABLED)
    set_target_properties(Parameters PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET Parameters PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

target_link_libraries(Parameters PUBLIC Manager Eigen3::Eigen)
# KS: ROOT ≤ 6.30 requires explicit linking of MathMore (if built, not sure why through...)
if(ROOT_VERSION VERSION_LESS 6.30)
    if(TARGET ROOT::MathMore)
        target_link_libraries(Parameters PUBLIC ROOT::MathMore)
    else()
        message(WARNING "ROOT::MathMore target not found; continuing without it")
    endif()
endif()

target_link_libraries(Parameters PRIVATE MaCh3Warnings)

target_include_directories(Parameters PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Parameters
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Parameters
)

add_library(MaCh3::Parameters ALIAS Parameters)

