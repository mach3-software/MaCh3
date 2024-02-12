set(HEADERS
    splineBase.h
    splineInterface.h
    splineFDBase.h
    SplineMonolith.h
)

add_library(Splines SHARED
    splineBase.cpp
    splineFDBase.cpp
    SplineMonolith.cpp
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:gpuSplineUtils.cu>
)


if(NOT ${CPU_ONLY})
    set_target_properties(Splines PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET Splines PROPERTY CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES})
endif()

if(MaCh3_DEBUG_ENABLED)
    target_compile_definitions(Splines PUBLIC CUDA_ERROR_CHECK)
endif()

set_target_properties(Splines PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Splines)

target_link_libraries(Splines Covariance MaCh3CompilerOptions)

target_include_directories(Splines PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:$ENV{CUDAPATH}/include>>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:$ENV{CUDAPATH}/common/inc>>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:$ENV{CUDAPATH}/samples/common/inc>>
)

install(TARGETS Splines
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/splines)

add_library(MaCh3::Splines ALIAS Splines)

