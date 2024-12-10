ROOT_GENERATE_DICTIONARY(MaCh3SplineDict
    ${CMAKE_CURRENT_SOURCE_DIR}/SplineCommon.h
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/MaCh3SplinesLinkDef.h
    OPTIONS -p)

set(HEADERS
    SplineBase.h
    splineFDBase.h
    SplineMonolith.h
    SplineCommon.h
    SplineStructs.h
    gpuSplineUtils.cuh
)

add_library(Splines SHARED
    SplineBase.cpp
    splineFDBase.cpp
    SplineMonolith.cpp
    #MaCh3SplineDict.cxx #KS: This will have to be sperate library because root like C-style casting too much
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:gpuSplineUtils.cu>
)

if(NOT CPU_ONLY)
    set_target_properties(Splines PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET Splines PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

set_target_properties(Splines PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Splines)

target_link_libraries(Splines PUBLIC Manager MaCh3CompilerOptions MaCh3GPUCompilerOptions)
target_link_libraries(Splines PRIVATE MaCh3Warnings)

target_include_directories(Splines PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libMaCh3SplineDict_rdict.pcm DESTINATION lib/)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libMaCh3SplineDict.rootmap DESTINATION lib/)

install(TARGETS Splines
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/splines)

add_library(MaCh3::Splines ALIAS Splines)

