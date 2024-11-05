

if(CUDA)
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
        $<$<NOT:$<BOOL:${CPU_ONLY}>>:gpuSplineUtils.cu>
    )

    set_target_properties(Splines PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET Splines PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})

elseif(USE_FPGA)
    set(HEADERS
        SplineBase.h
        splineFDBase.h
        SplineMonolith.h
        SplineCommon.h
        SplineStructs.h
        oneAPISplineUtils.dp.hpp
    )
    
    add_library(Splines SHARED
        SplineBase.cpp
        splineFDBase.cpp
        SplineMonolith.cpp
        oneAPISplineUtils.dp.cpp
    )
    set_source_files_properties(oneAPISplineUtils.dp.cpp PROPERTIES COMPILE_FLAGS "-fsycl -fintelfpga")

elseif(CPU_ONLY)
    set(HEADERS
        SplineBase.h
        splineFDBase.h
        SplineMonolith.h
        SplineCommon.h
        SplineStructs.h
    )

    add_library(Splines SHARED
        SplineBase.cpp
        splineFDBase.cpp
        SplineMonolith.cpp
    )
 
endif()

ROOT_GENERATE_DICTIONARY(MaCh3SplineDict
    ${CMAKE_CURRENT_SOURCE_DIR}/SplineCommon.h
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/MaCh3SplinesLinkDef.h
    OPTIONS -p)

set_target_properties(Splines PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Splines)

target_link_libraries(Splines Manager MaCh3CompilerOptions)

target_include_directories(Splines PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS Splines
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/splines)

add_library(MaCh3::Splines ALIAS Splines)

