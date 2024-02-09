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
    set_property(TARGET Splines PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)
endif()


set_target_properties(Splines PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Splines)

target_link_libraries(Splines Covariance)
target_include_directories(Splines PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Splines
        EXPORT mach3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/splines)

add_library(MACH3::Splines ALIAS Splines)


#############################
# For testing purposes...
#cuda_compile_ptx(
#  cuda_ptx_files
#  ${Splines_cuda_files}
#  OPTIONS -Xptxas -dlcm=ca,-allow-expensive-optimizations=true,-fmad=true,-O3,-warn-lmem-usage,-warn-spills -Xcompiler -fpic,-Wall,-Wextra,-Werror,-c -lcudart -L${CUDA_TOOLKIT_ROOT_DIR}/lib64 -I${CUDA_TOOLKIT_ROOT_DIR}/include -I${CUDA_TOOLKIT_ROOT_DIR}/common/inc -lcudart -I${CUDA_TOOLKIT_ROOT_DIR}/#samples/common/inc
#)
# later options: -DCMAKE_CUDA_FLAGS

#add_custom_target(ptx ALL
#    DEPENDS ${cuda_ptx_files} 
#    SOURCES ${Splines_cuda_files} )

