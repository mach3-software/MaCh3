set(Splines_header_files 
	splineBase.h 
	splineInterface.h 
	splineFDBase.h 
        splineMonolith.h)

set(Splines_implementation_files 
	splineBase.cpp 
	splineFDBase.cpp
	SplineMonolith.cpp
        gpuSplineUtils.cu)
# Note: SplineMonolith is for GPU


set(Splines_implementation_files_CPU 
	splineBase.cpp 
	splineFDBase.cpp)

if(NOT ${CPU_ONLY})

  add_library(Splines SHARED ${Splines_implementation_files})
  target_include_directories(Splines PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
  add_library(gpuSpline gpuSplineUtils.cu)
  set_target_properties(Splines PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  if(${CUDA_ERROR_CHECK})
	target_compile_options(Splines PUBLIC -DCUDA_ERROR_CHECK)
  endif()
  set_property(TARGET gpuSpline Splines PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)
  target_link_libraries(Splines ThrowParms gpuSpline)

else()

  add_library(Splines SHARED ${Splines_implementation_files_CPU})
  target_include_directories(Splines PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
  target_link_libraries(Splines ThrowParms)

endif()

set_target_properties(Splines PROPERTIES 
  PUBLIC_HEADER "${Splines_header_files}"
  EXPORT_NAME Splines)

install(TARGETS Splines
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/splines)

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

