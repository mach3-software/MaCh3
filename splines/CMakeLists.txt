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
  target_include_directories(Splines PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
  target_compile_features(Splines PUBLIC)
  add_library(gpuSpline gpuSplineUtils.cu)
  set_target_properties(Splines PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_property(TARGET gpuSpline Splines PROPERTY CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)
  add_dependencies(Splines throwParms gpuSpline)

else()

  add_library(Splines SHARED ${Splines_implementation_files_CPU})
  target_include_directories(Splines PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
  target_compile_features(Splines PUBLIC)
  add_dependencies(Splines throwParms)

endif()

install(TARGETS Splines DESTINATION lib)

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

