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
  $<$<BOOL:${MaCh3_GPU_ENABLED}>:gpuSplineUtils.cuh>     
  $<$<BOOL:${USE_FPGA}>:${CMAKE_CURRENT_SOURCE_DIR}/oneAPIUtils.dp.hpp>
)

add_library(Splines SHARED
  SplineBase.cpp
  splineFDBase.cpp
  SplineMonolith.cpp
  $<$<BOOL:${MaCh3_GPU_ENABLED}>:gpuSplineUtils.cu>
  $<$<BOOL:${USE_FPGA}>:${CMAKE_CURRENT_SOURCE_DIR}/oneAPIUtils.dp.cpp>
)


if (USE_FPGA)

    if(COMPILE_EMULATOR)
      message(STATUS "Compiling for FPGA_EMULATOR")
      target_compile_options(Splines PUBLIC
	-DUSE_FPGA
	-DFPGA_EMULATOR
	-fsycl
	-fintelfpga
      )
      target_link_options(Splines PRIVATE
	-fsycl
	-fintelfpga
	-Xstarget=${FPGA_DEVICE}
      )
    elseif(COMPILE_HARDWARE)
      message(STATUS "Compiling for FPGA_HARDWARE")
      target_compile_options(Splines PUBLIC
	-DUSE_FPGA
	-DFPGA_HARDWARE
	-fsycl
	-fintelfpga
	-Xshardware
	-Xstarget=${FPGA_DEVICE}

      )
      target_link_options(Splines PRIVATE
	-fsycl
	-fintelfpga
	-Xshardware
	-Xstarget=${FPGA_DEVICE}
      )
    elseif(COMPILE_REPORT)
      message(STATUS "Compiling for FPGA REPORT")
      target_compile_options(Splines PUBLIC
	-DUSE_FPGA
	-DFPGA_HARDWARE
	-fsycl
	-fintelfpga
	-Xshardware
	-Xstarget=${FPGA_DEVICE}
	-fsycl-link=early
      )
      target_link_options(Splines PRIVATE
	-fsycl
	-fintelfpga
	-Xshardware
	-Xstarget=${FPGA_DEVICE}
	# -fsycl-link=early
      )
    endif()
endif()

if(MaCh3_GPU_ENABLED)
  set_target_properties(Splines PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
      #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
  set_property(TARGET Splines PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

    

set_target_properties(Splines PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Splines)

target_link_libraries(Splines PUBLIC Covariance MaCh3CompilerOptions MaCh3GPUCompilerOptions)
target_link_libraries(Splines PRIVATE MaCh3Warnings)

target_include_directories(Splines PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
    /opt/intel/oneapi/2025.0/include
)

install(TARGETS Splines
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/splines)

add_library(MaCh3::Splines ALIAS Splines)


# KS: Now we generate root dictionary, we can't add it to above as ROOT is bit dumb and things will fail due to MaCh3Warnings etc
add_library(SplineDict SHARED MaCh3SplineDict.cxx)
target_include_directories(SplineDict PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(SplineDict PUBLIC Splines MaCh3CompilerOptions MaCh3GPUCompilerOptions)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libMaCh3SplineDict_rdict.pcm DESTINATION lib/)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libMaCh3SplineDict.rootmap DESTINATION lib/)
