#KS: CUDA Prob is terrible so you cannot compile both beam and atmos at the same time as you get nvlink error: Multiple definition. This leads to issues further down the line but can't be bothered right now. Solution would be to use extern rather then symbols.

if (USE_PROB3)
    cmessage(FATAL_ERROR "OscClass only defined using CUDAProb3. USE_PROB3 is defined and indicates Prob3++/ProbGPU is being expected to be used.")
endif()


######################Atmos######################
if (NOT CPU_ONLY)
  set(HEADERS
      ${CUDAProb3_SOURCE_DIR}/atmoscudapropagator.cuh
  )

  add_library(Atmos SHARED
      ${CUDAProb3_SOURCE_DIR}/atmoscudapropagator.cu
  )

  set_target_properties(Atmos PROPERTIES
      PUBLIC_HEADER "${HEADERS}"
      EXPORT_NAME Atmos)

  target_link_libraries(Atmos Manager MaCh3CompilerOptions)

  target_include_directories(Atmos PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
      $<BUILD_INTERFACE:${CUDAProb3_SOURCE_DIR}/>
      $<INSTALL_INTERFACE:include>
  )

  if(NOT CPU_ONLY)
      set_target_properties(Atmos PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
      #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
      set_property(TARGET Atmos PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
  endif()

  install(TARGETS Atmos
          EXPORT MaCh3-targets
          LIBRARY DESTINATION lib/
          PUBLIC_HEADER DESTINATION include/OscClass
          )

  add_library(MaCh3::Atmos ALIAS Atmos)
endif()
######################OscClass######################

set(HEADERS
    OscClass_CUDAProb3.h
    ${CUDAProb3_SOURCE_DIR}/beamcudapropagator.cuh
    ${CUDAProb3_SOURCE_DIR}/atmoscudapropagator.cuh
    ${CUDAProb3_SOURCE_DIR}/constants.hpp
    ${CUDAProb3_SOURCE_DIR}/hpc_helpers.cuh
    ${CUDAProb3_SOURCE_DIR}/beamcpupropagator.hpp
    ${CUDAProb3_SOURCE_DIR}/atmoscpupropagator.hpp
    ${CUDAProb3_SOURCE_DIR}/propagator.hpp
    ${CUDAProb3_SOURCE_DIR}/types.hpp
    ${CUDAProb3_SOURCE_DIR}/math.hpp
    ${CUDAProb3_SOURCE_DIR}/physics.hpp
    ${CUDAProb3_SOURCE_DIR}/cpupropagator.hpp
    ${CUDAProb3_SOURCE_DIR}/cuda_unique.cuh
)

add_library(OscClass SHARED
    OscClass_CUDAProb3.cpp
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:${CUDAProb3_SOURCE_DIR}/beamcudapropagator.cu>>
)

set_target_properties(OscClass PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME OscClass)


target_link_libraries(OscClass Manager MaCh3CompilerOptions)
if (NOT CPU_ONLY)
  target_link_libraries(OscClass Atmos)
endif()

target_include_directories(OscClass PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<BUILD_INTERFACE:${CUDAProb3_SOURCE_DIR}/>
    $<INSTALL_INTERFACE:include>
)

if(NOT CPU_ONLY)
    set_target_properties(OscClass PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET OscClass PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

install(TARGETS OscClass
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/OscClass
        )

add_library(MaCh3::OscClass ALIAS OscClass)

