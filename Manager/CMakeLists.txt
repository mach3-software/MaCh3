set(HEADERS
    Manager.h
    MaCh3Modes.h
    YamlHelper.h
    TablePrinter.h
    MaCh3Logger.h
    Monitor.h
    MaCh3Exception.h
    gpuUtils.cuh
    Random.h
    Core.h
)

add_library(Manager SHARED
    Manager.cpp
    MaCh3Modes.cpp
    Monitor.cpp
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:gpuUtils.cu>
)

set_target_properties(Manager PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Manager)

if(MaCh3_GPU_ENABLED)
    set_target_properties(Manager PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #KS: In cmake 3.18 this is no longer needed https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html#cmake-cuda-architectures
    set_property(TARGET Manager PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES}) 
endif()

#If compiling with GPU is not enabled MaCh3GPUCompilerOptions will be empty
target_link_libraries(Manager PUBLIC MaCh3CompilerOptions MaCh3GPUCompilerOptions yaml-cpp spdlog NuOscillator ROOT::Tree ROOT::Hist ROOT::Physics)
target_link_libraries(Manager PRIVATE MaCh3Warnings)

# KS: To use std::filesystem for older gcc version this is needed
# See https://en.cppreference.com/cpp/filesystem#Notes
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(Manager PUBLIC stdc++fs)
endif()

target_include_directories(Manager PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Manager
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Manager)

add_library(MaCh3::Manager ALIAS Manager)
target_link_libraries(MaCh3 INTERFACE Manager)


install(TARGETS yaml-cpp
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/yaml-cpp
        )

install(DIRECTORY ${yaml-cpp_SOURCE_DIR}/include/
    DESTINATION include/
)

install(TARGETS spdlog
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/spdlog
        )

install(DIRECTORY ${spdlog_SOURCE_DIR}/include/
    DESTINATION include/
)
