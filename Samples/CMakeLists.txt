set(HEADERS
    BinningHandler.h
    SampleHandlerInterface.h
    SampleHandlerBase.h
    SampleHandlerBase.txx
    SampleHandlerFunctional.h
    SampleStructs.h
    SampleInfo.h
    EventInfo.h
    OscillationHandler.h
)

add_library(Samples SHARED
    BinningHandler.cpp
    SampleHandlerInterface.cpp
    SampleHandlerBase.cpp
    SampleHandlerReweight.cpp
    OscillationHandler.cpp
    Gradient/LikelihoodGradient.cpp
    Gradient/GradientCheck.cpp
    Reweight/ReweightModel.cpp
    Reweight/EventGroups.cpp
    Reweight/CpuReweight.cpp
    Reweight/ReweightOptions.cpp
    Reweight/LinearShiftBuilder.cpp
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:Reweight/GPU/CudaReweight.cpp>
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:Reweight/GPU/CudaReweightDevice.cu>
    Oscillation/OscillationGradient.cpp
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:Oscillation/OscillationGradientDevice.cu>
)

target_link_libraries(Samples PUBLIC Splines NuOscillator)
if(MaCh3_GPU_ENABLED)
  set_property(TARGET Samples PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

#HW The Apple GPU backend: Objective-C++ against the Metal framework, compiled only when asked for
if(APPLE AND MaCh3_METAL_ENABLED)
  enable_language(OBJCXX)
  target_sources(Samples PRIVATE Reweight/GPU/MetalReweight.mm)
  set_source_files_properties(Reweight/GPU/MetalReweight.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
  target_link_libraries(Samples PRIVATE "-framework Metal" "-framework Foundation")
  target_compile_definitions(Samples PUBLIC MaCh3_METAL)
  cmessage(STATUS "Building the Metal reweight backend")
endif()
target_link_libraries(Samples PRIVATE MaCh3Warnings)

target_include_directories(Samples PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(Samples PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Samples)

install(TARGETS Samples
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Samples)
#HW The gradient sampling code, by topic. PUBLIC_HEADER would flatten these folders, so they are
#installed with their folders kept (Reweight/GPU included), which is how they are included
install(DIRECTORY Reweight/ DESTINATION include/Samples/Reweight FILES_MATCHING PATTERN "*.h")
install(DIRECTORY Gradient/ DESTINATION include/Samples/Gradient FILES_MATCHING PATTERN "*.h")
install(DIRECTORY Oscillation/ DESTINATION include/Samples/Oscillation FILES_MATCHING PATTERN "*.h")


add_library(MaCh3::Samples ALIAS Samples)
target_link_libraries(MaCh3 INTERFACE Samples)
