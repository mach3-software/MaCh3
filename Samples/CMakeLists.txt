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
    LikelihoodGradient.h
    GradientCheck.h
    ReweightModel.h
    EventGroups.h
    ReweightBackend.h
    CpuReweight.h
    FastEventAccess.h
    LinearShifts.h
    LinearShiftBuilder.h
    ReweightOptions.h
    WeightParameterSnapshot.h
    OscillationProbability.h
    OscillationGradient.h
    MetalReweight.h
    CudaReweight.h
    CudaReweightDevice.h
    OscillationGradientDevice.h
)

add_library(Samples SHARED
    BinningHandler.cpp
    SampleHandlerInterface.cpp
    SampleHandlerBase.cpp
    OscillationHandler.cpp
    LikelihoodGradient.cpp
    GradientCheck.cpp
    ReweightModel.cpp
    EventGroups.cpp
    CpuReweight.cpp
    SampleHandlerReweight.cpp
    ReweightOptions.cpp
    LinearShiftBuilder.cpp
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:CudaReweight.cpp>
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:CudaReweightDevice.cu>
    $<$<BOOL:${MaCh3_GPU_ENABLED}>:OscillationGradientDevice.cu>
    OscillationGradient.cpp
)

target_link_libraries(Samples PUBLIC Splines NuOscillator)
if(MaCh3_GPU_ENABLED)
  set_property(TARGET Samples PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

#HW The Apple GPU backend: Objective-C++ against the Metal framework, compiled only when asked for
if(APPLE AND MaCh3_METAL_ENABLED)
  enable_language(OBJCXX)
  target_sources(Samples PRIVATE MetalReweight.mm)
  set_source_files_properties(MetalReweight.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
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


add_library(MaCh3::Samples ALIAS Samples)
target_link_libraries(MaCh3 INTERFACE Samples)
