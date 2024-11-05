set(HEADERS
    samplePDFBase.h
    samplePDFFDBase.h
    Structs.h
    ShiftFunctors.h
    FDMCStruct.h
)

if (CUDA)
    add_library(SamplePDF SHARED
        samplePDFBase.cpp
        samplePDFFDBase.cpp
        Structs.cpp
        $<$<NOT:$<BOOL:${CPU_ONLY}>>:probGpu.cu>
        $<$<NOT:$<BOOL:${CPU_ONLY}>>:kdeGpu.cu>
        $<$<NOT:$<BOOL:${CPU_ONLY}>>:kdeGpu2.cu>
    )
    set_target_properties(SamplePDF PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET SamplePDF PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})

else()
    add_library(SamplePDF SHARED
        samplePDFBase.cpp
        samplePDFFDBase.cpp
        Structs.cpp
        $<$<NOT:$<BOOL:${CPU_ONLY}>>:probGpu.cu>
        $<$<NOT:$<BOOL:${CPU_ONLY}>>:kdeGpu.cu>
        $<$<NOT:$<BOOL:${CPU_ONLY}>>:kdeGpu2.cu>
    )

endif()

target_link_libraries(SamplePDF Splines Covariance MaCh3CompilerOptions)

if(USE_PROB3)
    target_link_libraries(SamplePDF Prob3plusplus)
else()
    target_link_libraries(SamplePDF OscClass)
endif()

target_include_directories(SamplePDF PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(SamplePDF PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME SamplePDF)

install(TARGETS SamplePDF
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/samplePDF)

add_library(MaCh3::SamplePDF ALIAS SamplePDF)


