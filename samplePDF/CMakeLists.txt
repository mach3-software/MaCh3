set(HEADERS
    samplePDFBase.h
    samplePDFFDBase.h
    interfacePDFEbE.h
    Structs.h
    samplePDFInterface.h
    ShiftFunctors.h
    FDMCStruct.h
)

add_library(SamplePDF SHARED
    samplePDFBase.cpp
    samplePDFFDBase.cpp
    interfacePDFEbE.cpp
    Structs.cpp
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:probGpu.cu>
)

target_link_libraries(SamplePDF Splines NuOscillator Covariance MaCh3CompilerOptions)

target_include_directories(SamplePDF PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:$ENV{CUDAPATH}/include>>
    $<$<NOT:$<BOOL:${CPU_ONLY}>>:$<BUILD_INTERFACE:$ENV{CUDAPATH}/samples/common/inc>>
)

if(NOT CPU_ONLY)
    set_target_properties(SamplePDF PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET SamplePDF PROPERTY CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES})
    target_compile_definitions(SamplePDF PUBLIC GPU_ON)
endif()

set_target_properties(SamplePDF PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME SamplePDF)

install(TARGETS SamplePDF
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/samplePDF)

add_library(MaCh3::SamplePDF ALIAS SamplePDF)

install(TARGETS NuOscillator
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
		PUBLIC_HEADER DESTINATION include/NuOscillator
        )

		#install(TARGETS Oscillator
		#        EXPORT MaCh3-targets
		#	    LIBRARY DESTINATION lib/	
		#		PUBLIC_HEADER DESTINATION include/Oscillator)
		#
		install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/Oscillator/
		  DESTINATION include/Oscillator
		)
		#
		#install(TARGETS OscProbCalcer
		#        EXPORT MaCh3-targets
		#	    LIBRARY DESTINATION lib/	
		#		PUBLIC_HEADER DESTINATION include/OscProbCalcer)
		#
		install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/OscProbCalcer/
		  DESTINATION include/OscProbCalcer
		)

		install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/Constants/
		  DESTINATION include/Constants
		)

		#install(TARGETS NuOscillatorInterface
		#EXPORT MaCh3-targets
		#LIBRARY DESTINATION lib/	
		#PUBLIC_HEADER DESTINATION include/OscProbCalcer)
