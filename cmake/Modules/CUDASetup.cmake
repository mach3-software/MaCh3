enable_language(CUDA)

cmessage("I have enabled CUDA!!!")

if(NOT DEFINED CUDA_SAMPLES)
  cmessage(FATAL_ERROR "When using CUDA, CUDA_SAMPLES must be defined to point to the CUDAToolkit samples directory (should contain common/helper_functions.h).")
endif()

find_package(CUDAToolkit)

#add_compile_options("$<$<COMPILE_LANGUAGE:CUDA>:-I${CUDA_SAMPLES}/common/inc>")

EXECUTE_PROCESS( COMMAND uname -m OUTPUT_VARIABLE OS_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT "x86_64 " STREQUAL "${OS_ARCH} ")
	cmessage(FATAL_ERROR "This build currently only support x86_64 target arches, determined the arch to be: ${OS_ARCH}")
endif()

EXECUTE_PROCESS( COMMAND ${CMAKE_SOURCE_DIR}/cmake/cudaver.sh --major OUTPUT_VARIABLE CUDA_MAJOR_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

cmessage(STATUS "CUDA_MAJOR_VERSION: ${CUDAToolkit_VERSION}")


add_compile_definitions(CUDA)


#  add_compile_definitions(DEBUG_DUMP DEBUG_CUDA_ND280 DEBUG)

#if(CUDAToolkit_VERSION GREATER_EQUAL 11)
#  add_compile_options( 
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-arch=sm_52>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_35,code=sm_35>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_52,code=sm_52>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_60,code=sm_60>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_61,code=sm_61>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_70,code=sm_70>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_75,code=sm_75>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_80,code=sm_80>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_86,code=sm_86>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_86,code=compute_86>"
#  	)
#elseif(CUDAToolkit_VERSION LESS_EQUAL 9)
#	add_compile_options(
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_30,code=sm_30>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_32,code=sm_32>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_35,code=sm_35>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_35,code=compute_35>"
#  	)
#else()
#	add_compile_options(
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_30,code=sm_30>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_32,code=sm_32>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_35,code=sm_35>"
#  	"$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_35,code=compute_35>"
 # 	)
#endif()
  
   add_compile_options(
  	"$<$<COMPILE_LANGUAGE:CUDA>:-prec-sqrt=false;-use_fast_math;-O3;-Werror;cross-execution-space-call;-w>"
  	"$<$<COMPILE_LANGUAGE:CUDA>:-Xptxas=-allow-expensive-optimizations=true;-Xptxas=-fmad=true;-Xptxas=-O3;>"
  	"$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fpic;-Xcompiler=-O3;-Xcompiler=-Wall;-Xcompiler=-Wextra;-Xcompiler=-Werror;-Xcompiler=-Wno-error=unused-parameter>"
  	)
