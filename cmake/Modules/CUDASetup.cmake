enable_language(CUDA)
add_compile_definitions(CUDA)

#KS: If Debug is not defined disable it by default
if(NOT DEFINED MaCh3_DEBUG_ENABLED)
  SET(MaCh3_DEBUG_ENABLED FALSE)
endif()

EXECUTE_PROCESS( COMMAND uname -m OUTPUT_VARIABLE OS_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT "x86_64 " STREQUAL "${OS_ARCH} ")
    cmessage(FATAL_ERROR "This build currently only support x86_64 target arches, determined the arch to be: ${OS_ARCH}")
endif()

cmessage(STATUS "CUDA_MAJOR_VERSION: ${CUDAToolkit_VERSION}")

# Call nvidia-smi using execute_process
execute_process(
    COMMAND nvidia-smi
    RESULT_VARIABLE NVSMI_RESULT
)

# Check the result of nvidia-smi command
if(NVSMI_RESULT EQUAL 0)
    cmessage(STATUS "nvidia-smi command executed successfully.")
else()
    cmessage(WARNING "Failed to execute nvidia-smi command.")
endif()

#KS: Allow user to define CMAKE_CUDA_ARCHITECTURES
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  #KS: See this for more info https://cmake.org/cmake/help/latest/prop_tgt/CUDA_ARCHITECTURES.html
  if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.23 )
    set(CMAKE_CUDA_ARCHITECTURES all )
    #KS: Consider using native, requires cmake 3.24... will be terrible for containers but should results in more optimised code
    #set(CMAKE_CUDA_ARCHITECTURES native )
  else()
    #KS: Apparently with newer cmake and GPU
    set(CMAKE_CUDA_ARCHITECTURES 35 52 60 61 70 75 80 86)
  endif()
endif()

include_directories(
    $ENV{CUDAPATH}/include
    $ENV{CUDAPATH}/samples/common/inc>
)

# Join elements of the list with spaces
string(REPLACE ";" " " CUDA_ARCHITECTURES_STR "${CMAKE_CUDA_ARCHITECTURES}")

# Output the message with spaces between numbers
cmessage(STATUS "Using following CUDA architectures: ${CUDA_ARCHITECTURES_STR}")

if(NOT MaCh3_DEBUG_ENABLED)
    target_compile_options(MaCh3CompilerOptions INTERFACE
        "$<$<COMPILE_LANGUAGE:CUDA>:-prec-sqrt=false;-use_fast_math;-O3;-Werror;cross-execution-space-call;-w>"
        "$<$<COMPILE_LANGUAGE:CUDA>:-Xptxas=-allow-expensive-optimizations=true;-Xptxas=-fmad=true;-Xptxas=-O3;>"
        "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fpic;-Xcompiler=-O3;-Xcompiler=-Wall;-Xcompiler=-Wextra;-Xcompiler=-Werror;-Xcompiler=-Wno-error=unused-parameter>"
    )
else()
#CWret: -g and -G for debug flags to use cuda-gdb; slows stuff A LOT
#-pxtas-options=-v, -maxregcount=N
    target_compile_options(MaCh3CompilerOptions INTERFACE
        "$<$<COMPILE_LANGUAGE:CUDA>:-prec-sqrt=false;-use_fast_math;-Werror;cross-execution-space-call;-w>"
        "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-g;>"
        "$<$<COMPILE_LANGUAGE:CUDA>:-Xptxas=-dlcm=ca;-Xptxas=-warn-lmem-usage;-Xptxas=-warn-spills;-Xptxas=-v;-Xcompiler=-Wall;-Xcompiler=-Wextra;-Xcompiler=-Werror;-Xcompiler=-Wno-error=unused-parameter>"
    )

    target_compile_definitions(MaCh3CompilerOptions INTERFACE "$<$<COMPILE_LANGUAGE:CUDA>:CUDA_ERROR_CHECK>")
endif()



target_link_options(MaCh3CompilerOptions INTERFACE -I$ENV{CUDAPATH}/lib64 -I$ENV{CUDAPATH}/include -I$ENV{CUDAPATH}/common/inc -I$ENV{CUDAPATH}/samples/common/inc)

if(NOT DEFINED ENV{CUDAPATH})
    cmessage(FATAL_ERROR "CUDAPATH environment variable is not defined. Please set it to the root directory of your CUDA installation.")
endif()

#if(NOT DEFINED CUDA_SAMPLES)
#  cmessage(FATAL_ERROR "When using CUDA, CUDA_SAMPLES must be defined to point to the CUDAToolkit samples directory (should contain common/helper_functions.h).")
#endif()

#HW: If we're doing GPU stuff, we need the CUDA helper module
if(BUILTIN_CUDASAMPLES_ENABLED )
    cmessage("Adding Builtin CUDA Sample Library")
    CPMAddPackage(
        NAME cuda-samples
        GITHUB_REPOSITORY "NVIDIA/cuda-samples"
        GIT_TAG v12.3
    SOURCE_SUBDIR Common
    DOWNLOAD_ONLY YES
    )
# Now we add the library
    if(cuda-samples_ADDED)
         add_library(cuda-samples INTERFACE IMPORTED)
         target_include_directories(cuda-samples INTERFACE ${cuda-samples_SOURCE_DIR})
    endif()
endif()

#KS: Keep this for backward compatibility

#-Xptxas "-dlcm=ca -v --allow-expensive-optimizations=true -fmad=true"
#-Xcompiler "-fpic" -c
#-prec-sqrt=false -use_fast_math

#CWRET comment: could change arch here. also, prec-div is only used in (seg+khig)/2; could replace by *0.5
# -prec-sqrt is not needed (no sqrts in program code!), -use_fast_math forces intrinsic
#  tried this on GTX 660 and no speed increase, buu and P6 6c
#  don't use fastmath!

#-Xptxas -dlcm=ca increases L1 cache from 16kB to 48kB for Fermi
#-Xptxas -dlcm=cg caches global memory into shared memory and turns off the L1 cache entirely
#

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
  
