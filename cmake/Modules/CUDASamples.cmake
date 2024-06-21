# Variable to hold the paths once found
set(CMAKE_CUDA_SAMPLES_PATH "")
set(CUDASAMPLES_FOUND FALSE)

#KS: List of directories to search for helper_functions.h and cuda_runtime.h
# Please expand if needed
set(CUDA_SAMPLE_SEARCH_PATHS
    "${CUDAToolkit_TARGET_DIR}/samples/common/inc/"
)

# Loop over the directories and append to CMAKE_CUDA_SAMPLES_PATH if headers are found
foreach(dir ${CUDA_SAMPLE_SEARCH_PATHS})
    if(EXISTS "${dir}/helper_functions.h")
        list(APPEND CMAKE_CUDA_SAMPLES_PATH ${dir})
        set(CUDASAMPLES_FOUND TRUE)
    endif()
endforeach()

#HW: If we're doing GPU stuff, we need the CUDA helper module
if(NOT CUDASAMPLES_FOUND)
    cmessage(WARNING "Adding Built-in CUDA Sample Library, it make take some time")
    CPMAddPackage(
        NAME cuda-samples
        GITHUB_REPOSITORY "NVIDIA/cuda-samples"
        GIT_TAG v12.3
    SOURCE_SUBDIR Common
    DOWNLOAD_ONLY YES
    )
    list(APPEND CMAKE_CUDA_SAMPLES_PATH ${cuda-samples_SOURCE_DIR}/Common)
endif()

cmessage(STATUS "Using the following CUDA samples paths: ${CMAKE_CUDA_SAMPLES_PATH}")
target_include_directories(MaCh3CompilerOptions INTERFACE ${CMAKE_CUDA_SAMPLES_PATH})
