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
target_include_directories(MaCh3GPUCompilerOptions INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CUDA_SAMPLES_PATH}>"
    "$<INSTALL_INTERFACE:include>"
)

# KS: Perform fancy CUDA Benchmarking
DefineEnabledRequiredSwitch(MaCh3_GPU_BENCHMARK FALSE)
if(MaCh3_GPU_BENCHMARK)
    cmessage(STATUS "Building CUDA Benchmark")

    # KS: Define directories to iterate over, might be useful to expand
    set(CUDA_SAMPLES_DIRS
        "deviceQuery"
        "bandwidthTest"
    )

    # KS: Iterate over each directory
    foreach(sample_dir ${CUDA_SAMPLES_DIRS})
        # Define source and destination directories
        set(SRC_DIR "${CMAKE_CUDA_SAMPLES_PATH}/../Samples/1_Utilities/${sample_dir}")
        set(DST_DIR "${CMAKE_BINARY_DIR}/GPU_Benchmark/")

        # CW: Copy over the provided nvidia utility
        # CW: Often we can't write to the CUDA install directory, so let's build it here
        file(COPY ${SRC_DIR} DESTINATION ${DST_DIR})

        # KS: Change directory to copied sample
        set(SAMPLE_DIR "${CMAKE_BINARY_DIR}/GPU_Benchmark/${sample_dir}")

        # Modify Makefile path
        set(MAKEFILE_PATH "${SAMPLE_DIR}/Makefile")

        # CW: Patch the litle hard-coded NVIDIA makefile
        execute_process(
            COMMAND sed -i "s,../../../Common,${CMAKE_CUDA_SAMPLES_PATH},g" ${MAKEFILE_PATH}
            RESULT_VARIABLE SED_RESULT
        )

        # Add custom target to run make
        add_custom_target(run_${sample_dir} ALL
            COMMAND make
            WORKING_DIRECTORY ${SAMPLE_DIR}
        )

        # Add custom target to run sample
        add_custom_target(run_${sample_dir}_exec ALL
            COMMAND ./${sample_dir}
            WORKING_DIRECTORY ${SAMPLE_DIR}
            DEPENDS run_${sample_dir}
        )
    endforeach(sample_dir)
endif()
