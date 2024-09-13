################################## Oscillation ################################

# Oscillation calculation
# In the future which osc calc we use might be set with a flag
set(MaCh3_Oscillator_ENABLED "")
LIST(APPEND MaCh3_Oscillator_ENABLED "CUDAProb3Linear")

set(DAN_USE_GPU 0)
set(DAN_USE_MULTITHREAD 0)
set(DAN_DOUBLE 0)

IsTrue(MaCh3_GPU_ENABLED DAN_USE_GPU)
IsTrue(MaCh3_MULTITHREAD_ENABLED DAN_USE_MULTITHREAD)

IsTrue(MaCh3_LOW_MEMORY_STRUCTS_ENABLED DAN_DOUBLE)
SwitchLogic(DAN_DOUBLE)

# Get the compile options for MaCh3CompilerOptions
get_target_property(compile_options MaCh3CompilerOptions INTERFACE_COMPILE_OPTIONS)

# Join the compile options list into a space-separated string
string(REPLACE ";" " " compile_options_string "${compile_options}")

#KS this may look hacky however CPM isn't build for passing stuff like this. If CMAKE_CUDA_ARCHITECTURES is passed CPM it will be string not list. Thus we convert it to list
set(CMAKE_CUDA_ARCHITECTURES_STRING ${CMAKE_CUDA_ARCHITECTURES})
string(REPLACE " " ";" CMAKE_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")

#Try adding Oscillator Class
CPMAddPackage(
  NAME NuOscillator
    VERSION 0.0
    GITHUB_REPOSITORY "dbarrow257/NuOscillator"
    GIT_TAG "feature_Cmake"
    OPTIONS
    "UseGPU ${DAN_USE_GPU}"
    "UseMultithreading ${DAN_USE_MULTITHREAD}"
    "UseDoubles ${DAN_DOUBLE}"

    "UseCUDAProb3 0"
    "UseCUDAProb3Linear 1"
    "UseProbGPULinear 0"
    "UseProb3ppLinear 0"

    "NuOscillator_Compiler_Flags ${compile_options_string}"
    "CMAKE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES_STRING}"
)

if(NOT TARGET NuOscillator)
  cmessage(FATAL_ERROR "Expecting dependency NuOscillator")
endif()

install(TARGETS NuOscillator
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/NuOscillator)

install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/Oscillator/
        DESTINATION include/Oscillator)

install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/OscProbCalcer/
        DESTINATION include/OscProbCalcer)

install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/Constants/
        DESTINATION include/Constants)
