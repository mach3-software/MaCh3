################################## Oscillation ################################

# Oscillation calculation
DefineEnabledRequiredSwitch(CUDAProb3Linear_ENABLED FALSE)
DefineEnabledRequiredSwitch(CUDAProb3_ENABLED FALSE)
DefineEnabledRequiredSwitch(ProbGPULinear_ENABLED FALSE)
DefineEnabledRequiredSwitch(Prob3ppLinear_ENABLED FALSE)
DefineEnabledRequiredSwitch(NuFastLinear_ENABLED FALSE)
DefineEnabledRequiredSwitch(OscProb_ENABLED FALSE)

#KS: If all Oscillators are turned off then enable CUDAProb3Linear_ENABLED
if (NOT CUDAProb3Linear_ENABLED AND
    NOT CUDAProb3_ENABLED AND
    NOT ProbGPULinear_ENABLED AND
    NOT Prob3ppLinear_ENABLED AND
    NOT NuFastLinear_ENABLED AND
    NOT OscProb_ENABLED)
    set(NuFastLinear_ENABLED TRUE)
endif()

#KS: Save which oscillators are being used
set(MaCh3_Oscillator_ENABLED "")
if(CUDAProb3Linear_ENABLED)
  LIST(APPEND MaCh3_Oscillator_ENABLED "CUDAProb3Linear")
endif()
if(CUDAProb3_ENABLED)
  LIST(APPEND MaCh3_Oscillator_ENABLED "CUDAProb3")
endif()
if(ProbGPULinear_ENABLED)
  LIST(APPEND MaCh3_Oscillator_ENABLED "ProbGPULinear")
endif()
if(Prob3ppLinear_ENABLED)
  LIST(APPEND MaCh3_Oscillator_ENABLED "Prob3ppLinear")
endif()
if(NuFastLinear_ENABLED)
  LIST(APPEND MaCh3_Oscillator_ENABLED "NuFast")
endif()
if(OscProb_ENABLED)
  LIST(APPEND MaCh3_Oscillator_ENABLED "OscProb")
endif()

#NuOscillator uses 1/0 instead of true/false thus use conversion
IsTrue(CUDAProb3Linear_ENABLED USE_CUDAProb3Linear)
IsTrue(CUDAProb3_ENABLED USE_CUDAProb3)
IsTrue(ProbGPULinear_ENABLED USE_ProbGPULinear)
IsTrue(Prob3ppLinear_ENABLED USE_Prob3ppLinear)
IsTrue(NuFastLinear_ENABLED USE_NuFastLiner)
IsTrue(OscProb_ENABLED USE_OscProb)

#Also additional flags
IsTrue(MaCh3_GPU_ENABLED DAN_USE_GPU)
IsTrue(MaCh3_MULTITHREAD_ENABLED DAN_USE_MULTITHREAD)
IsTrue(MaCh3_LOW_MEMORY_STRUCTS_ENABLED DAN_DOUBLE)
SwitchLogic(DAN_DOUBLE)

# Get the CPU compile options for MaCh3CompilerOptions
get_target_property(cpu_compile_options MaCh3CompilerOptions INTERFACE_COMPILE_OPTIONS)

# Join the compile options list into a space-separated string
string(REPLACE ";" " " cpu_compile_options_string "${cpu_compile_options}")


#KS: This may seem hacky, but when CMAKE_CUDA_ARCHITECTURES is passed, it's treated as a string rather than a list. Since CMake uses semi-colon-delimited strings to represent lists, we convert it to a proper list to handle CUDA architectures correctly.
set(CMAKE_CUDA_ARCHITECTURES_STRING ${CMAKE_CUDA_ARCHITECTURES})
#string(REPLACE " " ";" CMAKE_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
string(REPLACE " " ";" CMAKE_CUDA_ARCHITECTURES_STRING "${CMAKE_CUDA_ARCHITECTURES}")

#Try adding Oscillator Class
CPMAddPackage(
  NAME NuOscillator
    VERSION 1.1.0
    GITHUB_REPOSITORY "dbarrow257/NuOscillator"
    GIT_TAG "dbarrow257/feature/BetterThrowing"
    GIT_SHALLOW YES
    OPTIONS
    "UseGPU ${DAN_USE_GPU}"
    "UseMultithreading ${DAN_USE_MULTITHREAD}"
    "UseDoubles ${DAN_DOUBLE}"

    "UseCUDAProb3Linear ${USE_CUDAProb3Linear}"
    "UseCUDAProb3 ${USE_CUDAProb3}"
    "UseProbGPULinear ${USE_ProbGPULinear}"
    "UseProb3ppLinear ${USE_Prob3ppLinear}"
    "UseNuFASTLinear  ${USE_NuFastLiner}"
    "UseOscProb ${USE_OscProb}"

    "NuOscillator_Compiler_Flags ${cpu_compile_options_string}"
    "CMAKE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES_STRING}"
    "CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD}"
)

if(NOT TARGET NuOscillator)
  cmessage(FATAL_ERROR "Expecting dependency NuOscillator")
endif()

install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/Oscillator/
        DESTINATION include/Oscillator)

install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/OscProbCalcer/
        DESTINATION include/OscProbCalcer)

install(DIRECTORY ${NUOSCILLATOR_SOURCE_DIR}/Constants/
        DESTINATION include/Constants)
