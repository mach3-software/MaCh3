################################## Oscillation ################################
# KS: Define the list of oscillator options
set(OSCILLATOR_OPTIONS
    CUDAProb3Linear
    CUDAProb3
    ProbGPULinear
    Prob3ppLinear
    NuFASTLinear
    NuSQUIDSLinear
    OscProb
    GLoBESLinear
)
# KS: Tells whether all oscillators were disabled
set(ALL_OSCILLATORS_DISABLED TRUE)

# KS: Oscillation calculation
foreach(option IN LISTS OSCILLATOR_OPTIONS)
    DefineEnabledRequiredSwitch(${option}_ENABLED FALSE)
    # If at least one oscillator is enabled then set ALL_OSCILLATORS_DISABLED to false
    if(${option}_ENABLED)
      set(ALL_OSCILLATORS_DISABLED FALSE)
  endif()
endforeach()

#KS: If all Oscillators are turned off then enable NuFastLinear_ENABLED and CUDAProb3_ENABLED
if(ALL_OSCILLATORS_DISABLED)
    set(NuFASTLinear_ENABLED TRUE)
    set(CUDAProb3_ENABLED TRUE)
endif()

#KS: Save which oscillators are being used
set(MaCh3_Oscillator_ENABLED "")
foreach(option IN LISTS OSCILLATOR_OPTIONS)
  if(${option}_ENABLED)
    LIST(APPEND MaCh3_Oscillator_ENABLED ${option})
  endif()

  #NuOscillator uses 1/0 instead of true/false thus use conversion
  IsTrue(${option}_ENABLED USE_${option})
endforeach()

# Prepare the options for CPMAddPackage
set(NUOSCILLATOR_OPTIONS "")
foreach(option IN LISTS OSCILLATOR_OPTIONS)
    list(APPEND NUOSCILLATOR_OPTIONS "Use${option} ${USE_${option}}")
endforeach()

#Also additional flags
IsTrue(MaCh3_GPU_ENABLED DAN_USE_GPU)
IsTrue(MaCh3_MULTITHREAD_ENABLED DAN_USE_MULTITHREAD)
IsTrue(MaCh3_LOW_MEMORY_STRUCTS_ENABLED DAN_DOUBLE)
SwitchLogic(DAN_DOUBLE)

# Diable GPU at NuOsc
if(NOT MaCh3_NuOsc_GPU_ENABLED)
  set(DAN_USE_GPU 0)
endif()
# Get the CPU compile options for MaCh3CompilerOptions
get_target_property(cpu_compile_options MaCh3CompilerOptions INTERFACE_COMPILE_OPTIONS)

# Join the compile options list into a space-separated string
string(REPLACE ";" " " cpu_compile_options_string "${cpu_compile_options}")

#KS: This may seem hacky, but when CMAKE_CUDA_ARCHITECTURES is passed, it's treated as a string rather than a list. Since CMake uses semi-colon-delimited strings to represent lists, we convert it to a proper list to handle CUDA architectures correctly.
set(CMAKE_CUDA_ARCHITECTURES_STRING ${CMAKE_CUDA_ARCHITECTURES})
string(REPLACE " " ";" CMAKE_CUDA_ARCHITECTURES_STRING "${CMAKE_CUDA_ARCHITECTURES}")

if(NOT DEFINED MaCh3_NuOscillatorBranch)
  set(MaCh3_NuOscillatorBranch "v1.4.4")
endif()

#Try adding Oscillator Class
CPMAddPackage(
  NAME NuOscillator
    GITHUB_REPOSITORY "dbarrow257/NuOscillator"
    GIT_TAG ${MaCh3_NuOscillatorBranch}
    GIT_SHALLOW YES
    OPTIONS
    "UseGPU ${DAN_USE_GPU}"
    "UseMultithreading ${DAN_USE_MULTITHREAD}"
    "UseDoubles ${DAN_DOUBLE}"

    ${NUOSCILLATOR_OPTIONS}

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
