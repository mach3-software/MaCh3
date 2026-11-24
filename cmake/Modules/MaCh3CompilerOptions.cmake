# KS: Compile and link options for more see https://github.com/cpp-best-practices/cppbestpractices/tree/master
add_library(MaCh3Warnings INTERFACE)
target_compile_options(MaCh3Warnings INTERFACE
  -Wextra                 # Enable extra warning flags
  -Wall                   # Enable all standard warning flags
  -Wshadow                # Warn when a variable declaration shadows one from an outer scope
  -Wuninitialized         # Warn about uninitialized variables
  -Wnon-virtual-dtor      # Warn when a class with virtual functions has a non-virtual destructor
  -Woverloaded-virtual    # Warn when a function declaration hides a virtual function from a base class
  -Wformat=2              # Warn on security issues around functions that format output (ie printf)
  -Wunused                # Warn on anything being unused
  -Wredundant-decls       # Warn about multiple declarations of the same entity. Useful for code cleanup.
  -Wstrict-aliasing=2     # Helps detect potential aliasing issues that could lead to undefined behavior.
  -Wnull-dereference      # Warn if a null dereference is detected (only in GCC >= 6.0)
  -Wold-style-cast        # Warn for c-style casts
  -Wconversion            # Warn on type conversions that may lose data
  -Wformat-security       # Warn on functions that are potentially insecure for formatting
  -Walloca                # Warn if `alloca` is used, as it can lead to stack overflows
  -Wswitch-enum           # Warn if a `switch` statement on an enum does not cover all values
  -pedantic               # Enforce strict ISO compliance (all versions of GCC, Clang >= 3.2)
  -Wcast-align            # Warn when a pointer cast may result in stricter alignment requirements
  #-Wdouble-promotion     # Warn when float values are implicitly promoted to double
  #-Wfloat-equal          # Warn if floating-point values are compared directly
  #-Wpadded               # Warn when padding is added to a structure or class for alignment
)
# KS Some compiler options are only available in GCC, in case we move to other compilers we will have to expand this
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(MaCh3Warnings INTERFACE
    -Wlogical-op            # Warn about logical operations being used where bitwise were probably wanted (only in GCC)
    -Wduplicated-cond       # Warn if if / else chain has duplicated conditions (only in GCC >= 6.0)
    -Wduplicated-branches   # Warn if if / else branches have duplicated code (only in GCC >= 7.0)
    -Wuseless-cast          # Warn if you perform a cast to the same type (only in GCC >= 4.8)
    )
  # KS: For older gcc it complaints that we should have override final, while we follow modern approach and have just final
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 10)
    target_compile_options(MaCh3Warnings INTERFACE
      -Wsuggest-override   # Warn if a virtual function override is not marked 'override' (only in GCC >= 5.1)
    )
  endif()
endif()

# KS: Clang and IntelLLVM is super picky, it is almost impossible to make it work with Werror
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  set(MaCh3_WERROR_ENABLED FALSE)
endif()

if(MaCh3_WERROR_ENABLED)
target_compile_options(MaCh3Warnings INTERFACE
  -Werror                # Treat Warnings as Errors
)
endif()

add_library(MaCh3CompilerOptions INTERFACE)
target_link_libraries(MaCh3CompilerOptions INTERFACE MaCh3CompileDefinitions)
target_compile_options(MaCh3CompilerOptions INTERFACE
    -g                   # Generate debug information
)

#If DEBUG_LEVEL was defined but MaCh3_DEBUG_ENABLED not, enable debug flag
if(DEFINED DEBUG_LEVEL)
  set(MaCh3_DEBUG_ENABLED TRUE)
else()
  #If MaCh3 debug was enable but level not, set it to 1. In very rare cases we want to go beyond 1.
  if(MaCh3_DEBUG_ENABLED)
    set(DEBUG_LEVEL 1)
  endif()
endif()

#KS: If Debug add debugging compile flag if not add optimisation for speed
if(MaCh3_DEBUG_ENABLED)
  target_compile_options(MaCh3CompilerOptions INTERFACE
    -O0                     # Turn off any optimisation to have best debug experience
    -fno-omit-frame-pointer # Keep frame pointers: essential for accurate stack traces
    #-fsanitize=address     # Enable AddressSanitizer
  )
  # KS: Disable AVX only on x86 platforms
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i.86")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      target_compile_options(MaCh3CompilerOptions INTERFACE
          -mno-avx  # Disable AVX instructions to simplify debugging and improve reproducibility
      )
    endif()
  endif()
  target_compile_definitions(MaCh3CompileDefinitions INTERFACE MACH3_DEBUG=${DEBUG_LEVEL})
  cmessage(STATUS "Enabling DEBUG with Level: \"${DEBUG_LEVEL}\"")
else()
  #https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
  target_compile_options(MaCh3CompilerOptions INTERFACE
    -O3                                   # Optimize code for maximum speed
    -funroll-loops                        # Unroll loops where possible for performance
    -fno-math-errno                       # Don’t set errno after math calls (can speed up hot loops with math)
    # KS: Below could help with performance but require validations and benchmarking
    #-fdevirtualize-at-ltrans             # Perform additional virtual function devirtualisation during link-time optimisation
    #-fno-signed-zeros                    # Allow optimisations that don't preserve the distinction between +0.0 and -0.0
    #-fno-trapping-math                   # Assume FP ops don’t trap (helps vectorisation)
    #--param=max-vartrack-size=100000000  # Set maximum size of variable tracking data to avoid excessive memory usage
  )
  # KS Some compiler options are only available in GCC, in case we move to other compilers we will have to expand this
  #KS: Consider in future __attribute__((always_inline)) see https://indico.cern.ch/event/386232/sessions/159923/attachments/771039/1057534/always_inline_performance.pdf
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(MaCh3CompilerOptions INTERFACE
      $<$<COMPILE_LANGUAGE:CXX>:-flto=auto>   # Enable link-time optimization, don't apply it to CUDA or other languages!
      -finline-limit=100000000                # Increase the limit for inlining functions to improve performance
      #-fipa-pta                              # Pointer analysis for better intraprocedural optimization, increases compilation time
    )
  #KS: add Link-Time Optimization (LTO)
  target_link_libraries(MaCh3CompilerOptions INTERFACE -flto=auto)
  endif()
endif()

if(MaCh3_NATIVE_ENABLED)
  target_compile_options(MaCh3CompilerOptions INTERFACE
    -march=native    #Generate instructions for the host CPU architecture.
    -mtune=native    #Optimize instruction scheduling and pipeline usage for host CPU.
  )
endif()

#Add Multithread flags
if(MaCh3_MULTITHREAD_ENABLED)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    #ETA: intel icpx requires this I believe
    find_package(OpenMP REQUIRED)
    target_compile_options(MaCh3CompilerOptions INTERFACE -qopenmp)
    target_link_libraries(MaCh3CompilerOptions INTERFACE OpenMP::OpenMP_CXX)
  # KS: Not a clue why clang need different...
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(MaCh3CompilerOptions INTERFACE -fopenmp)
    target_link_libraries(MaCh3CompilerOptions INTERFACE omp)
  else()
    target_compile_options(MaCh3CompilerOptions INTERFACE -fopenmp)
    target_link_libraries(MaCh3CompilerOptions INTERFACE gomp)
  endif()
  target_compile_definitions(MaCh3CompileDefinitions INTERFACE MULTITHREAD)
  # KS: This will pass variable to doxygen and tell whether documentation should be with MP functions
  set(DOXYGEN_PREDEFINED "MULTITHREAD")
else()
  set(DOXYGEN_PREDEFINED "")
endif()

if(MaCh3_GPU_ENABLED)
  target_compile_definitions(MaCh3GPUCompilerOptions INTERFACE CUDA)
  target_compile_definitions(MaCh3GPUCompilerOptions INTERFACE GPU_ON)
endif()

if(MaCh3_LOW_MEMORY_STRUCTS_ENABLED)
  target_compile_definitions(MaCh3CompileDefinitions INTERFACE _LOW_MEMORY_STRUCTS_)
endif()
