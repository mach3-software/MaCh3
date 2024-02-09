set(HEADERS OscClass_CUDAProb3.h
)

add_library(OscClass SHARED
    OscClass_CUDAProb3.cpp
)

set_target_properties(OscClass PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME OscClass)


if (USE_PROB3)
   message(FATAL_ERROR "OscClass only defined using CUDAProb3. USE_PROB3 is defined and indicates Prob3++/ProbGPU is being expected to be used.")
else()
   target_link_libraries(OscClass Manager CUDAProb3Beam CUDAProb3Atmos)
endif()
target_include_directories(OscClass PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS CUDAProb3Atmos
        EXPORT mach3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/CUDAProb3Atmos
        )

install(TARGETS CUDAProb3Beam
        EXPORT mach3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/CUDAProb3Beam
        )

install(TARGETS OscClass
        EXPORT mach3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/OscClass)

add_library(MACH3::OscClass ALIAS OscClass)
