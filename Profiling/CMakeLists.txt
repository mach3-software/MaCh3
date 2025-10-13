
set(HEADERS
    Profiling.h
)

add_library(Profiling SHARED
    ${tracy_SOURCE_DIR}/public/TracyClient.cpp
)

set_target_properties(Profiling PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Profiling
)

if(MaCh3_SAMPLED_PROFILING_ENABLED OR MaCh3_INSTRUMENTED_PROFILING_ENABLED)
    target_compile_definitions(Profiling PUBLIC TRACY_ENABLE)

    if(MaCh3_INSTRUMENTED_PROFILING_ENABLED)
        target_compile_definitions(Profiling PUBLIC MaCh3_INSTRUMENTED_PROFILING_ENABLED)
    endif()

    set_property( TARGET Profiling PROPERTY POSITION_INDEPENDENT_CODE ON )
    target_link_libraries(Profiling PUBLIC Tracy::TracyClient)

    install(TARGETS TracyClient
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/tracy
    )

endif()

install(TARGETS Profiling
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Profiling)

add_library(MaCh3::Profiling ALIAS Profiling)
