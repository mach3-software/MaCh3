set(HEADERS
    BinningHandler.h
    SampleHandlerInterface.h
    SampleHandlerBase.h
    SampleStructs.h
    HistogramUtils.h
    SampleInfo.h
    EventInfo.h
    OscillationHandler.h
)

add_library(Samples SHARED
    BinningHandler.cpp
    SampleHandlerInterface.cpp
    SampleHandlerBase.cpp
    HistogramUtils.cpp
    OscillationHandler.cpp
)

target_link_libraries(Samples PUBLIC Splines NuOscillator)
target_link_libraries(Samples PRIVATE MaCh3Warnings)

target_include_directories(Samples PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(Samples PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Samples)

install(TARGETS Samples
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Samples)

install(DIRECTORY detail DESTINATION include/Samples)

add_library(MaCh3::Samples ALIAS Samples)
target_link_libraries(MaCh3 INTERFACE Samples)
