set(HEADERS
    IOUtils.h
    HistogramUtils.h
    BinningUtils.h
    PhysicsUtils.h
)

add_library(Utils SHARED
    IOUtils.cpp
    HistogramUtils.cpp
    BinningUtils.cpp
    PhysicsUtils.cpp
)

set_target_properties(Utils PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Utils
)

target_link_libraries(Utils PUBLIC MaCh3::Manager)
target_link_libraries(Utils PRIVATE MaCh3Warnings)

install(TARGETS Utils
    EXPORT MaCh3-targets
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/Utils
)

add_library(MaCh3::Utils ALIAS Utils)
target_link_libraries(MaCh3 INTERFACE Utils)

# Add new folders
add_subdirectory(Plotting)
