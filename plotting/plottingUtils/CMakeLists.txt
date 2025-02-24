set( HEADERS
    plottingManager.h
    inputManager.h
    styleManager.h
    plottingUtils.h
)

add_library( Plotting SHARED
    plottingManager.cpp
    inputManager.cpp
    styleManager.cpp
    plottingUtils.cpp
)

set_target_properties( Plotting PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Plotting
)

target_link_libraries(Plotting PUBLIC MaCh3::Manager MaCh3::Fitter)
target_link_libraries(Plotting PRIVATE MaCh3Warnings)

## to be compiled into python module needs to be compiled as position independent library
if( MaCh3_PYTHON_ENABLED )
    set_property( TARGET Plotting PROPERTY POSITION_INDEPENDENT_CODE ON )
endif()

install(TARGETS Plotting
    EXPORT MaCh3-targets
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/plotting
)

add_library(MaCh3::Plotting ALIAS Plotting)
