
## build the base level library targets
foreach(lib
    inputManager
    styleManager
    plottingUtils
    )
    add_library( ${lib} ${lib}.cpp )
    target_link_libraries( ${lib} ${ROOT_LIBRARIES} MaCh3::Manager )

    set_property( TARGET ${lib} PROPERTY POSITION_INDEPENDENT_CODE ON )

endforeach(lib)


## build the higher level plottingManager
add_library( plottingManager plottingManager.cpp plottingManager.h )
set_property( TARGET plottingManager PROPERTY POSITION_INDEPENDENT_CODE ON )
target_link_libraries( plottingManager inputManager styleManager )


## chuck it all together into the Plotting library
add_library( Plotting INTERFACE ) 
target_link_libraries( Plotting INTERFACE inputManager plottingManager styleManager plottingUtils )