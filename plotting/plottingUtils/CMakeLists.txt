
## build the base level library targets
foreach(lib
    inputManager
    styleManager
    plottingUtils
    )
    add_library( ${lib} ${lib}.cpp )
    target_link_libraries( ${lib} ${ROOT_LIBRARIES} MaCh3::Manager )

endforeach(lib)


## build the higher level plottingManager
add_library( plottingManager plottingManager.cpp plottingManager.h )
target_link_libraries( plottingManager inputManager styleManager )


## chuck it all together into the Plotting library
add_library( Plotting INTERFACE ) 
target_link_libraries( Plotting INTERFACE inputManager plottingManager styleManager plottingUtils )
