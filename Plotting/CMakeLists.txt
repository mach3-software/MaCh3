add_subdirectory(plottingUtils)

foreach(app
    GetPostfitParamPlots
    PredictivePlotting
    PlotLLH
    MatrixPlotter
    PlotSigmaVariation
    )
    add_executable(${app} ${app}.cpp )
    target_link_libraries(${app} MaCh3::Plotting MaCh3Warnings)
    set_target_properties(${app} PROPERTIES
        BUILD_RPATH "${CMAKE_BINARY_DIR}/lib"
        INSTALL_RPATH "@loader_path/../lib"
    )

    install(TARGETS ${app} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach(app)
