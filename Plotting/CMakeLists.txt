add_subdirectory(PlottingUtils)

foreach(app
    GetPostfitParamPlots
    PredictivePlotting
    PlotLLH
    MatrixPlotter
    PlotSigmaVariation
    )
    add_executable(${app} ${app}.cpp )
    target_link_libraries(${app} MaCh3::Plotting MaCh3Warnings)

    install(TARGETS ${app} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach(app)
