
add_subdirectory(plottingUtils)

foreach(app
    GetPostfitParamPlots
    PlotLLH
    MatrixPlotter
    )
    add_executable(${app} ${app}.cpp )
    target_link_libraries(${app} MaCh3::Plotting)

    install(TARGETS ${app} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach(app)
