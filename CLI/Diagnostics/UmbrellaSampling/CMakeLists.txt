add_custom_target(UmbreallaSamplingApps)

foreach(app
    UmbrellaPlotting
    UmbrellaSolver
    )
    add_executable(${app} ${app}.cpp)
    add_dependencies(UmbreallaSamplingApps ${app})
    target_link_libraries(${app} MaCh3::All MaCh3Warnings)

    install(TARGETS ${app} DESTINATION bin)
endforeach(app)
