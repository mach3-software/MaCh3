add_custom_target(DiagApps)

foreach(app
    DiagMCMC
    RHat
    GetPenaltyTerm
    )
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} MaCh3::All)

    add_dependencies(DiagApps ${app})

endforeach(app)

install(FILES
    Diag_Config.yaml
    DESTINATION ${CMAKE_BINARY_DIR}/Diagnostics)

install(DIRECTORY Diagnostics_utils DESTINATION ${CMAKE_BINARY_DIR}/Diagnostics)
