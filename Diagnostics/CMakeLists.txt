add_custom_target(DiagApps)

foreach(app
    DiagMCMC
    RHat
    RHat_HighMem
    GetPenaltyTerm
    ProcessMCMC
    CombineMaCh3Chains
    PlotMCMCDiag
    SmearChain
    ReweightMCMC
    )
    add_executable(${app} ${app}.cpp)
    add_dependencies(DiagApps ${app})
    target_link_libraries(${app} MaCh3::All MaCh3Warnings)
    set_target_properties(${app} PROPERTIES
        BUILD_RPATH "${CMAKE_BINARY_DIR}/lib"
        INSTALL_RPATH "@loader_path/../lib"
    )

    install(TARGETS ${app} DESTINATION bin)
endforeach(app)
