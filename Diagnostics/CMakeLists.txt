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
    )
    add_executable(${app} ${app}.cpp)
    add_dependencies(DiagApps ${app})
    target_link_libraries(${app} MaCh3::All MaCh3Warnings)

    install(TARGETS ${app} DESTINATION bin)
endforeach(app)
