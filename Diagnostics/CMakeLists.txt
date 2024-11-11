add_custom_target(DiagApps)

foreach(app
    DiagMCMC
    RHat
    GetPenaltyTerm
    ProcessMCMC
    CombineMaCh3Chains
    )
    add_executable(${app} ${app}.cpp)
    add_dependencies(DiagApps ${app})
    target_link_libraries(${app} MaCh3::All MaCh3Warnings)

    install(TARGETS ${app} DESTINATION bin)
endforeach(app)

install(DIRECTORY Diagnostics_utils DESTINATION Diagnostics USE_SOURCE_PERMISSIONS)
