add_subdirectory(modules)
# set(CLI_MODULES
#     ProcessMCMCPlugin.cpp
#     DiagMCMCPlugin.cpp
#     GetPenaltyTermPlugin.cpp
# )

# list(TRANSFORM CLI_MODULES PREPEND "modules/")



add_executable(mach3
    mach3.cpp
    MaCh3Program.cpp
    # ${CLI_MODULES}
)

target_link_libraries(mach3 PRIVATE MaCh3::CLIModules MaCh3::All)

install(TARGETS mach3 DESTINATION bin)
