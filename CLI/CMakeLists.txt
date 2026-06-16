add_subdirectory(API)
add_subdirectory(Modules)


add_executable(mach3
    mach3.cpp
    MaCh3Program.cpp
)

target_link_libraries(mach3 PRIVATE MaCh3::CLIModules MaCh3::All)

install(TARGETS mach3 DESTINATION bin)
