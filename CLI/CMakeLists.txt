add_subdirectory(API)
add_subdirectory(Modules)

add_subdirectory(Diagnostics)
add_subdirectory(Plotting)

add_executable(mach3
    mach3.cpp
    MaCh3Program.cpp
)

target_link_libraries(mach3 PRIVATE MaCh3::CLIModules MaCh3::All)
target_link_libraries(mach3 PRIVATE MaCh3Warnings)

# libdl is only available as a separate library on Unix/Linux systems.
if(UNIX AND NOT APPLE)
    target_link_libraries(mach3 PRIVATE dl)
endif()
install(TARGETS mach3 DESTINATION bin)
