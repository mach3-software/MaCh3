    add_executable(mach3 mach3.cpp MaCh3Program.cpp)
    target_link_libraries(mach3 PRIVATE MaCh3::CorePlugins MaCh3::All)

    install(TARGETS mach3 DESTINATION bin)
