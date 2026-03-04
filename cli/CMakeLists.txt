    add_executable(mach3 mach3.cpp)
    target_link_libraries(mach3 MaCh3::All MaCh3::Plugins)

    install(TARGETS mach3 DESTINATION bin)
