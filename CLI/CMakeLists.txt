# Fix build/platform consistency: the loader supports Linux/macOS, not Windows DLLs.
if(NOT CMAKE_SYSTEM_NAME MATCHES "^(Linux|Darwin)$")
    message(FATAL_ERROR "The MaCh3 CLI currently supports Linux and macOS only")
endif()

add_subdirectory(API)
add_subdirectory(Modules)

add_subdirectory(Diagnostics)
add_subdirectory(Plotting)

add_executable(mach3
    mach3.cpp
    MaCh3Program.cpp
)

target_link_libraries(mach3 PRIVATE MaCh3::CLIPlotModules MaCh3::CLIModules MaCh3::All)
target_link_libraries(mach3 PRIVATE MaCh3Warnings)

# Fix build/platform consistency: let CMake select the platform's dynamic-loader libraries.
target_link_libraries(mach3 PRIVATE ${CMAKE_DL_LIBS})
install(TARGETS mach3 DESTINATION bin)
