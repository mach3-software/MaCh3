add_library(CLIPlotModules SHARED
    LLHMapModule.cpp
    LLHModule.cpp
    PlotModule.cpp
)

target_link_libraries(CLIPlotModules PRIVATE MaCh3::All)
target_link_libraries(CLIPlotModules PRIVATE MaCh3Warnings)

target_include_directories(CLIPlotModules PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>)

install(TARGETS CLIPlotModules LIBRARY DESTINATION lib)

# set_target_properties(Plugins PROPERTIES
#     PUBLIC_HEADER "${HEADERS}"
#     EXPORT_NAME Plugins)


# install(TARGETS Plugins
#         EXPORT MaCh3-targets
#         LIBRARY DESTINATION lib/
#         PUBLIC_HEADER DESTINATION include/Plugins)

add_library(MaCh3::CLIPlotModules ALIAS CLIPlotModules)
