add_library(CLIPlotModules SHARED
    LLHModule.cpp
    PostFitParamModule.cpp
    DiagMCMCPlotModule.cpp
    MatrixModule.cpp
    SigmaVariationModule.cpp
    PredictiveModule.cpp
    PlotModule.cpp
)

target_link_libraries(CLIPlotModules PRIVATE MaCh3::All)
target_link_libraries(CLIPlotModules PRIVATE MaCh3Warnings)

target_include_directories(CLIPlotModules PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>)

install(TARGETS CLIPlotModules LIBRARY DESTINATION lib)

add_library(MaCh3::CLIPlotModules ALIAS CLIPlotModules)
