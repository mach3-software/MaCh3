add_subdirectory(Plot)

add_library(CLIModules SHARED
    ProcessMCMCModule.cpp
    DiagMCMCModule.cpp
    LLHMapModule.cpp
    DiagModule.cpp
    GetPenaltyTermModule.cpp
    RHatModule.cpp
    SmearModule.cpp
    ReweightModule.cpp
    CombineModule.cpp
    UmbrellaModule.cpp
    UmbrellaPlotModule.cpp
)

target_link_libraries(CLIModules PRIVATE MaCh3::All)
target_link_libraries(CLIModules PRIVATE MaCh3Warnings)

target_include_directories(CLIModules PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>)

install(TARGETS CLIModules LIBRARY DESTINATION lib)

add_library(MaCh3::CLIModules ALIAS CLIModules)
