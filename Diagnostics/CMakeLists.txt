set(HEADERS
    ProcessMCMCPlugin.hpp
)


add_library(Plugins STATIC
    ProcessMCMCPlugin.cpp
)

target_link_libraries(Plugins PRIVATE MaCh3::All)
target_include_directories(Plugins PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>)

# set_target_properties(Plugins PROPERTIES
#     PUBLIC_HEADER "${HEADERS}"
#     EXPORT_NAME Plugins)


# install(TARGETS Plugins
#         EXPORT MaCh3-targets
#         LIBRARY DESTINATION lib/
#         PUBLIC_HEADER DESTINATION include/Plugins)

add_library(MaCh3::Plugins ALIAS Plugins)

add_custom_target(DiagApps)

foreach(app
    DiagMCMC
    RHat
    RHat_HighMem
    GetPenaltyTerm
    ProcessMCMC
    CombineMaCh3Chains
    PlotMCMCDiag
    SmearChain
    ReweightMCMC
    PlotLLHMap
    )
    add_executable(${app} ${app}.cpp)
    add_dependencies(DiagApps ${app})
    target_link_libraries(${app} MaCh3::All MaCh3::Plugins MaCh3Warnings ROOT::ROOTDataFrame)  # dont need plugins in MaCh3::All at root scope as only used here. Once all transferred to plugins wont need MaCh3::All either

    install(TARGETS ${app} DESTINATION bin)
endforeach(app)

add_subdirectory(UmbrellaSampling)
