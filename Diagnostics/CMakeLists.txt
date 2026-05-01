# set(HEADERS
#     ProcessMCMCPlugin.hpp
# )


add_library(CorePlugins STATIC
    ProcessMCMCPlugin.cpp
    DiagMCMCPlugin.cpp
    GetPenaltyTermPlugin.cpp
)

target_link_libraries(CorePlugins PRIVATE MaCh3::All)
target_include_directories(CorePlugins PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>)

# set_target_properties(Plugins PROPERTIES
#     PUBLIC_HEADER "${HEADERS}"
#     EXPORT_NAME Plugins)


# install(TARGETS Plugins
#         EXPORT MaCh3-targets
#         LIBRARY DESTINATION lib/
#         PUBLIC_HEADER DESTINATION include/Plugins)

add_library(MaCh3::CorePlugins ALIAS CorePlugins)

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
    target_link_libraries(${app} MaCh3::CorePlugins MaCh3::All MaCh3Warnings ROOT::ROOTDataFrame)  # dont need plugins in MaCh3::All at root scope as only used here. Once all transferred to plugins wont need MaCh3::All either

    install(TARGETS ${app} DESTINATION bin)
endforeach(app)

add_subdirectory(UmbrellaSampling)
