# set(HEADERS
#     ProcessMCMCPlugin.hpp
# )


add_library(CLIModules SHARED
    ProcessMCMCModule.cpp
    DiagMCMCModule.cpp
    GetPenaltyTermModule.cpp
)

target_link_libraries(CLIModules PRIVATE MaCh3::All)
target_include_directories(CLIModules PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>)

install(TARGETS CLIModules LIBRARY DESTINATION lib)

# set_target_properties(Plugins PROPERTIES
#     PUBLIC_HEADER "${HEADERS}"
#     EXPORT_NAME Plugins)


# install(TARGETS Plugins
#         EXPORT MaCh3-targets
#         LIBRARY DESTINATION lib/
#         PUBLIC_HEADER DESTINATION include/Plugins)

add_library(MaCh3::CLIModules ALIAS CLIModules)