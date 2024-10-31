set(HEADERS
    Binning.h
)

add_library(Binning SHARED
    Binning.cpp
)

set_target_properties(Binning PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Binning)

target_link_libraries(Binning PUBLIC MaCh3CompilerOptions ROOT::Hist yaml-cpp spdlog)
target_link_libraries(Binning PRIVATE MaCh3Warnings)

target_include_directories(Binning PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS Binning
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Binning)

add_library(MaCh3::Binning ALIAS Binning)


add_executable(BinningTest BinningTest.cpp)
target_link_libraries(BinningTest PRIVATE MaCh3::Binning)