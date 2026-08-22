set(HEADERS
    FitterBase.h
    MCMCBase.h
    MR2T2.h
    DelayedMR2T2.h
    LikelihoodFit.h
    MinuitFit.h
    PSO.h
    MaCh3Factory.h
    PredictiveThrower.h
    MulticanonicalMCMCHandler.h
)

add_library(Fitters SHARED
    FitterBase.cpp
    MCMCBase.cpp
    MR2T2.cpp
    DelayedMR2T2.cpp
    LikelihoodFit.cpp
    $<$<BOOL:${MaCh3_MINUIT2_ENABLED}>:MinuitFit.cpp>
    PSO.cpp
    PredictiveThrower.cpp
    MulticanonicalMCMCHandler.cpp
    MaCh3Factory.cpp
)

set_target_properties(Fitters PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME Fitters)

target_link_libraries(Fitters PUBLIC MaCh3::Samples MaCh3::FitProcessing ROOT::ROOTDataFrame)
target_link_libraries(Fitters PRIVATE MaCh3Warnings)

target_include_directories(Fitters PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS Fitters
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Fitters)

add_library(MaCh3::Fitters ALIAS Fitters)
target_link_libraries(MaCh3 INTERFACE Fitters)
