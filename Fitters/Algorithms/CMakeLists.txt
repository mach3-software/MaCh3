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
    Hamiltonian/HamiltonianMC.cpp
    Hamiltonian/HamiltonianMetric.cpp
    Hamiltonian/NoUTurnMC.cpp
    Hybrid/HybridMCMC.cpp
    Hybrid/HybridSettings.cpp
    Hybrid/HybridBlocks.cpp
    Hybrid/HybridRandomWalk.cpp
    Hybrid/HybridDelayedAcceptance.cpp
    Hybrid/HybridFlipMoves.cpp
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
        PUBLIC_HEADER DESTINATION include/Fitters/Algorithms)
#HW The gradient samplers, one folder each. PUBLIC_HEADER would flatten them, so they are installed
#with their folders kept, which is how they are included
install(DIRECTORY Hamiltonian/ DESTINATION include/Fitters/Algorithms/Hamiltonian FILES_MATCHING PATTERN "*.h")
install(DIRECTORY Hybrid/ DESTINATION include/Fitters/Algorithms/Hybrid FILES_MATCHING PATTERN "*.h")

add_library(MaCh3::Fitters ALIAS Fitters)
target_link_libraries(MaCh3 INTERFACE Fitters)
