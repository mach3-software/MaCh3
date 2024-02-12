set(HEADERS ThrowParms.h
)

add_library(ThrowParms SHARED
    ThrowParms.cpp
)

set_target_properties(ThrowParms PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME ThrowParms)


target_link_libraries(ThrowParms ROOT::ROOT MaCh3CompilerOptions)
target_include_directories(ThrowParms PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS ThrowParms
        EXPORT MaCh3-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/throwParms)

add_library(MaCh3::ThrowParms ALIAS ThrowParms)
