set(Throwparms_implementation_files
	ThrowParms.cpp)

set(Throwparms_header_files
	ThrowParms.h)

add_library(ThrowParms SHARED ${Throwparms_implementation_files})
target_include_directories(ThrowParms PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(ThrowParms ROOT::ROOT)

set_target_properties(ThrowParms PROPERTIES 
  PUBLIC_HEADER "${Throwparms_header_files}"
  EXPORT_NAME ThrowParms)

install(TARGETS ThrowParms
    EXPORT MaCh3-target
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/throwParms)
