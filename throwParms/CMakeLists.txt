set(Throwparms_implementation_files
	ThrowParms.cpp)

set(Throwparms_header_files
	ThrowParms.h)

add_library(throwParms SHARED ${Throwparms_implementation_files})

target_link_libraries(throwParms)

install(TARGETS throwParms DESTINATION lib)

