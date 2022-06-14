if (NOT DEFINED ENV{ROOTSYS})
  cmessage (FATAL_ERROR "$ROOTSYS is not defined, please set up ROOT first.")
else()
  cmessage(STATUS "Using ROOT installed at $ENV{ROOTSYS}")
  set(CMAKE_ROOTSYS $ENV{ROOTSYS})
endif()

# Get cflags from ROOT
execute_process (COMMAND root-config
  --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS_RAW OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE ";" " " ROOT_CXX_FLAGS "${ROOT_CXX_FLAGS_RAW}")
# Get libdir from ROOT
execute_process (COMMAND root-config
  --libdir OUTPUT_VARIABLE ROOT_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get ldflags from ROOT
execute_process (COMMAND root-config
  --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get version from ROOT
execute_process (COMMAND root-config
  --version OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get features from ROOT
execute_process (COMMAND root-config
  --features OUTPUT_VARIABLE ROOT_FEATURES OUTPUT_STRIP_TRAILING_WHITESPACE)

LIST(APPEND CMAKE_LINK_DIRS ${ROOT_LIBDIR})
LIST(APPEND CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})
LIST(APPEND CMAKE_LIBS ${ROOT_LIBS})
LIST(APPEND CMAKE_EXE_LINKER_FLAGS ${ROOT_LIBS})

LIST(APPEND ROOT_LIBS
  Core
  RIO
  Net
  Hist
  Graf
  Graf3d
  Gpad
  Tree
  Rint
  Postscript
  Matrix
  Physics
  MathCore)

cmessage ( STATUS "[ROOT]: root-config --version: ${ROOT_VERSION} ")
cmessage ( STATUS "[ROOT]: root-config --cflags : ${ROOT_CXX_FLAGS} ")
cmessage ( STATUS "[ROOT]: root-config --libs   : ${ROOT_LD_FLAGS} ")
#cmessage ( STATUS "[ROOT]: Libs                 : ${ROOT_LIBS} ")


