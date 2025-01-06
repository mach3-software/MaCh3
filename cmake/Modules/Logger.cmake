######################### Logging options #########################
set(SPDLOG_VERSION 1.11.0)
CPMAddPackage(
    NAME spdlog
    VERSION ${SPDLOG_VERSION}
    GITHUB_REPOSITORY gabime/spdlog
    OPTIONS
      "SPDLOG_BUILD_PIC ON"
      "SPDLOG_MASTER_PROJECT ON"
)
if(NOT TARGET spdlog::spdlog)
  cmessage(FATAL_ERROR "MaCh3 Expected dependency target: spdlog::spdlog")
endif()

#From: https://github.com/gabime/spdlog/blob/62302019babd8fdf63c1c6dc4c9faa48d441a7f0/include/spdlog/spdlog.h#L284C1-L291C20
# define SPDLOG_ACTIVE_LEVEL to one of those (before including spdlog.h):
# SPDLOG_LEVEL_TRACE,
# SPDLOG_LEVEL_DEBUG,
# SPDLOG_LEVEL_INFO,
# SPDLOG_LEVEL_WARN,
# SPDLOG_LEVEL_ERROR,
# SPDLOG_LEVEL_CRITICAL,
# SPDLOG_LEVEL_OFF

if(NOT DEFINED LOG_LEVEL)
  if(MaCh3_DEBUG_ENABLED)
    set(LOG_LEVEL "DEBUG")
  else()
    set(LOG_LEVEL "INFO")
  endif()
endif()

## EM: Check the specified log level is valid
set(VALID_LOG_OPTIONS OFF CRITICAL ERROR WARN INFO DEBUG TRACE)
list(FIND VALID_LOG_OPTIONS ${LOG_LEVEL} index)
if(${index} GREATER -1)
  cmessage(STATUS "LOG LEVEL: ${LOG_LEVEL}")
  target_compile_definitions(MaCh3CompilerOptions INTERFACE SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_${LOG_LEVEL})
else()
  cmessage(FATAL_ERROR "Invalid log level specified: ${LOG_LEVEL} \n Should be one of: ${VALID_LOG_OPTIONS}")
endif()

# KS: If logger is set to off many functions which sole purpose is to print message will not print anything. Thus variable will not be used and our very picky WErrror will throw errors
if(LOG_LEVEL STREQUAL "OFF" OR LOG_LEVEL STREQUAL "CRITICAL" OR LOG_LEVEL STREQUAL "ERROR" OR LOG_LEVEL STREQUAL "WARN")
  target_compile_options(MaCh3Warnings INTERFACE
    -Wno-error=unused-variable
    -Wno-error=unused-parameter
  )
endif()
